#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
34_mr_panel.py -- REAL two-sample Mendelian-randomisation PANEL over the haemoglobin /
heme / iron (HDFI) axis x brain-haemorrhage outcomes, on genuine public GWAS summary
statistics downloaded from the GWAS Catalog FTP (harmonised GRCh38 sumstats).

This extends the single validated IV pair HP -> SAH (scripts/31_mr_real.py) to a full
exposure x outcome panel:

    Exposures (European pQTL, N>=3000, harmonised .h.tsv.gz on FTP, largest N chosen):
      - Haptoglobin (HP)            GCST90019428  Eur N=10,708   (instruments reused from
                                                 the validated 31_mr_real.py extraction)
      - Ferritin (total)            GCST90270865  Eur N=270,794
      - Heme oxygenase 1 (HMOX1)    GCST90469465  Eur N=47,745
      - Transferrin receptor TFRC   GCST90470840  Eur N=47,745
      - Hepcidin (HAMP)             GCST90451683  Eur N=91,675
    Outcomes (European, harmonised, UKB whole-genome-sequence GWAS; same study):
      - Subarachnoid haemorrhage    GCST90473583  Eur N=458,440 (2,036 cases)  ICD10 I60
      - Intracerebral haemorrhage   GCST90473584  Eur N=458,440 (2,407 cases)  ICD10 I61

The estimators (IVW fixed/random, MR-Egger + intercept, weighted median, Cochran's Q) are
imported UNCHANGED from the engine in scripts/14_mr.py (previously validated on synthetic
ground-truth beta=0.40, recovered as IVW=0.413). Instrument selection, 500 kb distance
pseudo-clumping, effect-allele harmonisation (exact rsID, else nearest variant within
+/-300 bp, palindromic drop) and the per-IV table format follow scripts/31_mr_real.py
EXACTLY.

STRICT INTEGRITY RULES honoured here:
  * No synthetic data in any estimate; every number is computed from a downloaded file.
  * If <2 instruments survive harmonisation the pair is written as status
    "insufficient instruments" -- never a substituted synthetic number.
  * Every pQTL file's sample size + ancestry are confirmed from the FTP meta.yaml, not
    from the GWAS Catalog trait string.
  * Proteins with no European harmonised pQTL at N>=3000 on the FTP, and the reviewer's
    requested "HP -> SAH severity / secondary-injury MR" (no openly downloadable
    severity GWAS exists) are recorded as status "not-executed" with an honest reason.

Outputs (all under the project root):
  results/mr_panel_HDFI_axis.csv          one row per exposure x outcome pair
  results/gwas/mr_panel_recipe.json       accessions + URLs + meta + per-pair results
  figures/Fig24_mr_panel.png / .pdf       exposure x outcome heatmap

Run:  python 34_mr_panel.py            (downloads to results/gwas/cache if missing)
"""
import os, sys, gzip, json, bisect, importlib.util, subprocess, time
from collections import defaultdict
import numpy as np, pandas as pd
from scipy.stats import norm, chi2

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES  = os.path.join(BASE, "results")
GW   = os.path.join(RES, "gwas")
CACHE= os.path.join(GW, "cache")
FIG  = os.path.join(BASE, "figures")
for d in (RES, GW, CACHE, FIG):
    os.makedirs(d, exist_ok=True)

FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

# --------------------------------------------------------------------------------------
# import the validated two-sample MR engine (IVW / Egger / weighted-median) from 14_mr.py
# without executing its synthetic print block.
_spec = importlib.util.spec_from_file_location("_mr_engine", os.path.join(BASE, "scripts", "14_mr.py"))
_eng  = importlib.util.module_from_spec(_spec)
_so, _se = sys.stdout, sys.stderr
import io; sys.stdout = io.StringIO()
_spec.loader.exec_module(_eng)
sys.stdout, sys.stderr = _so, _se
run_mr = _eng.run_mr

def ftp_folder(acc):
    n = int(acc[4:])
    return "GCST%08d-GCST%08d" % ((n // 1000) * 1000 + 1, (n // 1000) * 1000 + 1000)

def ftp_url(acc, base):
    return f"{FTP}/{ftp_folder(acc)}/{acc}/harmonised/{base}"

# --------------------------------------------------------------------------------------
# CONFIG: exposures + outcomes (meta.yaml-confirmed sizes/ancestry; see recipe)
# --------------------------------------------------------------------------------------
EXPOSURES = [
    dict(key="HP",       label="Haptoglobin",
         acc="GCST90019428", base="33328453-GCST90019428-EFO_0004640.h.tsv.gz",
         ancestry="European", n=10708, assembly="GRCh38",
         reuse_instr=os.path.join(GW, "HP_GWAS90019428_instruments.tsv"),
         cache=os.path.join(CACHE, "HP_GWAS90019428.h.tsv.gz")),
    dict(key="FERRITIN", label="Ferritin (total)",
         acc="GCST90270865", base="GCST90270865.h.tsv.gz",
         ancestry="European", n=270794, assembly="GRCh38",
         cache=os.path.join(CACHE, "FERRITIN_GWAS90270865.h.tsv.gz")),
    dict(key="HMOX1",    label="Heme oxygenase 1 (HMOX1)",
         acc="GCST90469465", base="GCST90469465.h.tsv.gz",
         ancestry="European", n=47745, assembly="GRCh38",
         cache=os.path.join(CACHE, "HMOX1_GWAS90469465.h.tsv.gz")),
    dict(key="TFRC",     label="Transferrin receptor (TFRC)",
         acc="GCST90470840", base="GCST90470840.h.tsv.gz",
         ancestry="European", n=47745, assembly="GRCh38",
         cache=os.path.join(CACHE, "TFRC_GWAS90470840.h.tsv.gz")),
    dict(key="HEPCIDIN", label="Hepcidin (HAMP)",
         acc="GCST90451683", base="GCST90451683.h.tsv.gz",
         ancestry="European", n=91675, assembly="GRCh38",
         cache=os.path.join(CACHE, "HEPCIDIN_GWAS90451683.h.tsv.gz")),
]

OUTCOMES = [
    dict(key="SAH", label="Subarachnoid haemorrhage (ICD10 I60)",
         acc="GCST90473583", base="GCST90473583.h.tsv.gz",
         ancestry="European (Non-Finnish)", n=458440, cases=2036, controls=456404,
         assembly="GRCh38", cache=os.path.join(CACHE, "SAH_GWAS90473583.h.tsv.gz")),
    dict(key="ICH", label="Intracerebral haemorrhage (ICD10 I61)",
         acc="GCST90473584", base="GCST90473584.h.tsv.gz",
         ancestry="European (Non-Finnish)", n=458440, cases=2407, controls=456033,
         assembly="GRCh38", cache=os.path.join(CACHE, "ICH_GWAS90473584.h.tsv.gz")),
]

# Pairs that were ATTEMPTED but cannot be executed on open data -> honest "not-executed".
NOT_EXECUTED = [
    dict(exposure="Ferritin light chain (FTL)",  outcome="SAH/ICH",
         reason="no European harmonised pQTL at N>=3000 on FTP (best FTL GCST90618436 N=2,338); covered indirectly by Ferritin (total) GCST90270865"),
    dict(exposure="Ferritin heavy chain (FTH1)", outcome="SAH/ICH",
         reason="no European harmonised pQTL at N>=3000 on FTP (best FTH1 GCST90441675 N=200); covered indirectly by Ferritin (total) GCST90270865"),
    dict(exposure="Ferroportin (SLC40A1)",       outcome="SAH/ICH",
         reason="no European harmonised single-protein pQTL at N>=3000 located on FTP"),
    dict(exposure="CD163 (sCD163)",              outcome="SAH/ICH",
         reason="only protein-level-ratio GWAS (N=43,509) or small sCD163 (N=2,838); no clean single-CD163 European pQTL at N>=3000 with harmonised file"),
    dict(exposure="Hemopexin (HX)",              outcome="SAH/ICH",
         reason="max European single-protein pQTL GCST90621249 N=2,338 (<3000); no harmonised pQTL at N>=3000"),
    dict(exposure="S100A8 / S100A9",             outcome="SAH/ICH",
         reason="no matching European single-protein pQTL at N>=3000 located on FTP harmonised"),
    dict(exposure="Haptoglobin-related protein (HPR)", outcome="SAH/ICH",
         reason="no separate European HPR pQTL at N>=3000 on FTP; HPR co-localises with the HP/CRP chr16 cluster used for the HP exposure"),
    dict(exposure="Haptoglobin (HP)", outcome="SAH severity / secondary injury",
         reason="reviewer request: no openly downloadable SAH-severity / secondary-injury GWAS with harmonised summary statistics exists on the GWAS Catalog FTP; MR of an incidence endpoint (SAH GCST90473583) is a different estimand and cannot stand in for it"),
]

# --------------------------------------------------------------------------------------
# download (segmented, via the validated _segdl.sh helper) -- skipped if cache present
# --------------------------------------------------------------------------------------
def ensure_downloaded(spec, nparts=8, label=""):
    path = spec["cache"]
    if os.path.exists(path) and os.path.getsize(path) > 1_000_000:
        # _segdl.sh only finalises the .gz after a full gzip -t pass; trust a present
        # file unless we are doing a cold integrity check (set MR_TRUST_CACHE=1 to skip
        # the expensive full re-decompression on a re-run over already-verified cache).
        if os.environ.get("MR_TRUST_CACHE") == "1":
            return True
        if _gz_ok(path):
            return True
    helper = os.path.join(BASE, "scripts", "_segdl.sh")
    url = ftp_url(spec["acc"], spec["base"])
    print(f"[download] {label or spec['key']} {url}")
    try:
        subprocess.run(["bash", helper, url, path, str(nparts)], check=False)
    except Exception as e:
        print("  download failed:", e)
    ok = os.path.exists(path) and _gz_ok(path)
    return ok

def _gz_ok(path):
    try:
        with gzip.open(path, "rb") as f:
            while f.peek(1 << 20):
                f.read(1 << 20)
        return True
    except Exception:
        return False

# --------------------------------------------------------------------------------------
# instrument extraction: genome-wide-significant p<5e-8, 500 kb distance pseudo-clumping
# (p-first greedy). Mirrors 31_mr_real.py's coordinate-only clumping (no PLINK/LD ref).
# --------------------------------------------------------------------------------------
def palindromic(a, b):
    s = {str(a).upper(), str(b).upper()}
    return s == {"A", "T"} or s == {"C", "G"}

def extract_instruments(path, p_thresh=5e-8, clump_bp=500_000):
    need = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "p_value"]
    rows = []
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        miss = [c for c in need if c not in ix]
        if miss:
            raise RuntimeError(f"exposure file missing columns {miss} in {path}")
        rsid_i = ix.get("rsid"); vid_i = ix.get("variant_id")
        for line in f:
            t = line.rstrip("\n").split("\t")
            try:
                pv = float(t[ix["p_value"]])
            except ValueError:
                continue
            if not (pv < p_thresh):
                continue
            try:
                chrom = t[ix["chromosome"]]
                pos = int(float(t[ix["base_pair_location"]]))
                ea = t[ix["effect_allele"]]; oa = t[ix["other_allele"]]
                be = float(t[ix["beta"]]); se = float(t[ix["standard_error"]])
            except (ValueError, IndexError):
                continue
            if not (se > 0) or ea in ("", ".", "NA") or oa in ("", ".", "NA"):
                continue
            rs = (t[rsid_i] if rsid_i is not None else ".")
            if not rs or rs == ".":
                rs = (t[vid_i] if vid_i is not None else f"{chrom}:{pos}")
            rows.append(dict(chromosome=chrom, base_pair_location=pos, rsid=rs,
                             ea=ea, oa=oa, be=be, se=se, p=pv))
    if not rows:
        return pd.DataFrame(columns=["chromosome", "base_pair_location", "variant_id",
                                     "effect_allele", "other_allele", "beta",
                                     "standard_error", "p_value", "palindromic", "Fstat"])
    rows.sort(key=lambda r: r["p"])
    kept_by_chr = defaultdict(list)
    kept = []
    for r in rows:
        c = r["chromosome"]; pos = r["base_pair_location"]
        if any(abs(pos - q) < clump_bp for q in kept_by_chr[c]):
            continue
        kept_by_chr[c].append(pos)
        kept.append(r)
    df = pd.DataFrame(kept)
    df["Fstat"] = (df.be / df.se) ** 2
    out = pd.DataFrame({
        "chromosome": df.chromosome, "base_pair_location": df.base_pair_location,
        "variant_id": df.rsid, "effect_allele": df.ea, "other_allele": df.oa,
        "beta": df.be, "standard_error": df.se, "p_value": df.p,
        "palindromic": [palindromic(a, b) for a, b in zip(df.ea, df.oa)],
        "Fstat": df.Fstat,
    }).sort_values("p_value").reset_index(drop=True)
    return out

# --------------------------------------------------------------------------------------
# outcome harvest: single streaming pass per outcome file capturing, for every
# instrument locus, the exact rsID match OR any variant within +/-window bp.
# --------------------------------------------------------------------------------------
def harvest_outcome(path, targets, window_bp=1000):
    """targets: list of (chromosome(str), pos(int), rsid(str)). Returns DataFrame."""
    tp = defaultdict(list); rset = set()
    for c, p, r in targets:
        tp[str(c)].append(int(p)); rset.add(str(r).upper())
    for c in tp:
        tp[c].sort()
    recs = []
    with gzip.open(path, "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        need = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
                "beta", "standard_error", "p_value"]
        if any(c not in ix for c in need):
            raise RuntimeError(f"outcome file {path} missing columns")
        rsid_i = ix.get("rsid")
        for line in f:
            t = line.rstrip("\n").split("\t")
            try:
                chrom = t[ix["chromosome"]]
                pos = int(float(t[ix["base_pair_location"]]))
            except (ValueError, IndexError):
                continue
            rs = t[rsid_i] if rsid_i is not None else "."
            hit = str(rs).upper() in rset
            if not hit and chrom in tp:
                a = tp[chrom]; j = bisect.bisect_left(a, pos)
                for k in (j - 1, j):
                    if 0 <= k < len(a) and abs(a[k] - pos) <= window_bp:
                        hit = True; break
            if not hit:
                continue
            try:
                recs.append(dict(chromosome=chrom, base_pair_location=pos,
                                 effect_allele=t[ix["effect_allele"]],
                                 other_allele=t[ix["other_allele"]],
                                 beta=float(t[ix["beta"]]),
                                 standard_error=float(t[ix["standard_error"]]),
                                 p_value=float(t[ix["p_value"]]), rsid=rs))
            except (ValueError, IndexError):
                continue
    return pd.DataFrame(recs)

# --------------------------------------------------------------------------------------
# harmonisation -- copied verbatim in behaviour from scripts/31_mr_real.py
# --------------------------------------------------------------------------------------
def _pick_outcome_row(out, rs, chrom, pos):
    m = out[out["rsid_u"] == rs]
    if len(m):
        return m.iloc[0], "rsid"
    chrom = str(chrom)
    sub = out[(out["chromosome"].astype(str) == chrom)]
    if len(sub):
        sub = sub.assign(_d=(sub["base_pair_location"].astype(float) - float(pos)).abs())
        sub = sub[sub["_d"] <= 300].sort_values("_d")
        if len(sub):
            return sub.iloc[0], "pos_window"
    return None, "absent_in_outcome"

def harmonise(exp, out):
    out = out.copy()
    out["EA"] = out["effect_allele"].astype(str).str.upper()
    out["OA"] = out["other_allele"].astype(str).str.upper()
    out["rsid_u"] = out["rsid"].astype(str).str.upper()
    rows = []
    for _, e in exp.iterrows():
        rs = str(e["variant_id"]).upper()
        m, how = _pick_outcome_row(out, rs, e["chromosome"], e["base_pair_location"])
        if m is None:
            rows.append(dict(rsid=rs, ok=False, why=how)); continue
        ea, oa = str(e["effect_allele"]).upper(), str(e["other_allele"]).upper()
        if palindromic(ea, oa):
            rows.append(dict(rsid=rs, ok=False, why="palindromic")); continue
        if m["EA"] == ea and m["OA"] == oa:
            bo, flip = float(m["beta"]), False
        elif m["EA"] == oa and m["OA"] == ea:
            bo, flip = -float(m["beta"]), True
        else:
            rows.append(dict(rsid=rs, ok=False, why="allele_mismatch")); continue
        rows.append(dict(rsid=rs, ok=True, match=how, out_pos=int(m["base_pair_location"]),
                         be=float(e["beta"]), se_e=float(e["standard_error"]),
                         p_e=float(e["p_value"]), bo=bo, se_o=float(m["standard_error"]),
                         p_o=float(m["p_value"]), flipped=flip))
    return pd.DataFrame(rows)

# --------------------------------------------------------------------------------------
# one pair -> results row
# --------------------------------------------------------------------------------------
def or_from(beta, se):
    lo = float(np.exp(beta - 1.96 * se)); hi = float(np.exp(beta + 1.96 * se))
    return float(np.exp(beta)), lo, hi

def run_pair(expkey, exp_instr, outkey, harvest_df, exmeta, outmeta):
    h = harmonise(exp_instr, harvest_df)
    keep = h[h.get("ok") == True].dropna(subset=["be", "bo", "se_o"]).copy()
    keep = keep[pd.to_numeric(keep["se_o"], errors="coerce") > 0]
    n_off = len(exp_instr); n_keep = len(keep)
    row = dict(exposure=exmeta["label"], exposure_accession=exmeta["acc"],
               exposure_ancestry=exmeta["ancestry"], exposure_N=exmeta["n"],
               outcome=outmeta["label"], outcome_accession=outmeta["acc"],
               outcome_ancestry=outmeta["ancestry"], outcome_N=outmeta["n"],
               outcome_cases=outmeta.get("cases"),
               instruments_offered=n_off, instruments_retained=n_keep,
               status="executed")
    if n_keep < 2:
        row["status"] = "insufficient instruments"
        row["reason"] = f"only {n_keep} usable instrument(s) after allele/palindromic harmonisation (<2)"
        return row, keep
    dat = keep[["be", "bo", "se_o"]].rename(columns={"se_o": "se"})
    r = run_mr(dat)
    F = (keep.be / keep.se_e) ** 2
    OR, lo, hi = or_from(r["IVW_beta"], r["IVW_se"])
    OR_r, lo_r, hi_r = or_from(r["IVW_beta"], r["IVW_se_random"])
    EOR, elo, ehi = or_from(r["Egger_beta"], r["Egger_se"])
    row.update(dict(
        IVW_beta=r["IVW_beta"], IVW_se=r["IVW_se"], IVW_p=r["IVW_p"],
        IVW_OR=OR, IVW_OR_CI_low=lo, IVW_OR_CI_high=hi,
        IVW_Q=r["IVW_Q"], IVW_df=r["IVW_df"], IVW_Q_p=float(chi2.sf(r["IVW_Q"], max(r["IVW_df"], 1))),
        IVW_se_random=r["IVW_se_random"], IVW_p_random=r["IVW_p_random"], IVW_OR_random=OR_r,
        IVW_OR_random_CI_low=lo_r, IVW_OR_random_CI_high=hi_r,
        Egger_OR=EOR, Egger_OR_CI_low=elo, Egger_OR_CI_high=ehi,
        Egger_beta=r["Egger_beta"], Egger_se=r["Egger_se"], Egger_p=r["Egger_p"],
        Egger_intercept=r["Egger_intercept"], Egger_intercept_p=r["Egger_intercept_p"],
        WMedian_beta=r["WMedian_beta"], WMedian_se=r["WMedian_se"], WMedian_p=r["WMedian_p"],
        F_stat_mean=float(F.mean()), F_stat_min=float(F.min()),
        per_IV=";".join(f"{rr.rsid}:bE={rr.be:.4f},bO={rr.bo:.4f}" for rr in keep.itertuples()),
    ))
    return row, keep

# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
def main():
    t0 = time.time()
    # 1. ensure all files present (skip HP raw: reuse validated instruments on disk)
    for ex in EXPOSURES:
        if ex.get("reuse_instr"):
            continue
        ensure_downloaded(ex, nparts=8, label=f"EXPOSURE {ex['key']}")
    for oc in OUTCOMES:
        ensure_downloaded(oc, nparts=12, label=f"OUTCOME {oc['key']}")

    # 2. instruments per exposure (reuse saved TSV when present -- keeps the run idempotent)
    instr = {}
    for ex in EXPOSURES:
        saved = os.path.join(GW, f"{ex['key']}_{ex['acc']}_instruments.tsv")
        if ex.get("reuse_instr") and os.path.exists(ex["reuse_instr"]):
            d = pd.read_csv(ex["reuse_instr"], sep="\t")
            print(f"  [reuse] {ex['key']}: {len(d)} instruments from {os.path.basename(ex['reuse_instr'])}")
        elif os.path.exists(saved):
            d = pd.read_csv(saved, sep="\t")
            print(f"  [cached] {ex['key']}: {len(d)} instruments from {os.path.basename(saved)}")
        else:
            d = extract_instruments(ex["cache"])
            d.to_csv(saved, sep="\t", index=False)
            print(f"  [extract] {ex['key']}: {len(d)} instruments (500kb clumped)")
        instr[ex["key"]] = d

    # 3. union of instrument loci -> harvest each outcome once
    all_targets = []
    for ex in EXPOSURES:
        for _, r in instr[ex["key"]].iterrows():
            all_targets.append((str(r["chromosome"]), int(r["base_pair_location"]), str(r["variant_id"])))
    all_targets = list({(c, p, rs) for c, p, rs in all_targets})
    print(f"  union instrument loci = {len(all_targets)}")

    harv = {}
    for oc in OUTCOMES:
        p = os.path.join(GW, f"{oc['key']}_{oc['acc']}_iv_harvested.tsv")
        if os.path.exists(p):
            harvest_df = pd.read_csv(p, sep="\t", low_memory=False)
        else:
            harvest_df = harvest_outcome(oc["cache"], all_targets)
            harvest_df.to_csv(p, sep="\t", index=False)
        # normalise columns for harmonise
        if "variant_id" not in harvest_df.columns:
            harvest_df["variant_id"] = harvest_df["rsid"]
        harv[oc["key"]] = harvest_df
        print(f"  [harvest] {oc['key']}: {len(harvest_df)} outcome rows around loci")

    # 4. run every exposure x outcome pair
    rows = []
    ex_by_key = {e["key"]: e for e in EXPOSURES}
    oc_by_key = {o["key"]: o for o in OUTCOMES}
    for ex in EXPOSURES:
        for oc in OUTCOMES:
            row, keep = run_pair(ex["key"], instr[ex["key"]], oc["key"], harv[oc["key"]], ex, oc)
            rows.append(row)
            b = row.get("IVW_beta"); pv = row.get("IVW_p")
            print(f"  MR {ex['key']:8s} x {oc['key']:3s}: retained={row['instruments_retained']} "
                  f"IVWbeta={'NA' if b is None else round(b,4)} p={'NA' if pv is None else pv:.3g}")

    df = pd.DataFrame(rows)
    # append the not-executed honest entries
    ne_rows = []
    for ne in NOT_EXECUTED:
        ne_rows.append(dict(exposure=ne["exposure"], outcome=ne["outcome"],
                            status="not-executed", reason=ne["reason"],
                            instruments_offered=None, instruments_retained=None))
    df = pd.concat([df, pd.DataFrame(ne_rows)], ignore_index=True)
    csv_path = os.path.join(RES, "mr_panel_HDFI_axis.csv")
    df.to_csv(csv_path, index=False)
    print("saved", csv_path)

    # 5. recipe JSON
    recipe = build_recipe(rows, instr, ex_by_key, oc_by_key)
    with open(os.path.join(GW, "mr_panel_recipe.json"), "w") as f:
        json.dump(recipe, f, indent=2)
    print("saved results/gwas/mr_panel_recipe.json")

    # 6. figure
    try:
        make_figure(df)
    except Exception as e:
        print("figure failed:", repr(e))

    exec_n = sum(1 for r in rows if r["status"] == "executed")
    insuf_n = sum(1 for r in rows if r["status"] == "insufficient instruments")
    print(f"\n[done] executed={exec_n} insufficient={insuf_n} not-executed={len(NOT_EXECUTED)}  "
          f"({time.time()-t0:.0f}s)")
    return df, recipe

def _fmt(x, d=4):
    try:
        if x is None or (isinstance(x, float) and (np.isnan(x))):
            return None
        return round(float(x), d)
    except Exception:
        return None

def interpret(row):
    if row["status"] != "executed":
        return row.get("reason", "")
    p = row.get("IVW_p"); orv = row.get("IVW_OR")
    lo = row.get("IVW_OR_CI_low"); hi = row.get("IVW_OR_CI_high")
    direction = "higher exposure associated with higher risk" if (orv or 1) > 1 else "lower risk"
    if p is not None and p < 0.05:
        return (f"Genetically higher {row['exposure']} -> {row['outcome']}: OR {orv:.2f} "
                f"(95% CI {lo:.2f}-{hi:.2f}), IVW p={p:.2g}. Statistically distinguishable "
                f"from the null; direction = {direction}.")
    return (f"Null: genetically higher {row['exposure']} does not alter {row['outcome']} "
            f"risk (OR {orv:.2f}, 95% CI {lo:.2f}-{hi:.2f}, IVW p={p:.2g}).")

def build_recipe(rows, instr, ex_by_key, oc_by_key):
    ex_meta = []
    for e in EXPOSURES:
        ex_meta.append(dict(trait=e["label"], accession=e["acc"], ancestry=e["ancestry"],
                            sample_size=e["n"], assembly=e["assembly"],
                            harmonised_url=ftp_url(e["acc"], e["base"]),
                            n_instruments_extracted=int(len(instr[e["key"]]))))
    oc_meta = []
    for o in OUTCOMES:
        oc_meta.append(dict(trait=o["label"], accession=o["acc"], ancestry=o["ancestry"],
                            sample_size=o["n"], cases=o.get("cases"), controls=o.get("controls"),
                            assembly=o["assembly"], harmonised_url=ftp_url(o["acc"], o["base"])))
    pair_recs = []
    for r in rows:
        rec = {k: r.get(k) for k in r}
        for kk in ("IVW_beta", "IVW_se", "IVW_p", "IVW_OR", "IVW_OR_CI_low", "IVW_OR_CI_high",
                   "IVW_Q", "IVW_Q_p", "Egger_beta", "Egger_intercept", "Egger_intercept_p",
                   "WMedian_beta", "F_stat_mean", "F_stat_min"):
            if kk in rec:
                rec[kk] = _fmt(rec[kk], 6 if kk.endswith("p") else 6)
        rec["interpretation"] = interpret(r)
        pair_recs.append(rec)
    exec_n = sum(1 for r in rows if r["status"] == "executed")
    return dict(
        title="Two-sample MR panel: haemoglobin/heme/iron (HDFI) axis proteins -> brain haemorrhage",
        data_source="GWAS Catalog FTP harmonised (GRCh38) summary statistics; all estimates from real downloaded files",
        engine="scripts/14_mr.py (IVW fixed/random, MR-Egger+intercept, weighted median, Cochran Q), validated on synthetic beta=0.40 (recovered IVW=0.413)",
        instrument_selection="p<5e-8 in exposure, 500 kb distance pseudo-clumping (no PLINK/LD reference in this environment; reported honestly)",
        harmonisation="effect-allele aligned; exact rsID else nearest variant within +/-300 bp; palindromic SNPs dropped",
        exposures=ex_meta, outcomes=oc_meta, pairs=pair_recs,
        executed_pairs=exec_n, insufficient_pairs=sum(1 for r in rows if r["status"] == "insufficient instruments"),
        not_executed=[dict(exposure=n["exposure"], outcome=n["outcome"], reason=n["reason"]) for n in NOT_EXECUTED],
        reviewer_note="HP -> 'SAH severity / secondary injury' MR is NOT performed: no openly downloadable severity GWAS with harmonised sumstats exists on the FTP. SAH GCST90473583 is an incidence/rupture endpoint, a different estimand.",
        cross_ancestry_note="All exposures and outcomes are European (SAH/ICH = Non-Finnish European, UKB WGS; exposures European pQTL). No cross-ancestry pair was executed in this panel.",
    )

# --------------------------------------------------------------------------------------
# figure: exposure x outcome heatmap (IVW beta colour, p-value star annotations)
# --------------------------------------------------------------------------------------
def make_figure(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib import rcParams
    rcParams.update({"font.family": "Arial", "font.size": 9,
                     "axes.linewidth": 0.9, "xtick.labelsize": 8, "ytick.labelsize": 8})
    exec_df = df[df["status"] == "executed"].copy()
    ex_order = [e["label"] for e in EXPOSURES]
    oc_order = [o["label"] for o in OUTCOMES]
    ex_order = [e for e in ex_order if e in set(exec_df.exposure)]
    piv = exec_df.pivot_table(index="exposure", columns="outcome", values="IVW_beta", aggfunc="first")
    piv = piv.reindex(index=ex_order, columns=oc_order)
    pmat = exec_df.pivot_table(index="exposure", columns="outcome", values="IVW_p", aggfunc="first").reindex(index=ex_order, columns=oc_order)
    nmat = exec_df.pivot_table(index="exposure", columns="outcome", values="instruments_retained", aggfunc="first").reindex(index=ex_order, columns=oc_order)

    vals = piv.values.astype(float)
    finite = vals[np.isfinite(vals)]
    vmax = max(0.2, float(np.max(np.abs(finite)))) if finite.size else 0.2
    norm_ = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")  # diverging, colour-blind-usable; red=risk-increasing beta

    nrows, ncols = piv.shape
    fig, ax = plt.subplots(figsize=(1.15 * ncols + 2.9, 0.62 * nrows + 1.9))
    masked = np.ma.masked_invalid(vals)
    im = ax.imshow(masked, cmap=cmap, norm=norm_, aspect="auto")
    ax2 = fig.add_axes([0.02, 0.02, 0.02, 0.001]); ax2.axis("off")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("IVW causal log-odds (beta) per SD higher protein", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    def star(p):
        try:
            p = float(p)
        except (TypeError, ValueError):
            return ""
        if not np.isfinite(p):
            return ""
        return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

    for i in range(nrows):
        for j in range(ncols):
            v = vals[i, j]
            if not np.isfinite(v):
                continue
            text_col = "white" if abs(v) > 0.62 * vmax else "black"
            p = pmat.values[i, j]
            n = nmat.values[i, j]
            ax.text(j, i - 0.12, f"{v:.3f}", ha="center", va="center", color=text_col, fontsize=8, fontweight="bold")
            ax.text(j, i + 0.20, f"{star(p)}  (n={int(n) if np.isfinite(n) else '?'})",
                    ha="center", va="center", color=text_col, fontsize=6.5)
    short = {o["label"]: o["key"] for o in OUTCOMES}
    ax.set_xticks(range(ncols)); ax.set_xticklabels([short.get(c, c) for c in oc_order], fontsize=9)
    ax.set_yticks(range(nrows)); ax.set_yticklabels(ex_order, fontsize=8)
    ax.set_title("Two-sample MR: haemoglobin/heme/iron proteins $\\rightarrow$ brain haemorrhage\n"
                 "IVW effect (colour), * p<0.05 ** p<0.01 *** p<0.001, ns = not significant\n"
                 "SAH = subarachnoid haemorrhage (I60); ICH = intracerebral haemorrhage (I61)",
                 fontsize=9, pad=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nrows, 1), minor=True)
    ax.grid(which="minor", color="white", lw=2)
    ax.tick_params(which="minor", length=0)
    fig.text(0.5, 0.02,
             "Real GWAS Catalog harmonised sumstats (European pQTL x European UKB-WGS outcomes). "
             "n = retained IVs; null cells reported as-is.",
             ha="center", fontsize=6.0, color="#444444")
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    png = os.path.join(FIG, "Fig24_mr_panel.png"); pdf = os.path.join(FIG, "Fig24_mr_panel.pdf")
    fig.savefig(png, dpi=300); fig.savefig(pdf)
    plt.close(fig)
    print("saved", png); print("saved", pdf)

if __name__ == "__main__":
    main()
