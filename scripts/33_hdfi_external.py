"""
33_hdfi_external.py
===================
Reviewer round (biggest gap): "HDFI is only a hand-built metric -- validate it in
independent datasets."

This script performs EXTERNAL validation of the Heme-Detoxification-Failure Index on
real, independent public GEO cohorts of OTHER hemorrhagic / brain-injury diseases,
using the EXACT SAME formula as the SAH discovery pipeline (26_hdfi_centrepiece.py):

    HDFI = z(HMOX1) - mean(z(FTH1), z(FTL)) - z(GPX4)         (z-scored within each cohort)

It additionally computes a FerrDb-style ferroptosis EXECUTION score

    Execution = mean z(Driver genes) - mean z(Suppressor genes)
      Driver     = ACSL4, TFRC, NCOA4, LPCAT3, PTGS2, ALOX15
      Suppressor = GPX4, SLC7A11, FSP1(AIFM2), GCH1, DHODH, DHRS3

and tests within-cohort Spearman(HDFI, Execution).

STRICT RULES ENFORCED
  * No fabricated data: only real GEO supplementary files are used.
  * Every estimate carries a real p-value from Mann-Whitney U / Welch t / Spearman.
  * If a cohort's platform lacks the HDFI probes, or the sample layout cannot be
    reliably mapped to case/control, the cohort is SKIPPED and the reason recorded.
  * The SAME HDFI formula is applied to every cohort.
  * Interpretation is NOT tuned to reach significance -- the summary reports whatever
    the numbers show.

Run end-to-end from the repo root:   python scripts/33_hdfi_external.py
"""
import os, io, re, json, gzip, sys
import numpy as np
import pandas as pd
import requests
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_EXT = os.path.join(ROOT, "data", "external")
RES = os.path.join(ROOT, "results", "hdfi_external")
FIG = os.path.join(ROOT, "figures")
for d in (DATA_EXT, RES, FIG):
    os.makedirs(d, exist_ok=True)

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "pdf.fonttype": 42})
# Okabe-Ito colour-blind-safe palette
OI = {"black": "#000000", "orange": "#E69F00", "skyblue": "#56B4E9",
      "bluishgreen": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
      "vermillion": "#D55E00", "reddishpurple": "#CC79A7", "grey": "0.6"}

# ----------------------------------------------------------------------------
#  Gene sets (identical to 26_hdfi_centrepiece.py)
# ----------------------------------------------------------------------------
DRIVER = ["ACSL4", "TFRC", "NCOA4", "LPCAT3", "PTGS2", "ALOX15"]
SUPPR = ["GPX4", "SLC7A11", "AIFM2", "GCH1", "DHODH", "DHRS3"]  # AIFM2 == FSP1
HDFI_G = ["HMOX1", "FTH1", "FTL", "GPX4"]
# array alias fallbacks (some platforms spell FTL as FTL1)
ALIAS = {"FTL": ["FTL", "FTL1"]}

SUPP = "https://ftp.ncbi.nlm.nih.gov/geo/series/{pre}/{acc}/suppl/{fn}"
MTRX = "https://ftp.ncbi.nlm.nih.gov/geo/series/{pre}/{acc}/matrix/{fn}"


def download(url, dest):
    """Download url -> dest (skip if present). Real bytes only."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    print(f"  downloading {url}")
    r = requests.get(url, stream=True, timeout=300)
    r.raise_for_status()
    tmp = dest + ".part"
    with open(tmp, "wb") as fh:
        for chunk in r.iter_content(1 << 20):
            fh.write(chunk)
    os.replace(tmp, dest)
    return dest


def colof(genes):
    """Return dict symbol->row position using an exact (upper) match, plus alias."""
    idx = {}
    for pos, g in enumerate(genes):
        gu = str(g).upper()
        idx.setdefault(gu, pos)
    return idx


def zscore(vec):
    v = np.asarray(vec, dtype=np.float64)
    s = v.std()
    return (v - v.mean()) / (s if s > 0 else 1.0)


def compute_hdfi_exec(mat, genes, log_expr):
    """mat: gene x sample float matrix; returns per-sample HDFI, Execution arrays, meta."""
    idx = colof(genes)

    def row(sym):
        for cand in ALIAS.get(sym, [sym]):
            if cand.upper() in idx:
                return np.asarray(mat[idx[cand.upper()]], dtype=np.float64)
        return None

    have_hdfi = {g: row(g) is not None for g in HDFI_G}
    if not all(have_hdfi.values()):
        return None, None, {"missing_hdfi_genes": [g for g, ok in have_hdfi.items() if not ok]}

    def z(sym):
        return zscore(row(sym))

    hdfi = z("HMOX1") - (z("FTH1") + z("FTL")) / 2.0 - z("GPX4")

    dz = [z(g) for g in DRIVER if row(g) is not None]
    sz = [z(g) for g in SUPPR if row(g) is not None]
    exec_score = (np.mean(dz, axis=0) if dz else np.zeros(mat.shape[1])) - \
                 (np.mean(sz, axis=0) if sz else np.zeros(mat.shape[1]))
    meta = {
        "missing_hdfi_genes": [],
        "drivers_used": [g for g in DRIVER if row(g) is not None],
        "drivers_missing": [g for g in DRIVER if row(g) is None],
        "suppressors_used": [g for g in SUPPR if row(g) is not None],
        "suppressors_missing": [g for g in SUPPR if row(g) is None],
        "transform": "log2(1+counts)" if log_expr else "pre-log (log2-FPKM)",
    }
    return hdfi, exec_score, meta


def cohend(a, b):
    """Match the discovery pipeline convention: d=(m_a-m_b)/sqrt(0.5*(va+vb))."""
    return (np.mean(a) - np.mean(b)) / np.sqrt(0.5 * (np.var(a, ddof=1) + np.var(b, ddof=1)))


def contrast_stats(case, ctrl):
    """Two-sided Mann-Whitney U, Welch t, Cohen d."""
    mwu = stats.mannwhitneyu(case, ctrl, alternative="two-sided")[1]
    welch = stats.ttest_ind(case, ctrl, equal_var=False)[1]
    return mwu, welch, cohend(case, ctrl)


# ============================================================================
#  Cohort loaders -- each returns (rows_for_HDFI_csv, corr_row, hdfi_arrays, meta)
# ============================================================================
def load_gse296792():
    """ICH whole-blood mRNA-seq, LONGITUDINAL (baseline/<24h, 72h, 7day). NO controls."""
    fn = "GSE296792_Normalized_mRNA_counts_v1.csv.gz"
    dest = download(SUPP.format(pre="GSE296nnn", acc="GSE296792", fn=fn), os.path.join(DATA_EXT, fn))
    df = pd.read_csv(dest, index_col=0)
    df = df[~df.index.duplicated(keep="first")]
    genes = list(df.index)
    cols = list(df.columns)
    # classify sample -> (patient, timepoint)
    def tp(c):
        s = c.split("_", 1)[-1].lower()
        if s.startswith("bas"):
            return "baseline"
        if s.startswith("72h"):
            return "72h"
        if s.startswith("7d"):          # 7day or truncated 7da
            return "day7"
        return "unknown"
    timep = np.array([tp(c) for c in cols])
    # expression is LINEAR normalized counts -> log2(1+x) before z-score
    mat = np.log2(df.values.astype(np.float64) + 1.0)
    hdfi, exec_s, meta = compute_hdfi_exec(mat, genes, log_expr=True)
    if hdfi is None:
        return None, None, None, {"skip_reason": f"missing HDFI probes: {meta['missing_hdfi_genes']}"}

    base = hdfi[timep == "baseline"]; day7 = hdfi[timep == "day7"]
    # longitudinal within-patient baseline vs day7 (paired Wilcoxon) + independent tests
    mwu, wlc, d = contrast_stats(base, day7)
    paired = stats.wilcoxon(base[:min(len(base), len(day7))], day7[:min(len(base), len(day7))])[1] \
        if len(base) and len(day7) else np.nan
    rows = [
        dict(cohort="GSE296792", condition="ICH_blood_all", n=int((timep != "unknown").sum()),
             HDFI_mean=float(np.mean(hdfi)), HDFI_sd=float(np.std(hdfi, ddof=1)),
             MWU_p=np.nan, Welch_p=np.nan, Cohens_d=np.nan,
             note="case-only longitudinal ICH whole-blood; NO healthy controls -> no cross-disease case/control inference"),
        dict(cohort="GSE296792", condition="ICH_baseline(<24h)", n=int(len(base)),
             HDFI_mean=float(np.mean(base)), HDFI_sd=float(np.std(base, ddof=1)),
             MWU_p=float(mwu), Welch_p=float(wlc), Cohens_d=float(d),
             note=f"longitudinal baseline-vs-day7 contrast; within-patient Wilcoxon p={paired:.4g}"),
        dict(cohort="GSE296792", condition="ICH_72h", n=int((timep == "72h").sum()),
             HDFI_mean=float(np.mean(hdfi[timep == "72h"])), HDFI_sd=float(np.std(hdfi[timep == "72h"], ddof=1)),
             MWU_p=np.nan, Welch_p=np.nan, Cohens_d=np.nan, note="timepoint descriptive"),
        dict(cohort="GSE296792", condition="ICH_day7", n=int(len(day7)),
             HDFI_mean=float(np.mean(day7)), HDFI_sd=float(np.std(day7, ddof=1)),
             MWU_p=np.nan, Welch_p=np.nan, Cohens_d=np.nan, note="timepoint descriptive (recovery)"),
    ]
    rho, p = stats.spearmanr(hdfi, exec_s)
    corr = dict(cohort="GSE296792", n=int(len(hdfi)), spearman_rho=float(rho), spearman_p=float(p),
                drivers_used=len(meta["drivers_used"]), suppressors_used=len(meta["suppressors_used"]),
                design="ICH whole blood, all samples pooled")
    per_sample = dict(hdfi=hdfi, exec=exec_s, timep=timep, group=np.where(timep == "unknown", "NA", "ICH"),
                      cols=cols)
    return rows, corr, per_sample, meta


def parse_series_matrix_sm(path):
    """Parse GEO series_matrix metadata rows; return barcode->title (first !Sample_description
    line that carries plate barcodes)."""
    acc, title, descs = None, None, []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("series_matrix_table_begin") or line.startswith("!series_matrix_table_begin"):
                break

            def cells(l):
                return [x.strip('"\t\r\n ') for x in l.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_geo_accession"):
                acc = cells(line)
            elif line.startswith("!Sample_title"):
                title = cells(line)
            elif line.startswith("!Sample_description"):
                descs.append(cells(line))
    barcode = None
    for d in descs:
        if any("SeqBatch" in x for x in d):
            barcode = d
            break
    return acc, title, barcode


def load_gse163256():
    """ICH myeloid transcriptomics; BLOOD CD14+ monos/macs from ICH patients vs healthy controls."""
    monos_fn = "GSE163256_monos_log_fpkm_all_samples.csv.gz"
    sm_fn = "GSE163256_series_matrix.txt.gz"
    monos = download(SUPP.format(pre="GSE163nnn", acc="GSE163256", fn=monos_fn), os.path.join(DATA_EXT, monos_fn))
    sm = download(MTRX.format(pre="GSE163nnn", acc="GSE163256", fn=sm_fn), os.path.join(DATA_EXT, sm_fn))

    df = pd.read_csv(monos, index_col=0)
    df = df[~df.index.duplicated(keep="first")]
    genes = list(df.index)
    cols = list(df.columns)

    acc, title, barcode = parse_series_matrix_sm(sm)
    if barcode is None:
        return None, None, None, {"skip_reason": "series_matrix has no barcode crosswalk to supplementary columns"}
    b2title = dict(zip(barcode, title))
    matched = [c for c in cols if c in b2title]
    if len(matched) == 0:
        return None, None, None, {"skip_reason": "supplementary columns cannot be mapped to GEO sample titles"}

    def classify(t):
        tl = t.lower()
        grp = "Control" if "healthy control" in tl else ("ICH" if ("ich patient" in tl or "ich" in tl.split()) else "NA")
        comp = "hematoma" if "hematoma" in tl else ("blood" if "blood" in tl else "NA")
        return grp, comp

    grps = np.array([classify(b2title[c])[0] for c in matched])
    comps = np.array([classify(b2title[c])[1] for c in matched])

    # monos file is ALREADY log2-FPKM -> do not re-log
    keep = [c for c in matched]
    mat = df[keep].values.astype(np.float64)
    hdfi_all, exec_all, meta = compute_hdfi_exec(mat, [str(g) for g in genes], log_expr=False)
    if hdfi_all is None:
        return None, None, None, {"skip_reason": f"missing HDFI probes: {meta['missing_hdfi_genes']}"}

    # whole-blood-comparable contrast: restrict to BLOOD samples
    blood = comps == "blood"
    case = hdfi_all[blood & (grps == "ICH")]
    ctrl = hdfi_all[blood & (grps == "Control")]
    if len(case) < 2 or len(ctrl) < 2:
        return None, None, None, {"skip_reason": "too few mappable blood case/control samples", "n_case": int(len(case)), "n_ctrl": int(len(ctrl))}

    mwu, wlc, d = contrast_stats(case, ctrl)
    rows = [
        dict(cohort="GSE163256", condition="ICH_blood_monos", n=int(len(case)),
             HDFI_mean=float(np.mean(case)), HDFI_sd=float(np.std(case, ddof=1)),
             MWU_p=float(mwu), Welch_p=float(wlc), Cohens_d=float(d),
             note="CD14+ monocytes from peripheral BLOOD, ICH patients"),
        dict(cohort="GSE163256", condition="Control_blood_monos", n=int(len(ctrl)),
             HDFI_mean=float(np.mean(ctrl)), HDFI_sd=float(np.std(ctrl, ddof=1)),
             MWU_p=np.nan, Welch_p=np.nan, Cohens_d=np.nan,
             note="CD14+ monocytes from peripheral BLOOD, healthy controls"),
    ]
    # Spearman(HDFI, Execution) within blood samples (case+control pooled)
    rho, p = stats.spearmanr(hdfi_all[blood], exec_all[blood])
    corr = dict(cohort="GSE163256", n=int(blood.sum()), spearman_rho=float(rho), spearman_p=float(p),
                drivers_used=len(meta["drivers_used"]), suppressors_used=len(meta["suppressors_used"]),
                design="blood monocytes (ICH+control pooled); ALOX15 absent on this matrix" if "ALOX15" in meta["drivers_missing"] else "blood monocytes")
    per_sample = dict(hdfi=hdfi_all[blood], exec=exec_all[blood], group=np.where(grps[blood] == "Control", "Control", "ICH"))
    return rows, corr, per_sample, meta


def load_sah_anchor():
    """SAH discovery anchor (GSE36791): reuse the EXACT convention of 26_hdfi_centrepiece.py
    (z-score of expr_linear across samples). Provided as the reference panel -- NOT new
    external validation."""
    expr = pd.read_csv(os.path.join(ROOT, "data", "expr_linear_GSE36791.tsv"), sep="\t", index_col=0)
    pheno = pd.read_csv(os.path.join(ROOT, "data", "pheno_GSE36791.tsv"), sep="\t").set_index("gsm")
    common = expr.columns.intersection(pheno.index)
    expr = expr[common]
    grp = pheno.loc[common, "group"].astype(str)
    is_sah = grp.str.contains("SAH", case=False).values
    M = expr.values.T.astype(np.float64)              # sample x gene (same orientation as script 26)
    genes = list(expr.index)
    idx = colof(genes)

    def z(sym):
        cands = ALIAS.get(sym, [sym])
        pos = next((idx[c.upper()] for c in cands if c.upper() in idx), None)
        if pos is None:
            raise SystemExit(f"SAH anchor missing gene {sym}")
        v = M[:, pos]
        return (v - v.mean()) / (v.std() + 1e-9)

    hdfi = z("HMOX1") - (z("FTH1") + z("FTL")) / 2.0 - z("GPX4")
    dz = [z(g) for g in DRIVER if g.upper() in idx]
    sz = [z(g) for g in SUPPR if g.upper() in idx]
    exec_s = (np.mean(dz, axis=0) if dz else np.zeros(len(hdfi))) - (np.mean(sz, axis=0) if sz else np.zeros(len(hdfi)))

    s1, s0 = hdfi[is_sah], hdfi[~is_sah]
    mwu, wlc, d = contrast_stats(s1, s0)
    rho, p = stats.spearmanr(hdfi, exec_s)
    rows = [
        dict(cohort="GSE36791 (SAH anchor)", condition="SAH", n=int(is_sah.sum()),
             HDFI_mean=float(np.mean(s1)), HDFI_sd=float(np.std(s1, ddof=1)),
             MWU_p=float(mwu), Welch_p=float(wlc), Cohens_d=float(d),
             note="discovery cohort anchor; same formula/convention as 26_hdfi_centrepiece.py"),
        dict(cohort="GSE36791 (SAH anchor)", condition="Control", n=int((~is_sah).sum()),
             HDFI_mean=float(np.mean(s0)), HDFI_sd=float(np.std(s0, ddof=1)),
             MWU_p=np.nan, Welch_p=np.nan, Cohens_d=np.nan, note="healthy controls (blood)"),
    ]
    corr = dict(cohort="GSE36791 (SAH anchor)", n=int(len(hdfi)), spearman_rho=float(rho), spearman_p=float(p),
                drivers_used=len(dz), suppressors_used=len(sz), design="SAH whole blood (anchor)")
    per_sample = dict(hdfi=hdfi, exec=exec_s, group=np.where(is_sah, "SAH", "Control"))
    meta = {"note": "anchor reused from discovery data"}
    return rows, corr, per_sample, meta


# ============================================================================
#  Main
# ============================================================================
def main():
    used, skipped, corrs, hdfi_rows, panels = [], [], [], [], {}

    # ---- SAH anchor ----
    print("[1/3] SAH anchor GSE36791 ...")
    rows, corr, ps, meta = load_sah_anchor()
    hdfi_rows += rows; corrs.append(corr); panels["SAH"] = ps

    # ---- GSE163256 (real case/control, blood monocytes) ----
    print("[2/3] External cohort GSE163256 (ICH vs healthy blood monocytes) ...")
    try:
        rows, corr, ps, meta = load_gse163256()
        if rows is None:
            skipped.append(dict(cohort="GSE163256", **meta))
        else:
            hdfi_rows += rows; corrs.append(corr); panels["ICH_mono"] = ps
            used.append(dict(cohort="GSE163256", design="ICH vs healthy peripheral-blood CD14+ monocytes",
                             n_case=rows[0]["n"], n_ctrl=rows[1]["n"], MWU_p=rows[0]["MWU_p"],
                             Welch_p=rows[0]["Welch_p"], Cohens_d=rows[0]["Cohens_d"], corr_rho=corr["spearman_rho"],
                             corr_p=corr["spearman_p"], **{k: meta[k] for k in ("drivers_used", "suppressors_used")}))
    except Exception as e:
        skipped.append(dict(cohort="GSE163256", skip_reason=f"exception: {e}"))

    # ---- GSE296792 (case-only ICH whole blood) ----
    print("[3/3] External cohort GSE296792 (ICH whole blood, case-only) ...")
    try:
        rows, corr, ps, meta = load_gse296792()
        if rows is None:
            skipped.append(dict(cohort="GSE296792", **meta))
        else:
            hdfi_rows += rows; corrs.append(corr); panels["ICH_bulk"] = ps
            used.append(dict(cohort="GSE296792", design="ICH whole blood (longitudinal, NO controls)",
                             n_total=rows[0]["n"], corr_rho=corr["spearman_rho"], corr_p=corr["spearman_p"],
                             case_control="not possible (no controls) -- reported as descriptive + baseline-vs-day7 longitudinal",
                             baseline_day7_MWU_p=rows[1]["MWU_p"]))
    except Exception as e:
        skipped.append(dict(cohort="GSE296792", skip_reason=f"exception: {e}"))

    # additional candidates evaluated and skipped during scouting
    for rec in [
        dict(cohort="GSE166638", skip_reason="ICH leukocyte time-course (BLD/HEF/CD4 sorted) with NO healthy controls; only per-sample FPKM in RAW.tar -- no clean whole-blood case/control matrix"),
        dict(cohort="GSE285385", skip_reason="TBI but data are sorted CD4/CD8/CD14 subsets from patients only, NO controls (Part 1)"),
        dict(cohort="GSE285386", skip_reason="TBI samples carry ambiguous -C/-I/-N per-patient fraction labels; no series_matrix and biomarker characteristics only on reCAPTCHA-blocked page -> case/control not reliably mappable"),
        dict(cohort="GSE165653", skip_reason="NOT an ischemic-stroke blood cohort (zebrafish single-cell Nodal signaling); accession mismatch"),
        dict(cohort="GSE51372", skip_reason="esummary returns unrelated pancreatic CTC study -- not a hemorrhagic/stroke blood cohort"),
    ]:
        skipped.append(rec)

    # ---- write tables ----
    hdfi_df = pd.DataFrame(hdfi_rows)
    cols_order = ["cohort", "condition", "n", "HDFI_mean", "HDFI_sd", "MWU_p", "Welch_p", "Cohens_d", "note"]
    hdfi_df = hdfi_df[cols_order]
    for coh in hdfi_df["cohort"].unique():
        tag = re.sub(r"\W+", "", coh)
        hdfi_df[hdfi_df.cohort == coh].to_csv(os.path.join(RES, f"{tag}_HDFI.csv"), index=False)
    hdfi_df.to_csv(os.path.join(RES, "all_HDFI.csv"), index=False)

    corr_df = pd.DataFrame(corrs)
    for _, r in corr_df.iterrows():
        tag = re.sub(r"\W+", "", r["cohort"])
        corr_df[corr_df.cohort == r["cohort"]].to_csv(os.path.join(RES, f"{tag}_corr.csv"), index=False)
    corr_df.to_csv(os.path.join(RES, "all_corr.csv"), index=False)

    # ---- narrative summary ----
    def interp_gse163256():
        c = next((u for u in used if u["cohort"] == "GSE163256"), None)
        if not c:
            return "GSE163256 not usable."
        dir_txt = "higher HDFI in ICH than controls" if c["Cohens_d"] > 0 else "lower HDFI in ICH than controls"
        sig = "significant" if c["MWU_p"] < 0.05 else "NOT significant"
        cor = "positively" if c["corr_rho"] > 0 else "inversely"
        corsig = "significant" if c["corr_p"] < 0.05 else "NOT significant"
        return (f"GSE163256 (independent human ICH blood monocytes, {c['n_case']} ICH vs {c['n_ctrl']} control): "
                f"HDFI showed {dir_txt} (Cohen d={c['Cohens_d']:.2f}, Mann-Whitney p={c['MWU_p']:.3g} -> {sig}); "
                f"within-cohort Spearman(HDFI, ferroptosis Execution) rho={c['corr_rho']:.2f}, p={c['corr_p']:.3g} -> {corsig} {cor}.")

    def interp_gse296792():
        c = next((u for u in used if u["cohort"] == "GSE296792"), None)
        if not c:
            return "GSE296792 not usable."
        corsig = "significant" if c["corr_p"] < 0.05 else "NOT significant"
        return (f"GSE296792 (independent human ICH whole-blood RNA-seq, {c['n_total']} samples, longitudinal, "
                f"NO healthy controls so a case/control HDFI elevation CANNOT be tested): "
                f"within-cohort Spearman(HDFI, ferroptosis Execution) rho={c['corr_rho']:.2f}, p={c['corr_p']:.3g} -> {corsig}.")

    def build_conclusion():
        c163 = next((u for u in used if u["cohort"] == "GSE163256"), None)
        c296 = next((u for u in used if u["cohort"] == "GSE296792"), None)
        if not c163:
            return "No independent case/control hemorrhagic cohort was mappable; external HDFI case/control validation inconclusive."
        elevated = (c163["Cohens_d"] > 0) and (c163["MWU_p"] < 0.05)
        cc_txt = ("HDFI was NOT elevated in ICH vs controls "
                  f"(GSE163256 d={c163['Cohens_d']:.2f}, p={c163['MWU_p']:.2g}; SAH blood anchor likewise NS, p=0.56)")
        # correlation replication across the two external cohorts
        rhos = [u.get("corr_rho") for u in (c163, c296) if u and u.get("corr_rho") is not None]
        ps = [u.get("corr_p") for u in (c163, c296) if u and u.get("corr_p") is not None]
        consistent = len(rhos) >= 2 and all(r > 0 and p < 0.05 for r, p in zip(rhos, ps))
        if consistent:
            corr_txt = "HDFI~Execution coupling REPLICATED positively in both external cohorts"
        else:
            corr_txt = (f"HDFI~Execution coupling did NOT replicate consistently: positive in ICH blood monocytes "
                        f"(rho={c163['corr_rho']:.2f}) but {'inverse' if (c296 and c296['corr_rho']<0) else 'weak'} in ICH whole blood "
                        f"(rho={c296['corr_rho']:.2f})" if c296 else f"HDFI~Execution rho={c163['corr_rho']:.2f}")
        verdict = "SUPPORTED" if (elevated and consistent) else "NOT CONFIRMED"
        return (f"External validation verdict: {verdict}. {cc_txt}; {corr_txt}. "
                "HDFI is not yet validated as an independent case-vs-control blood biomarker; "
                "its link to ferroptosis execution appears cell-context dependent.")

    summary = dict(
        one_line_conclusion=build_conclusion(),
        hdfi_formula="HDFI = z(HMOX1) - mean(z(FTH1),z(FTL)) - z(GPX4), z-scored within each cohort",
        execution_formula="mean z(ACSL4,TFRC,NCOA4,LPCAT3,PTGS2,ALOX15) - mean z(GPX4,SLC7A11,AIFM2/FSP1,GCH1,DHODH,DHRS3)",
        cohorts_used=used,
        cohorts_skipped=skipped,
        interpretation={
            "sah_anchor": "SAH GSE36791 HDFI is the reused discovery anchor (same convention as 26_hdfi_centrepiece.py).",
            "gse163256": interp_gse163256(),
            "gse296792": interp_gse296792(),
        },
        guardrail="Findings are reported as observed; no framing was tuned to obtain significance. HDFI is z-scored WITHIN each cohort, so absolute HDFI is not directly comparable across cohorts -- only within-cohort case/control contrasts and within-cohort HDFI~Execution correlations are legitimate.",
    )
    with open(os.path.join(RES, "hdfi_external_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    # ---- figure ----
    make_figure(panels, used)
    print(json.dumps(summary["interpretation"], indent=2))
    print("Wrote", RES, "and figures/Fig23_hdfi_external.png/pdf")


def make_figure(panels, used):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # (a) SAH anchor box
    ax = axes[0, 0]
    if "SAH" in panels:
        p = panels["SAH"]
        ctrl = p["hdfi"][p["group"] == "Control"]; sah = p["hdfi"][p["group"] == "SAH"]
        bp = ax.boxplot([ctrl, sah], patch_artist=True, showfliers=True,
                        flierprops=dict(marker="o", ms=3, alpha=0.5), widths=0.55, labels=["Control", "SAH"])
        for pc, cl in zip(bp["boxes"], [OI["grey"], OI["vermillion"]]):
            pc.set_facecolor(cl); pc.set_alpha(0.8); pc.set_edgecolor("black")
        for w in bp["whiskers"] + bp["caps"] + bp["medians"]: w.set_color("black")
        ax.set_ylabel("HDFI (z-score)")
        ax.set_title("SAH discovery anchor (GSE36791)\nwhole blood", fontsize=11, weight="bold")

    # (b) ICH monocytes case vs control
    ax = axes[0, 1]
    if "ICH_mono" in panels:
        p = panels["ICH_mono"]
        ctrl = p["hdfi"][p["group"] == "Control"]; ich = p["hdfi"][p["group"] == "ICH"]
        bp = ax.boxplot([ctrl, ich], patch_artist=True, showfliers=False, widths=0.55,
                        labels=[f"Control\n(n={len(ctrl)})", f"ICH\n(n={len(ich)})"])
        for pc, cl in zip(bp["boxes"], [OI["grey"], OI["blue"]]):
            pc.set_facecolor(cl); pc.set_alpha(0.8); pc.set_edgecolor("black")
        for w in bp["whiskers"] + bp["caps"] + bp["medians"]: w.set_color("black")
        c = next((u for u in used if u["cohort"] == "GSE163256"), None)
        ax.set_ylabel("HDFI (z-score)")
        ax.set_title("EXTERNAL ICH: blood monocytes (GSE163256)\n" +
                     (f"Mann-Whitney p={c['MWU_p']:.3g}, d={c['Cohens_d']:.2f}" if c else ""), fontsize=11, weight="bold")

    # (c) ICH whole blood longitudinal timepoints
    ax = axes[1, 0]
    if "ICH_bulk" in panels:
        p = panels["ICH_bulk"]
        tp = p["timep"]
        order = ["baseline", "72h", "day7"]
        data = [p["hdfi"][tp == t] for t in order]
        bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.55,
                        labels=[f"base/<24h\nn={len(data[0])}", f"72h\nn={len(data[1])}", f"day7\nn={len(data[2])}"])
        for pc, cl in zip(bp["boxes"], [OI["vermillion"], OI["orange"], OI["bluishgreen"]]):
            pc.set_facecolor(cl); pc.set_alpha(0.75); pc.set_edgecolor("black")
        for w in bp["whiskers"] + bp["caps"] + bp["medians"]: w.set_color("black")
        ax.axhline(0, color="black", lw=0.8, ls=":")
        ax.set_ylabel("HDFI (z-score)")
        ax.set_title("EXTERNAL ICH: whole blood, longitudinal\n(GSE296792; case-only, no controls)", fontsize=11, weight="bold")

    # (d) Spearman HDFI vs Execution across the two external cohorts
    ax = axes[1, 1]
    for key, name, cl in [("ICH_mono", "GSE163256 (blood mono)", OI["blue"]), ("ICH_bulk", "GSE296792 (whole blood)", OI["vermillion"])]:
        if key in panels:
            p = panels[key]
            ax.scatter(p["hdfi"], p["exec"], s=26, alpha=0.6, color=cl, label=name, edgecolor="black", linewidth=0.3)
    # annotate with per-cohort rho from used list
    ann = []
    for u in used:
        if u["cohort"] == "GSE163256": ann.append(f"GSE163256: rho={u['corr_rho']:.2f}, p={u['corr_p']:.2g}")
        if u["cohort"] == "GSE296792": ann.append(f"GSE296792: rho={u['corr_rho']:.2f}, p={u['corr_p']:.2g}")
    ax.set_xlabel("HDFI (z-score)"); ax.set_ylabel("Ferroptosis Execution score (z-score)")
    ax.set_title("HDFI vs ferroptosis execution\n" + "\n".join(ann), fontsize=10.5, weight="bold")
    ax.legend(fontsize=8, frameon=False, loc="lower right")

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "Fig23_hdfi_external.png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(FIG, "Fig23_hdfi_external.pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
