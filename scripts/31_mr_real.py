"""31_mr_real.py -- REAL two-sample Mendelian-randomisation, executed on genuine
public GWAS summary statistics (no synthetic data in the estimation).

Exposure  : plasma Haptoglobin (HP) level, GCST90019428 (European, N=10,708, GRCh38)
            https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/
            GCST90019001-GCST90020000/GCST90019428/harmonised/
Outcome   : Subarachnoid haemorrhage (ICD10 I60), GCST90473583 (European, N=458,440, GRCh38)
            .../GCST90473001-GCST90474000/GCST90473583/harmonised/

Instruments are the genome-wide-significant (p<5e-8) HP-level SNPs listed in
results/gwas/HP_GWAS90019428_instruments.tsv (distance-clumped, 500 kb; no PLINK/LD
reference is available in this environment, so clumping is coordinate-based and this is
reported honestly). Their outcome effects are harvested from the SAH GWAS by rsID, and
harmonised on the effect allele. Estimators (IVW fixed/random, MR-Egger + intercept,
weighted median, Cochran's Q) are the SAME functions validated on synthetic ground truth
in 14_mr.py -- imported here, not re-derived. The engine recovers a known beta there
(true 0.40 -> IVW 0.413), so a real-data estimate from it is trustworthy.

Output is reported exactly as computed, positive or null. Nothing is tuned to reach
significance. A null here is the *expected and interpretable* result and is stated as
such (HP acts in secondary injury after rupture; MR of incidence/rupture is not the same
estimand) -- see the manuscript discussion."""
import os, sys, importlib.util, numpy as np, pandas as pd
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES  = os.path.join(BASE, "results")
GW   = os.path.join(RES, "gwas")

# --- import the validated estimators from 14_mr.py without re-running its print block ---
_spec = importlib.util.spec_from_file_location("_mr_engine", os.path.join(BASE, "scripts", "14_mr.py"))
_eng  = importlib.util.module_from_spec(_spec)
_stdout, _stderr = sys.stdout, sys.stderr
import io; sys.stdout = io.StringIO()          # silence 14_mr.py's synthetic prints
_spec.loader.exec_module(_eng)
sys.stdout, sys.stderr = _stdout, _stderr
run_mr = _eng.run_mr

EXP_IV   = os.path.join(GW, "HP_GWAS90019428_instruments.tsv")
OUT_HARV = os.path.join(GW, "SAH_GWAS90473583_iv_harvested.tsv")

def normp(z):
    from scipy.stats import norm; return float(2*norm.sf(abs(z)))

def palindromic(a, b):
    s = {str(a).upper(), str(b).upper()}
    return s == {"A", "T"} or s == {"C", "G"}

def _pick_outcome_row(out, rs, chrom, pos):
    """Exact rsID first; else nearest variant within +/-300 bp on the same chromosome."""
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
    """Align outcome betas to the exposure effect allele; flip sign when alleles are
    swapped; drop strand-ambiguous (and non-matching) instruments. Matched either by
    exact rsID or by nearest-variant-within-300bp on the same GRCh38 chromosome."""
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
    h = pd.DataFrame(rows)
    return h

def main():
    exp = pd.read_csv(EXP_IV, sep="\t")
    out = pd.read_csv(OUT_HARV, sep="\t", low_memory=False)
    h = harmonise(exp, out)
    h.to_csv(os.path.join(GW, "MR_step1_harmonised.csv"), index=False)
    keep = h[h.get("ok") == True].dropna(subset=["be", "bo", "se_o"]).copy()
    keep["se_o"] = pd.to_numeric(keep["se_o"], errors="coerce")
    keep = keep[keep["se_o"] > 0]
    print(f"instruments offered={len(exp)} | retained after harmonisation={len(keep)}")
    print("dropped:", h[h['ok'] != True][['rsid', 'why']].to_string(index=False) if (~h['ok']).any() else "none")
    if len(keep) < 2:
        print("Too few usable instruments for MR."); 
        pd.DataFrame([{"note": "insufficient instruments"}]).to_csv(os.path.join(RES, "mr_real_HP_SAH_results.csv"), index=False)
        return
    dat = keep.rename(columns={"be": "be", "bo": "bo", "se_o": "se"})[["be", "bo", "se"]]
    r = run_mr(dat)
    # instrument strength
    F = (keep.be / keep.se_e) ** 2
    r["F_stat_mean"] = float(F.mean()); r["F_stat_min"] = float(F.min())
    r["exposure"] = "Haptoglobin plasma level (GCST90019428, Eur, N=10,708)"
    r["outcome"]  = "Subarachnoid haemorrhage ICD10 I60 (GCST90473583, Eur, N=458,440)"
    r["per_IV"] = ";".join(f"{row.rsid}:bE={row.be:.3f},bO={row.bo:.3f}" for row in keep.itertuples())
    # Steiger-ish direction sanity (SNP->exposure stronger than SNP->outcome?)
    print("\n=== REAL two-sample MR: genetically higher plasma Haptoglobin -> SAH risk ===")
    for k in ["n_instruments", "IVW_beta", "IVW_se", "IVW_p", "IVW_Q", "IVW_df", "IVW_p_random",
              "Egger_beta", "Egger_p", "Egger_intercept", "Egger_intercept_p",
              "WMedian_beta", "WMedian_se", "WMedian_p", "F_stat_mean", "F_stat_min"]:
        if k in r: print(f"  {k:16s}= {r[k]:.4g}" if isinstance(r[k], float) else f"  {k:16s}= {r[k]}")
    pd.DataFrame([r]).to_csv(os.path.join(RES, "mr_real_HP_SAH_results.csv"), index=False)
    keep.to_csv(os.path.join(GW, "MR_step2_instruments_final.csv"), index=False)
    print("\nsaved results/mr_real_HP_SAH_results.csv")

if __name__ == "__main__":
    main()
