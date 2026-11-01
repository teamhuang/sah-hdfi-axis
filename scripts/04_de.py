"""04_de.py -- vectorized linear-model DE (per gene: expr ~ sex + group), BH-FDR.
limma-equivalent fixed effect (no empirical-Bayes shrinkage); documented as such.
Runs both cohorts; writes normalised log2 matrices + DEG tables + cross-cohort consistency."""
import os, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__)); BASE = os.path.dirname(HERE)
DAT, RES = os.path.join(BASE,"data"), os.path.join(BASE,"results")
os.makedirs(RES, exist_ok=True)

def qnorm(mat):
    """column-wise quantile normalise to common mean distribution (rank-based)."""
    ranks = np.argsort(np.argsort(mat, axis=0), axis=0)
    means = np.sort(mat, axis=0).mean(axis=1)
    return means[ranks]

def load(tag):
    e = pd.read_csv(os.path.join(DAT,f"expr_linear_{tag}.tsv"), sep="\t", index_col=0)
    p = pd.read_csv(os.path.join(DAT,f"pheno_{tag}.tsv"), sep="\t")
    p = p[p.gsm.isin(e.columns)].set_index("gsm").loc[e.columns]
    return e, p

def ols_de(lg, y_sex, y_grp):
    """Matrix OLS for all genes at once. Returns logFC(SAH vs Control), t, p."""
    # design: intercept, sex dummies (Female ref), grp(SAH vs Control ref)
    sexD = pd.get_dummies(y_sex, prefix="sex", drop_first=True).astype(float).values
    grp  = (y_grp == "SAH").astype(float)[:, None]      # 1=SAH,0=Control
    X = np.hstack([np.ones((lg.shape[1],1)), sexD, grp])       # n x p
    p = X.shape[1]
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ (X.T @ lg.T)                                 # p x n_genes
    resid = lg.T - X @ B
    dof = lg.shape[1] - p
    sigma2 = (resid**2).sum(axis=0) / dof
    grp_idx = X.shape[1]-1
    se = np.sqrt(np.clip(sigma2 * XtX_inv[grp_idx, grp_idx], 1e-300, None))
    logFC = B[grp_idx]; t = np.where(se>0, logFC/se, 0.0)
    pv = 2*stats.t.sf(np.abs(t), dof)
    pv = np.clip(np.nan_to_num(pv, nan=1.0, posinf=1.0, neginf=1.0), 0, 1)
    return logFC, t, pv, dof

def run(tag):
    e, p = load(tag)
    nn = int(e.isna().any(axis=1).sum())
    if nn: 
        e = e[~e.isna().any(axis=1)]
        print(f"{tag}: dropped {nn} genes with missing values")
    lg_raw = np.log2(e.values + 1.0)
    # filter: expressed in >=10% samples with log2>1
    keep = (lg_raw > 1).sum(axis=1) >= max(2, np.ceil(0.1*lg_raw.shape[1]))
    lg = qnorm(lg_raw[keep])
    genes = e.index.values[keep]
    sym = p["sex"].fillna("NA").values
    grp = p["group"].values
    logFC, t, pv, dof = ols_de(lg, sym, grp)
    padj = multipletests(pv, method="fdr_bh")[1]
    tt = pd.DataFrame({"gene":genes,"logFC":logFC,"AveExpr":lg.mean(axis=1),
                       "t":t,"P.Value":pv,"adj.P.Val":padj}).sort_values("P.Value")
    tt["n_SAH"]=(grp=="SAH").sum(); tt["n_Control"]=(grp=="Control").sum()
    sig = tt[(tt["adj.P.Val"]<0.05)&(tt["logFC"].abs()>1)]
    print(f"{tag}: n={lg.shape[1]} (SAH={ (grp=='SAH').sum()}, Ctrl={(grp=='Control').sum()}) genes={len(genes)} "
          f"| DEGs(FDR<.05,|lfc|>1)={len(sig)} up={ (sig.logFC>0).sum() } down={ (sig.logFC<0).sum() }")
    tt.to_csv(os.path.join(RES,f"DEG_{tag}.csv"), index=False)
    sig.to_csv(os.path.join(RES,f"DEGsig_{tag}.csv"), index=False)
    pd.DataFrame(genes).to_csv(os.path.join(DAT,f"genes_{tag}.txt"), index=False, header=False)
    np.save(os.path.join(DAT,f"lnorm_{tag}.npy"), lg)
    p.assign(gsm=p.index).to_csv(os.path.join(DAT,f"pheno_aligned_{tag}.csv"), index=False)
    return tt, sig

tt36, sig36 = run("GSE36791")
tt73, sig73 = run("GSE73378")
print("Top GSE36791 DEGs:")
print(tt36.head(12)[["gene","logFC","P.Value","adj.P.Val"]].to_string(index=False))

# cross-cohort consistency
m = tt36.merge(tt73, on="gene", suffixes=("_36","_73"))
mm = m[["logFC_36","logFC_73"]].dropna()
pear = stats.pearsonr(mm.logFC_36, mm.logFC_73)
sp   = stats.spearmanr(mm.logFC_36, mm.logFC_73)
conc = m[(m["adj.P.Val_36"]<0.05)&(m["logFC_36"].abs()>1)]
repro = conc[(conc["P.Value_73"]<0.05)&(np.sign(conc.logFC_73)==np.sign(conc.logFC_36))]
print(f"\ncommon genes={len(m)} | logFC Pearson r={pear[0]:.3f} (p={pear[1]:.2g}) "
      f"Spearman rho={sp.correlation:.3f} (p={sp.pvalue:.2g})")
print(f"GSE36791 DEGs testable in 73378={len(conc)}, concordant sig in 73378={len(repro)} "
      f"({100*len(repro)/max(1,len(conc)):.1f}%; chance~5%)")
pd.DataFrame({"metric":["n_common","pearson_r","pearson_p","spearman_rho","spearman_p",
                        "n_deg36_testable","n_reproduced","pct_reproduced","deg36","deg73"],
              "value":[len(m),pear[0],pear[1],sp.correlation,sp.pvalue,len(conc),len(repro),
                       100*len(repro)/max(1,len(conc)),len(sig36),len(sig73)]}
             ).to_csv(os.path.join(RES,"cross_cohort_consistency.csv"), index=False)
print("DONE 04_de")
