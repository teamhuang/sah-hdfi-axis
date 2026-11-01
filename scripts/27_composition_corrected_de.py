"""
27_composition_corrected_de.py
Reviewer round-7 #5 (★★★★★, "most important methodological addition"):
the whole-blood heme/HP signal could just be neutrophilia. Test whether HP,
HMOX1, the heme/iron module and the HDFI survive adjustment for leukocyte
composition.

Model per outcome y (gene z-score or module score):
    y ~ SAH + neutrophil_fraction + monocyte_fraction [+ T_fraction]
using marker-based deconvolution fractions (a CIBERSORT/MuSiC/BisqueRNA-style
reference-free surrogate: per-sample z-mean of canonical lineage markers,
softmax-normalised to pseudo-fractions). If the SAH coefficient stays
significant after adjustment, the signal is not merely composition.
Outputs: results/composition_adjusted_de.csv and figures/Fig19_composition_adjusted.png/pdf.
"""
import os, json
import numpy as np, pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"pdf.fonttype":42})
OI = {"grey":"#8C8C8C","verm":"#D55E00","blue":"#0072B2","green":"#009E73","orange":"#E69F00","purple":"#CC79A7"}

MARKERS = {
 "Neutrophil":["CSF3R","FCGR3B","CXCR2","ELANE","MPO","DEFA4","S100A8","S100A9","S100A12","CAMP","LTF","RETN"],
 "Monocyte":["CD14","FCN1","VCAN","LYZ","CCR2","MS4A7","ITGAM","S100A4","MNDA"],
 "T":["CD3D","CD3E","CD2","IL7R","CD28","TRAC","CD27","CCR7"],
 "B":["CD79A","CD79B","MS4A1","BANK1","CD19","IGHM"],
 "NK":["NKG7","GNLY","KLRD1","NCR1","PRF1","GZMB"],
 "Erythroid":["HBB","HBA1","HBA2","AHSP","ALAS2","GYPA","GYPB","SLC4A1"],
}
HEME_MODULE = ["HP","HMOX1","FTH1","FTL","SLC40A1","SLC25A37","MARC1","BLVRB","ALAS1","GYPB","NFE2"]

expr = pd.read_csv(os.path.join(DATA, "expr_linear_GSE36791.tsv"), sep="\t", index_col=0)
pheno = pd.read_csv(os.path.join(DATA, "pheno_GSE36791.tsv"), sep="\t").set_index("gsm")
common = expr.columns.intersection(pheno.index); expr = expr[common]
grp = pheno.loc[common, "group"].astype(str); sah = grp.str.contains("SAH", case=False).values.astype(float)
M = expr.values.T.astype(np.float32); genes = list(expr.index); gi = {g:i for i,g in enumerate(genes)}
Z = (M - M.mean(0))/(M.std(0)+1e-9)

def score(gs):
    c=[gi[g] for g in gs if g in gi]
    return Z[:, c].mean(1) if c else np.zeros(M.shape[0])

# pseudo-fractions: z-mean markers -> softmax across lineages
lin_scores = np.vstack([score(ms) for ms in MARKERS.values()])      # 6 x samples
e = np.exp(lin_scores - lin_scores.max(0, keepdims=True))
fracs = e / e.sum(0, keepdims=True)                                  # columns sum to 1
frac = pd.DataFrame(fracs.T, index=common, columns=list(MARKERS.keys()))
frac.to_csv(os.path.join(RES, "composition_pseudo_fractions.csv"))
frac_test = {ln: dict(delta=float(frac.loc[sah.astype(bool), ln].mean()-frac.loc[~sah.astype(bool), ln].mean())) for ln in frac.columns}

def adjusted(name, y, covars=("Neutrophil","Monocyte","T")):
    X = np.column_stack([sah] + [frac[c].values for c in covars])
    Xc = sm.add_constant(X)
    r0 = sm.OLS(y, sm.add_constant(sah)).fit()                 # unadjusted
    r1 = sm.OLS(y, Xc).fit()                                  # composition-adjusted
    return dict(outcome=name, beta_unadj=float(r0.params[1]), p_unadj=float(r0.pvalues[1]),
                beta_adj=float(r1.params[1]), p_adj=float(r1.pvalues[1]))

rows = []
for g in ["HP","HMOX1","FTH1","FTL","SLC40A1","MARC1","GPX4","SLC25A37","NFE2"]:
    if g in gi: rows.append(adjusted(g, Z[:, gi[g]]))
rows.append(adjusted("HEME_MODULE", score(HEME_MODULE)))
# HDFI outcome
hdfi = Z[:, gi["HMOX1"]] - (Z[:, gi["FTH1"]]+Z[:, gi["FTL"]])/2 - Z[:, gi["GPX4"]]
rows.append(adjusted("HDFI", hdfi))
df = pd.DataFrame(rows)
# BH across adjusted p-values
df["q_adj"] = multipletests(df["p_adj"], method="fdr_bh")[1]
df.to_csv(os.path.join(RES, "composition_adjusted_de.csv"), index=False)

surv = df[(df.q_adj < 0.1) & (df.beta_adj > 0)]["outcome"].tolist()
summary=dict(composition_shift=frac_test,
             n_outcomes=len(df),
             survive_adjustment_q0_10=surv,
             table=df.round(4).to_dict("records"))
json.dump(summary, open(os.path.join(RES,"composition_summary.json"),"w"), indent=2, default=str)

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 2, figsize=(13.5,5.2), gridspec_kw={"width_ratios":[1,1.1]})
# panel 1: pseudo-fraction shift SAH vs Ctrl
ax=axes[0]
labs=list(frac.columns); deltas=[frac_test[l]["delta"] for l in labs]
cols=[OI[c] for c in ["orange","green","blue","purple","grey","verm"]]
ax.bar(range(len(labs)), deltas, color=cols)
ax.axhline(0,color="k",lw=1); ax.set_xticks(range(len(labs)))
ax.set_xticklabels(labs, rotation=25, ha="right", fontsize=9)
ax.set_ylabel("\u0394 pseudo-fraction (SAH\u2212Ctrl)")
ax.set_title("Estimated blood composition shift", fontsize=11, weight="bold")

# panel 2: forest plot unadjusted vs adjusted p (log scale)
ax=axes[1]
y=np.arange(len(df))
ax.scatter(-np.log10(df.p_unadj.clip(lower=1e-16)), y, s=55, color=OI["grey"], label="unadjusted", zorder=3)
ax.scatter(-np.log10(df.p_adj.clip(lower=1e-16)), y, s=55, color=OI["verm"], label="composition-adjusted", zorder=4)
ax.axvline(-np.log10(0.05), color=OI["blue"], ls="--", lw=1.2)
ax.text(-np.log10(0.05)+0.05, len(df)-0.3, "p=0.05", color=OI["blue"], fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(df.outcome, fontsize=9)
ax.invert_yaxis(); ax.set_xlabel(r"$-\log_{10}$ SAH p-value")
ax.set_title("SAH effect before vs after composition adjustment", fontsize=11, weight="bold")
ax.legend(fontsize=8, frameon=False, loc="lower right")

fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig19_composition_adjusted.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG,"Fig19_composition_adjusted.pdf"), bbox_inches="tight", facecolor="white")
plt.close(fig)
print(df.round(4).to_string(index=False))
print("\nsurvive adjustment (q<0.1, positive):", surv)
print("composition shift:", {k: round(v["delta"],3) for k,v in frac_test.items()})
print("wrote Fig19 + composition_adjusted_de.csv")
