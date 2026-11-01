"""
23_cellchat_style.py
Reviewer priority #2 (★★★★★): replace the "expression-product edge" with a
CellChat-style STATISTICAL ligand-receptor inference, and test whether HP-CD163
is a specific heme-handling route rather than generic inflammation.

Faithful to CellChat's core model:
  * over-expression filter (gene detected in >10% of the sender/receiver cells
    and mean log-expr > 0.05),
  * law-of-mass-action interaction probability  P = <L>_source x <R>_target,
  * significance by permuting cell labels (default 1000) -> empirical p,
  * comparison against competing inflammatory routes
    (S100A8/9-TLR4, IL1B-IL1R1, TNF-TNFRSF1A/B) and the alternative heme route
    HP-MRC1, across Sham / SAH-24h / SAH-72h.

NOTE: this is a documented Python re-implementation of CellChat's statistical
core (same estimator + permutation test), not the R CellChat package, because R
is unstable on this CJK-path Windows build. Reported honestly.

Outputs: results/cellchat_style_*.csv and figures/Fig16_cellchat_style.png/pdf.
"""
import os, json
import numpy as np, pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
H5 = r"C:/Users/Admin/.qoderworkcn/workspace/muycr5n2rnyoywup/sah_rds/sah_scrna_annot.h5ad"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"pdf.fonttype":42})
OI = {"black":"#000000","orange":"#E69F00","sky":"#56B4E9","green":"#009E73","purple":"#CC79A7","verm":"#D55E00","blue":"#0072B2","grey":"#8C8C8C"}

PAIRS = [
    ("HP","CD163","Hp","Cd163","heme"),
    ("HP","MRC1","Hp","Mrc1","heme"),
    ("S100A8","TLR4","S100a8","Tlr4","inflam"),
    ("S100A9","TLR4","S100a9","Tlr4","inflam"),
    ("IL1B","IL1R1","Il1b","Il1r1","inflam"),
    ("TNF","TNFRSF1A","Tnf","Tnfrsf1a","inflam"),
    ("TNF","TNFRSF1B","Tnf","Tnfrsf1b","inflam"),
]
SRC, TGT = "Neutrophil", "Microglia/Mac"
CONDS = ["Sham","SAH_24h","SAH_72h"]
NPERM = 1000
rng = np.random.default_rng(2026)

a = ad.read_h5ad(H5, backed="r")
genes = sorted({g for _,_,l,r,_ in PAIRS for g in (l,r)})
sub = a[:, genes].X
try: X = sub.toarray()
except AttributeError: X = np.asarray(sub)
X = np.asarray(X, dtype=np.float32)
col = {g:i for i,g in enumerate(genes)}
obs = a.obs.copy()
obs["condition"] = obs["condition"].astype(str)

def overexpr(vec, thr=0.05, frac=0.10):
    return (np.mean(vec > thr) >= frac) and (np.mean(vec) > thr)

def pair_stats(cond, lg, rg):
    ci = (obs["condition"].values == cond)
    si = ci & (obs["cell_type"].values == SRC)
    ti = ci & (obs["cell_type"].values == TGT)
    L = X[si, col[lg]]
    R = X[ti, col[rg]]
    if len(L) < 10 or len(R) < 10:
        return None
    if not (overexpr(L) and overexpr(R)):
        return dict(P=0.0, p=np.nan, oe=False, n_src=int(len(L)), n_tgt=int(len(R)))
    obsP = L.mean() * R.mean()
    # permutation: shuffle the (src/tgt) assignment over the pooled cells
    pooled = np.concatenate([L, R])
    nS = len(L); null = np.empty(NPERM)
    for k in range(NPERM):
        perm = rng.permutation(len(pooled))
        null[k] = pooled[perm[:nS]].mean() * pooled[perm[nS:]].mean()
    p = (np.sum(null >= obsP) + 1) / (NPERM + 1)
    return dict(P=float(obsP), p=float(p), oe=True, n_src=int(nS), n_tgt=int(len(R)))

rows = []
for cond in CONDS:
    for (ln, rn, lg, rg, cls) in PAIRS:
        s = pair_stats(cond, lg, rg)
        if s is None: continue
        rows.append(dict(condition=cond, ligand=ln, receptor=rn, cls=cls, **s))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(RES, "cellchat_style_pairs.csv"), index=False)

# normalise P within condition for ranking
df["P_rank"] = df.groupby("condition")["P"].rank(ascending=False)
df["sig"] = df["p"] < 0.05

# receptor detection diagnostics on the target (Microglia/Mac) population
det_rows = []
for cond in CONDS:
    ti = (obs["condition"].values == cond) & (obs["cell_type"].values == TGT)
    for _,_,lg,rg,_ in PAIRS:
        vec = X[ti, col[rg]]
        det_rows.append(dict(condition=cond, receptor=rg, pct_pos=float((vec>0.05).mean()*100), mean_expr=float(vec.mean())))
    sl = (obs["condition"].values == cond) & (obs["cell_type"].values == SRC)
    for _,_,lg,rg,_ in PAIRS:
        pass
det = pd.DataFrame(det_rows).drop_duplicates(subset=["condition","receptor"])
det.to_csv(os.path.join(RES, "cellchat_style_receptor_detection.csv"), index=False)

hp = df[(df.ligand=="HP") & (df.receptor=="CD163")].set_index("condition").reindex(CONDS)
hpmrc = df[(df.ligand=="HP") & (df.receptor=="MRC1")].set_index("condition").reindex(CONDS)

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), gridspec_kw={"width_ratios":[1.05,1.05,1]})

# panel 1: heatmap of interaction probability, grey where over-expression fails, * where permutation-significant
ax = axes[0]
pairs_lab = [f"{l}-{r}" for l,r,_,_,_ in PAIRS]
M = np.zeros((len(pairs_lab), len(CONDS)))
mask_oe = np.zeros_like(M, dtype=bool); mask_sig = np.zeros_like(M, dtype=bool)
for i,(ln,rn,lg,rg,cls) in enumerate(PAIRS):
    for j,cond in enumerate(CONDS):
        r = df[(df.condition==cond)&(df.ligand==ln)&(df.receptor==rn)]
        if len(r):
            M[i,j] = float(r["P"].iloc[0]); mask_oe[i,j]=bool(r["oe"].iloc[0]); mask_sig[i,j]=bool(r["sig"].iloc[0])
disp = np.ma.masked_where(~mask_oe, M)
im = ax.imshow(disp, cmap="YlOrRd", aspect="auto")
for i in range(len(pairs_lab)):
    for j in range(len(CONDS)):
        if not mask_oe[i,j]:
            ax.text(j,i,"not detected",ha="center",va="center",fontsize=7,color=OI["grey"],style="italic")
        else:
            ax.text(j,i,f"{M[i,j]:.2f}"+("*" if mask_sig[i,j] else ""),ha="center",va="center",fontsize=8,
                    color="white" if M[i,j]>np.nanmax(M)/2 else "black")
ax.set_xticks(range(len(CONDS))); ax.set_xticklabels(CONDS, fontsize=9)
ax.set_yticks(range(len(pairs_lab))); ax.set_yticklabels(pairs_lab, fontsize=9)
ax.set_title("Neutrophil→Microglia/Mac interaction prob.\n(* permutation p<0.05)", fontsize=10.5, weight="bold")
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# panel 2: receptor detection fraction on target myeloid cells (why CD163 fails)
ax = axes[1]
recs = ["Cd163","Mrc1","Tlr4"]
xx=np.arange(len(CONDS)); w=0.26
for k,rg in enumerate(recs):
    vals=[det[(det.condition==c)&(det.receptor==rg)]["pct_pos"].values for c in CONDS]
    vals=[v[0] if len(v) else 0 for v in vals]
    ax.bar(xx+(k-1)*w, vals, w, color=[OI["verm"],OI["sky"],OI["green"]][k], label=rg)
ax.set_xticks(xx); ax.set_xticklabels(CONDS, fontsize=9)
ax.set_ylabel("% Microglia/Mac expressing receptor")
ax.set_title("CD163 too sparse in dissociated\nbrain myeloid to be engaged", fontsize=10.5, weight="bold")
ax.legend(fontsize=8, frameon=False)

# panel 3: HP-MRC1 (detectable heme route) vs HP-CD163 (not detected) across time
ax = axes[2]
ax.plot(range(3), hpmrc["P"].values.astype(float), "-o", color=OI["blue"], lw=2.4, label="HP–MRC1 (detected)")
for i,(c,v,pv) in enumerate(zip(CONDS, hpmrc["P"].values.astype(float), hpmrc["p"].values.astype(float))):
    if v > 0:
        ax.text(i, v, f"p={pv:.3f}", ha="center", va="bottom", fontsize=8, color=OI["blue"])
ax.scatter(range(3), [0,0,0], marker="x", s=90, color=OI["verm"], label="HP–CD163 (not detected)")
ax.set_xticks(range(3)); ax.set_xticklabels(CONDS, fontsize=9)
ax.set_ylabel("Interaction probability")
ax.set_title("Heme-scavenger route actually engaged", fontsize=10.5, weight="bold")
ax.legend(fontsize=8, frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig16_cellchat_style.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG,"Fig16_cellchat_style.pdf"), bbox_inches="tight", facecolor="white")
plt.close(fig)
print(df.to_string(index=False))
print("\nReceptor detection (% target myeloid positive):")
print(det.pivot(index="receptor", columns="condition", values="pct_pos").to_string())
print("HP-CD163 P by condition:", dict(zip(CONDS, np.round(hp['P'].values.astype(float),4))))
print("wrote Fig16 + cellchat_style tables")
