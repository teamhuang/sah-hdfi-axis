"""
24_ferroptosis_gsva.py
Reviewer priority #3 (★★★★★): turn the weak "ferroptosis NES 1.44" into a
FerrDb-categorised, GSVA-style score, and test the iron-handling -> neuronal-
damage chain spatially across cell types.

FerrDb categories (canonical members; measured subset only):
  Driver     : ACSL4 TFRC NCOA4 LPCAT3 PTGS2 ALOX15
  Suppressor : GPX4 SLC7A11 FSP1(AIFM2) GCH1 DHODH DHRS3
  Iron-metab : FTH1 FTL HMOX1 SLC40A1 STEAP3 HAMP ALAS2 BLVRB
Execution score = z(Driver) - z(Suppressor).
Detox-failure  = z(HMOX1) - z(FTH1+FTL) - z(GPX4).

Computed in BOTH human acute whole blood (GSE36791, SAH vs control) and the
mouse brain scRNA atlas (per cell type, SAH vs Sham), and the Driver score is
tested along the proposed Neutrophil -> Macrophage -> Neuron chain.
Outputs: results/ferroptosis_* and figures/Fig17_ferroptosis_gsva.png/pdf.
"""
import os, json
import numpy as np, pandas as pd
import anndata as ad
from scipy import stats
import matplotlib

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
H5 = r"H:/徐伟明论文/sah_rbc_bioinfo/_processed/sah_scrna_annot.h5ad"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"pdf.fonttype":42})
OI = {"black":"#000000","orange":"#E69F00","sky":"#56B4E9","green":"#009E73","purple":"#CC79A7","verm":"#D55E00","blue":"#0072B2","grey":"#8C8C8C"}

DRIVER = ["ACSL4","TFRC","NCOA4","LPCAT3","PTGS2","ALOX15"]
SUPPR  = ["GPX4","SLC7A11","FSP1","AIFM2","GCH1","DHODH","DHRS3"]
IRON   = ["FTH1","FTL","HMOX1","SLC40A1","STEAP3","HAMP","ALAS2","BLVRB"]

def zscore_sets(mat, genes, sets):
    """mat: samples x genes (log). returns dict name-> per-sample z-mean score."""
    Z = (mat - mat.mean(0)) / (mat.std(0) + 1e-9)
    idx = {g:i for i,g in enumerate(genes)}
    out = {}
    for name, s in sets.items():
        cols = [idx[g] for g in s if g in idx]
        out[name] = Z[:, cols].mean(1) if cols else np.zeros(mat.shape[0])
    return out

# ---------- BULK (human acute blood) ----------
expr = pd.read_csv(os.path.join(DATA, "expr_linear_GSE36791.tsv"), sep="\t", index_col=0)
pheno = pd.read_csv(os.path.join(DATA, "pheno_GSE36791.tsv"), sep="\t")
pheno = pheno.set_index(pheno.columns[0])
common = expr.columns.intersection(pheno.index)
expr = expr[common]
grp = pheno.loc[common, [c for c in pheno.columns if c.lower() == "group"][0]].astype(str)
is_sah = grp.str.contains("SAH", case=False).values
M = expr.values.T.astype(np.float32)   # samples x genes
genes = list(expr.index)               # gene symbols
bulk_sets = dict(Driver=DRIVER, Suppressor=SUPPR, Iron=IRON)  # zscore_sets keeps measured only
bs = zscore_sets(M, genes, bulk_sets)
bulk_df = pd.DataFrame(dict(Sample=common, SAH=is_sah,
                            Driver=bs["Driver"], Suppressor=bs["Suppressor"], Iron=bs["Iron"]))
bulk_df["Execution"] = bulk_df.Driver - bulk_df.Suppressor
# detox-failure on bulk
idx = {g:i for i,g in enumerate(genes)}
Z = (M - M.mean(0))/(M.std(0)+1e-9)
def zc(g): return Z[:, idx[g]] if g in idx else np.zeros(M.shape[0])
bulk_df["DetoxFailure"] = zc("HMOX1") - (zc("FTH1")+zc("FTL"))/2 - zc("GPX4")
bulk_df.to_csv(os.path.join(RES, "ferroptosis_bulk_scores.csv"), index=False)

bulk_tests = {}
for met in ["Driver","Suppressor","Iron","Execution","DetoxFailure"]:
    a1 = bulk_df.loc[bulk_df.SAH, met]; a0 = bulk_df.loc[~bulk_df.SAH, met]
    u,p = stats.mannwhitneyu(a1,a0,alternative="two-sided")
    d = (a1.mean()-a0.mean())/np.sqrt(0.5*(a1.var()+a0.var()))
    bulk_tests[met] = dict(delta=float(a1.mean()-a0.mean()), cohend=float(d), p=float(p))
rho_ei, p_ei = stats.spearmanr(bulk_df.Execution, bulk_df.Iron)
rho_di, p_di = stats.spearmanr(bulk_df.DetoxFailure, bulk_df.Iron)

# ---------- SINGLE-CELL (mouse brain) ----------
a = ad.read_h5ad(H5, backed="r")
mg = {  # mouse symbols
 "Driver":["Acsl4","Tfrc","Ncoa4","Lpcat3","Ptgs2","Alox15"],
 "Suppressor":["Gpx4","Slc7a11","Aifm2","Gch1","Dhodh","Dhrs3"],
 "Iron":["Fth1","Ftl1","Hmox1","Slc40a1","Steap3","Hamp","Alas2","Blvrb"],
}
allg = sorted({g for s in mg.values() for g in s} | {"Hmox1","Fth1","Ftl1","Gpx4"})
allg = [g for g in allg if g in set(map(str,a.var_names))]
sub = a[:, allg].X
try: X = sub.toarray()
except AttributeError: X = np.asarray(sub)
X = np.asarray(X, dtype=np.float32)
gidx = {g:i for i,g in enumerate(allg)}
Zc = (X - X.mean(0))/(X.std(0)+1e-9)
def zset(gl):
    cols=[gidx[g] for g in gl if g in gidx]
    return Zc[:,cols].mean(1) if cols else np.zeros(X.shape[0])
obs = a.obs.copy(); obs["condition"]=obs["condition"].astype(str); obs["is_sah"]=obs["is_sah"].astype(bool)
obs["FerroDriver"]=zset(mg["Driver"]); obs["FerroSuppr"]=zset(mg["Suppressor"]); obs["FerroIron"]=zset(mg["Iron"])
obs["FerroExec"]=obs.FerroDriver-obs.FerroSuppr
obs["DetoxFailure"]=Zc[:,gidx["Hmox1"]]-(Zc[:,gidx["Fth1"]]+Zc[:,gidx["Ftl1"]])/2-Zc[:,gidx["Gpx4"]]

chain = ["Neutrophil","Microglia/Mac","Neuron"]
chain_rows=[]
for ct in chain:
    m=obs[obs.cell_type==ct]
    for met in ["FerroDriver","FerroSuppr","FerroIron","FerroExec","DetoxFailure"]:
        s=m[m.is_sah][met]; h=m[~m.is_sah][met]
        _,p=stats.mannwhitneyu(s,h,alternative="two-sided")
        chain_rows.append(dict(cell_type=ct,metric=met,sham=h.mean(),sah=s.mean(),delta=s.mean()-h.mean(),p=p))
chain_df=pd.DataFrame(chain_rows); chain_df.to_csv(os.path.join(RES,"ferroptosis_scrna_chain.csv"),index=False)

# gradient of Execution across the chain (SAH only), Kruskal
sah=obs[obs.is_sah]
grp_scores=[sah[sah.cell_type==ct].FerroExec.values for ct in chain]
kw=stats.kruskal(*grp_scores)
# Spearman of Execution vs Iron across SAH cells
rho_sc,p_sc=stats.spearmanr(sah.FerroExec, sah.FerroIron)

summary=dict(bulk_tests=bulk_tests,
             bulk_exec_vs_iron=dict(rho=float(rho_ei),p=float(p_ei)),
             bulk_detox_vs_iron=dict(rho=float(rho_di),p=float(p_di)),
             chain_execution_mean_SAH={ct:float(sah[sah.cell_type==ct].FerroExec.mean()) for ct in chain},
             chain_kruskal_p=float(kw.pvalue),
             scrna_exec_vs_iron=dict(rho=float(rho_sc),p=float(p_sc)))
json.dump(summary, open(os.path.join(RES,"ferroptosis_summary.json"),"w"), indent=2, default=str)

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 3, figsize=(15.5,5.0), gridspec_kw={"width_ratios":[1,1,1.1]})

# panel 1: bulk SAH vs Ctrl boxplots of Driver, Suppressor, Execution
ax=axes[0]
mets=["Driver","Suppressor","Execution"]
xx=np.arange(len(mets)); w=0.35
for i,met in enumerate(mets):
    c0=bulk_df.loc[~bulk_df.SAH,met].values; c1=bulk_df.loc[bulk_df.SAH,met].values
    bp=ax.boxplot([c0,c1], positions=[xx[i]-w/2,xx[i]+w/2], widths=w, patch_artist=True, showfliers=False)
    for patch,col in zip(bp["boxes"],[OI["grey"],OI["verm"]]): patch.set_facecolor(col); patch.set_alpha(.8)
    for m in bp["medians"]: m.set_color("black")
ax.set_xticks(xx); ax.set_xticklabels([f"{m}\nSAH vs Ctrl" for m in mets], fontsize=9)
ax.axhline(0,color="k",lw=.8); ax.set_ylabel("FerrDb z-score")
ax.set_title("Human acute blood: FerrDb scores", fontsize=11, weight="bold")
ax.legend(handles=[plt.Rectangle((0,0),1,1,color=OI["grey"]),plt.Rectangle((0,0),1,1,color=OI["verm"])],
          labels=["Control","SAH"], fontsize=8, frameon=False)

# panel 2: scRNA Execution gradient along chain (SAH) with Kruskal p
ax=axes[1]
vals=[sah[sah.cell_type==ct].FerroExec.values for ct in chain]
bp=ax.boxplot(vals, patch_artist=True, showfliers=False, widths=.6)
for patch,col in zip(bp["boxes"],[OI["orange"],OI["blue"],OI["purple"]]): patch.set_facecolor(col); patch.set_alpha(.75)
for m in bp["medians"]: m.set_color("black")
ax.set_xticklabels(["Neutrophil","Macrophage","Neuron"], fontsize=9)
ax.set_ylabel("Ferroptosis execution score (SAH)")
ax.set_title(f"Neutrophil→Mac→Neuron chain\n(Kruskal p={kw.pvalue:.1e})", fontsize=11, weight="bold")

# panel 3: bulk Execution vs Iron scatter
ax=axes[2]
c0=~bulk_df.SAH.values
ax.scatter(bulk_df.Iron[c0], bulk_df.Execution[c0], s=28, color=OI["grey"], label="Control")
ax.scatter(bulk_df.Iron[~c0], bulk_df.Execution[~c0], s=28, color=OI["verm"], label="SAH")
# fit line on all
z=np.polyfit(bulk_df.Iron, bulk_df.Execution, 1); xs=np.linspace(bulk_df.Iron.min(),bulk_df.Iron.max(),10)
ax.plot(xs, np.polyval(z,xs), "--", color=OI["blue"], lw=1.5)
ax.set_xlabel("FerrDb iron-metabolism score"); ax.set_ylabel("Ferroptosis execution score")
ax.set_title(f"Iron load tracks ferroptosis\n(blood ρ={rho_ei:.2f}, p={p_ei:.1e})", fontsize=11, weight="bold")
ax.legend(fontsize=8, frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig17_ferroptosis_gsva.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG,"Fig17_ferroptosis_gsva.pdf"), bbox_inches="tight", facecolor="white")
plt.close(fig)
print(json.dumps(summary, indent=2, default=str))
print("wrote Fig17 + ferroptosis tables")
