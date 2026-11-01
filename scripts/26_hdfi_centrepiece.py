"""
26_hdfi_centrepiece.py
Reviewer round-7 #3 (★★★★★): make the Heme Detoxification Failure Index (HDFI)
the CENTREPIECE of the paper and show it PREDICTS neuronal ferroptosis.

    HDFI = z(HMOX1) - z(FTH1 + FTL) - z(GPX4)
      high  => haem catabolism (HMOX1) out-paces iron storage (ferritin) and
               lipid-peroxide repair (GPX4)  => pro-ferroptotic, not protective.

Computed in: (a) human acute blood GSE36791 (SAH vs control); (b) every cell
type of the mouse brain atlas (SAH vs Sham); (c) across time (Sham/24h/72h).
Key falsifiable test: does HDFI track a FerrDb ferroptosis-execution score?
  - within neurons (cell-level, n large)  and
  - sample-level microglia-HDFI vs neuron-ferroptosis (n=5, directional only)
plus the reviewer's Neuron-ferroptosis vs Microglia-HMOX1 coupling.
Outputs: results/hdfi_* and figures/Fig18_hdfi_centrepiece.png/pdf.
"""
import os, json
import numpy as np, pandas as pd
import anndata as ad
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
H5 = r"C:/Users/Admin/.qoderworkcn/workspace/muycr5n2rnyoywup/sah_rds/sah_scrna_annot.h5ad"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"pdf.fonttype":42})
OI = {"black":"#000000","orange":"#E69F00","sky":"#56B4E9","green":"#009E73","purple":"#CC79A7","verm":"#D55E00","blue":"#0072B2","grey":"#8C8C8C"}

DRIVER = ["Acsl4","Tfrc","Ncoa4","Lpcat3","Ptgs2","Alox15"]
SUPPR  = ["Gpx4","Slc7a11","Aifm2","Gch1","Dhodh","Dhrs3"]

# ---------- BULK HDFI (human) ----------
expr = pd.read_csv(os.path.join(DATA, "expr_linear_GSE36791.tsv"), sep="\t", index_col=0)
pheno = pd.read_csv(os.path.join(DATA, "pheno_GSE36791.tsv"), sep="\t").set_index("gsm")
common = expr.columns.intersection(pheno.index); expr = expr[common]
grp = pheno.loc[common, "group"].astype(str); is_sah = grp.str.contains("SAH", case=False).values
M = expr.values.T.astype(np.float32); genes = list(expr.index)
gi = {g:i for i,g in enumerate(genes)}
def zc(g):
    v = M[:, gi[g]] if g in gi else np.zeros(M.shape[0])
    return (v - v.mean())/(v.std()+1e-9)
bulk_hdfi = zc("HMOX1") - (zc("FTH1")+zc("FTL"))/2 - zc("GPX4")
bulk_df = pd.DataFrame(dict(Sample=common, SAH=is_sah, HDFI=bulk_hdfi))
s1 = bulk_df[bulk_df.SAH].HDFI; s0 = bulk_df[~bulk_df.SAH].HDFI
u, p = stats.mannwhitneyu(s1, s0, alternative="greater")   # one-sided: SAH higher
d = (s1.mean()-s0.mean())/np.sqrt(0.5*(s1.var()+s0.var()))
bulk_test = dict(delta=float(s1.mean()-s0.mean()), cohend=float(d), p_two=float(stats.mannwhitneyu(s1,s0,alternative='two-sided')[1]), p_greater=float(p))

# ---------- scRNA HDFI + neuron ferroptosis ----------
a = ad.read_h5ad(H5, backed="r")
need = sorted(set(["Hmox1","Fth1","Ftl1","Gpx4"]+DRIVER+SUPPR) & set(map(str,a.var_names)))
sub = a[:, need].X
try: X = sub.toarray()
except AttributeError: X = np.asarray(sub)
X = np.asarray(X, dtype=np.float32)
idx = {g:i for i,g in enumerate(need)}
Z = (X - X.mean(0))/(X.std(0)+1e-9)
def mz(gs): 
    c=[idx[g] for g in gs if g in idx]; return Z[:,c].mean(1) if c else np.zeros(X.shape[0])
obs = a.obs.copy(); obs["condition"]=obs["condition"].astype(str); obs["is_sah"]=obs["is_sah"].astype(bool)
obs["HDFI"] = Z[:,idx["Hmox1"]] - (Z[:,idx["Fth1"]]+Z[:,idx["Ftl1"]])/2 - Z[:,idx["Gpx4"]]
obs["FerroExec"] = mz(DRIVER) - mz([g for g in SUPPR if g in idx])
obs["HMOX1z"] = Z[:,idx["Hmox1"]]

# HDFI by cell type SAH vs Sham
focus = ["Microglia/Mac","Neutrophil","Neuron","Astrocyte","Oligodendrocyte"]
ct_rows=[]
for ct in focus:
    m=obs[obs.cell_type==ct]; s=m[m.is_sah]; h=m[~m.is_sah]
    _,pp=stats.mannwhitneyu(s.HDFI,h.HDFI,alternative="two-sided")
    ct_rows.append(dict(cell_type=ct, sham=h.HDFI.mean(), sah=s.HDFI.mean(), delta=s.HDFI.mean()-h.HDFI.mean(), p=pp))
ct_df=pd.DataFrame(ct_rows); ct_df.to_csv(os.path.join(RES,"hdfi_by_celltype.csv"),index=False)

# time: HDFI by condition (all cells + neuron)
time_rows=[]
for cond in ["Sham","SAH_24h","SAH_72h"]:
    m=obs[obs.condition==cond]
    time_rows.append(dict(condition=cond, all_HDFI=m.HDFI.mean(), neuron_HDFI=m[m.cell_type=="Neuron"].HDFI.mean(),
                          microglia_HDFI=m[m.cell_type=="Microglia/Mac"].HDFI.mean(), neuron_FerroExec=m[m.cell_type=="Neuron"].FerroExec.mean()))
time_df=pd.DataFrame(time_rows); time_df.to_csv(os.path.join(RES,"hdfi_time.csv"),index=False)

# KEY PREDICTION 1: within neurons, HDFI vs ferroptosis execution (cell-level)
neu = obs[obs.cell_type=="Neuron"]
rho_n, p_n = stats.spearmanr(neu.HDFI, neu.FerroExec)
# KEY PREDICTION 2: sample-level microglia-HDFI vs neuron-ferroptosis (only 3 conditions -> descriptive)
samp = (obs.groupby("sample", observed=True)
        .apply(lambda g: pd.Series(dict(
            mg_hdfi=g[g.cell_type=="Microglia/Mac"].HDFI.mean(),
            neu_ferro=g[g.cell_type=="Neuron"].FerroExec.mean(),
            neu_hdfi=g[g.cell_type=="Neuron"].HDFI.mean())), include_groups=False)
        .reset_index())
rho_s, p_s = stats.spearmanr(samp.mg_hdfi, samp.neu_ferro)
# reviewer coupling: neuron ferroptosis vs microglia HMOX1 (sample-level, n=3 conditions)
mg_hmox = obs[obs.cell_type=="Microglia/Mac"].groupby("sample", observed=True).HMOX1z.mean()
neu_ferro_s = obs[obs.cell_type=="Neuron"].groupby("sample", observed=True).FerroExec.mean()
j = pd.concat([mg_hmox.rename("mg_hmox1"), neu_ferro_s.rename("neu_ferro")], axis=1).dropna()
rho_c, p_c = stats.spearmanr(j.mg_hmox1, j.neu_ferro)

summary=dict(bulk_HDFI=bulk_test, neuron_HDFI_vs_ferroexec=dict(rho=float(rho_n),p=float(p_n),n=int(len(neu))),
             samplelevel_mgHDFI_vs_neuronFerro=dict(rho=float(rho_s),p=float(p_s),n=int(len(samp))),
             samplelevel_neuronFerro_vs_microgliaHMOX1=dict(rho=float(rho_c),p=float(p_c),n=int(len(j))),
             hdfi_by_celltype=ct_rows, hdfi_time=time_rows)
json.dump(summary, open(os.path.join(RES,"hdfi_summary.json"),"w"), indent=2, default=str)

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 4, figsize=(19,4.6))

# panel 1: bulk HDFI SAH vs Ctrl
ax=axes[0]
ax.boxplot([s0.values,s1.values], patch_artist=True, showfliers=False, widths=.55,
           boxprops=dict(facecolor="white"))
for patch,col in zip(ax.patches,[OI["grey"],OI["verm"]]): patch.set_facecolor(col); patch.set_alpha(.8)
ax.set_xticklabels(["Control","SAH"]); ax.set_ylabel("HDFI (blood)")
ax.set_title(f"Blood HDFI SAH vs Ctrl\n\u0394={bulk_test['delta']:.2f}, p={bulk_test['p_greater']:.3f}", fontsize=10.5, weight="bold")

# panel 2: HDFI by cell type (delta)
ax=axes[1]
cols=[OI["blue"],OI["orange"],OI["purple"],OI["green"],OI["sky"]]
ax.bar(range(len(ct_df)), ct_df["delta"], color=cols)
ax.axhline(0,color="k",lw=1); ax.set_xticks(range(len(ct_df)))
ax.set_xticklabels([c.replace("/","\n/") for c in ct_df.cell_type], fontsize=8)
ax.set_ylabel("\u0394 HDFI (SAH\u2212Sham)"); ax.set_title("HDFI shift by cell type", fontsize=10.5, weight="bold")

# panel 3: KEY - within-neuron HDFI vs ferroptosis execution
ax=axes[2]
rng=np.random.default_rng(2026); ii=rng.choice(len(neu), size=min(3000,len(neu)), replace=False)
nn=neu.iloc[ii]
ax.scatter(nn.HDFI[~nn.is_sah], nn.FerroExec[~nn.is_sah], s=5, alpha=.3, color=OI["grey"], label="Sham")
ax.scatter(nn.HDFI[nn.is_sah], nn.FerroExec[nn.is_sah], s=5, alpha=.3, color=OI["verm"], label="SAH")
zfit=np.polyfit(neu.HDFI, neu.FerroExec, 1); xs=np.linspace(neu.HDFI.min(),neu.HDFI.max(),10)
ax.plot(xs, np.polyval(zfit,xs), "--", color=OI["blue"], lw=1.6)
ax.set_xlabel("Neuron HDFI"); ax.set_ylabel("Neuron ferroptosis execution")
ax.set_title(f"HDFI predicts neuronal ferroptosis\n(cell-level \u03c1={rho_n:.2f}, n={len(neu)})", fontsize=10.5, weight="bold")
ax.legend(fontsize=8, frameon=False)

# panel 4: temporal trajectory of neuron HDFI and neuron ferroptosis execution
ax=axes[3]
td = time_df  # rows Sham/SAH_24h/SAH_72h
xx=np.arange(len(td))
ax.plot(xx, td.neuron_HDFI.values, "-o", color=OI["verm"], lw=2.2, label="Neuron HDFI")
ax.plot(xx, td.neuron_FerroExec.values, "-s", color=OI["blue"], lw=2.2, label="Neuron ferroptosis exec")
ax.axhline(0,color=OI["grey"],lw=.8)
ax.set_xticks(xx); ax.set_xticklabels([c.replace("SAH_","SAH-") for c in td.condition], fontsize=9)
ax.set_ylabel("z-score (per cell mean)")
ax.set_title(f"HDFI & ferroptosis track over time\n(microglia-HMOX1\u2194neuron-ferro \u03c1={rho_c:.2f}, n=3 cond)", fontsize=10.5, weight="bold")
ax.legend(fontsize=8, frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig18_hdfi_centrepiece.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG,"Fig18_hdfi_centrepiece.pdf"), bbox_inches="tight", facecolor="white")
plt.close(fig)
print(json.dumps({k:summary[k] for k in ["bulk_HDFI","neuron_HDFI_vs_ferroexec","samplelevel_mgHDFI_vs_neuronFerro","samplelevel_neuronFerro_vs_microgliaHMOX1"]}, indent=2, default=str))
print("wrote Fig18 + hdfi tables")
