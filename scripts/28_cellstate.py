"""
28_cellstate.py
Reviewer round-7 #6.1: analyse CELL STATE, not just cell type.

Macrophage/microglia states (argmax of z-scored marker modules, cross-checked
with KMeans on the 3 module scores):
    homeostatic      P2ry12 Tmem119 Cx3cr1 Sall1 Tgfbr2 Hexb
    iron-handling    Hmox1 Fth1 Ftl1 Slc40a1 Blvrb Spi1
    ferroptosis-prone Acsl4 Ptgs2 Slc7a11 Tfrc Ncoa4  (low Gpx4)
Neutrophil states: normal vs heme-stressed (Hp Hmox1 Fth1 S100a8 S100a9 Lcn2).

Test whether SAH shifts the state composition (chi-square / proportion).
Outputs: results/cellstate_*.csv and figures/Fig20_cellstate.png/pdf.
"""
import os, json
import numpy as np, pandas as pd
import anndata as ad
from scipy import stats
from sklearn.cluster import KMeans
import matplotlib

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
H5 = r"H:/徐伟明论文/sah_rbc_bioinfo/_processed/sah_scrna_annot.h5ad"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],"pdf.fonttype":42})
OI = {"grey":"#8C8C8C","verm":"#D55E00","blue":"#0072B2","green":"#009E73","orange":"#E69F00","purple":"#CC79A7","sky":"#56B4E9"}

HOME = ["P2ry12","Tmem119","Cx3cr1","Sall1","Tgfbr2","Hexb","C1qa"]
IRON = ["Hmox1","Fth1","Ftl1","Slc40a1","Blvrb","Spi1"]
FERO = ["Acsl4","Ptgs2","Slc7a11","Tfrc","Ncoa4"]
HEMEN = ["Hp","Hmox1","Fth1","S100a8","S100a9","Lcn2"]

a = ad.read_h5ad(H5, backed="r")
need = sorted(set(HOME+IRON+FERO+HEMEN+["Gpx4"]) & set(map(str,a.var_names)))
sub = a[:, need].X
try: X = sub.toarray()
except AttributeError: X = np.asarray(sub)
X = np.asarray(X, dtype=np.float32)
idx = {g:i for i,g in enumerate(need)}
Z = (X - X.mean(0))/(X.std(0)+1e-9)
def mz(gs):
    c=[idx[g] for g in gs if g in idx]; return Z[:,c].mean(1) if c else np.zeros(X.shape[0])

obs = a.obs.copy(); obs["condition"]=obs["condition"].astype(str); obs["is_sah"]=obs["is_sah"].astype(bool)
obs["s_home"]=mz(HOME); obs["s_iron"]=mz(IRON); obs["s_fero"]=mz(FERO)-mz(["Gpx4"])
obs["s_hemeneut"]=mz(HEMEN)

# ---- macrophage states ----
mg = obs[obs.cell_type=="Microglia/Mac"].copy()
S = mg[["s_home","s_iron","s_fero"]].values
mg["state_argmax"] = np.array(["homeostatic","iron-handling","ferroptosis-prone"])[S.argmax(1)]
km = KMeans(n_clusters=3, n_init=10, random_state=2026).fit(S)
# map km labels to nearest centroid state
cent = km.cluster_centers_
lab2state = cent.argmax(1)
mg["state_km"] = pd.Series(lab2state[km.labels_]).map({0:"homeostatic",1:"iron-handling",2:"ferroptosis-prone"}).values

def prop_test(df, col):
    tab = pd.crosstab(df["is_sah"], df[col])
    chi2,p,dof,_ = stats.chi2_contingency(tab)
    # SAH-vs-Sham share of the damage states
    share = df.groupby("is_sah")[col].value_counts(normalize=True).unstack()
    return dict(chi2=float(chi2), p=float(p),
                sham_share={s: float(share.loc[False,s]) if s in share.columns else np.nan for s in ["homeostatic","iron-handling","ferroptosis-prone"]},
                sah_share ={s: float(share.loc[True ,s]) if s in share.columns else np.nan for s in ["homeostatic","iron-handling","ferroptosis-prone"]})
mg_arg = prop_test(mg, "state_argmax")
mg_km  = prop_test(mg, "state_km")
mg[["sample","condition","cell_type","is_sah","state_argmax","state_km","s_home","s_iron","s_fero"]].to_csv(os.path.join(RES,"cellstate_macrophage.csv"), index=False)

# ---- neutrophil heme-stress state ----
nu = obs[obs.cell_type=="Neutrophil"].copy()
thr = nu.s_hemeneut.median()
nu["state"] = np.where(nu.s_hemeneut > thr, "heme-stressed", "normal")
ntab = pd.crosstab(nu["is_sah"], nu["state"])
nchi2, np_ = stats.chi2_contingency(ntab)[:2]
nshare = nu.groupby("is_sah")["state"].value_counts(normalize=True).unstack()
neut = dict(chi2=float(nchi2), p=float(np_),
            sham_heme_stressed=float(nshare.loc[False,"heme-stressed"]) if "heme-stressed" in nshare.columns else np.nan,
            sah_heme_stressed=float(nshare.loc[True,"heme-stressed"]) if "heme-stressed" in nshare.columns else np.nan,
            mean_score_sham=float(nu[~nu.is_sah].s_hemeneut.mean()), mean_score_sah=float(nu[nu.is_sah].s_hemeneut.mean()))
_, mwu_p = stats.mannwhitneyu(nu[nu.is_sah].s_hemeneut, nu[~nu.is_sah].s_hemeneut, alternative="greater")
neut["mwu_p_greater"]=float(mwu_p)

summary=dict(macrophage_argmax=mg_arg, macrophage_kmeans=mg_km, neutrophil=neut)
json.dump(summary, open(os.path.join(RES,"cellstate_summary.json"),"w"), indent=2, default=str)

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 3, figsize=(15,4.8))
states=["homeostatic","iron-handling","ferroptosis-prone"]; scol=[OI["blue"],OI["orange"],OI["verm"]]

# panel 1: macrophage state share argmax
ax=axes[0]
xs=np.arange(len(states)); w=0.38
sh=[mg_arg["sham_share"][s] for s in states]; sa=[mg_arg["sah_share"][s] for s in states]
ax.bar(xs-w/2, np.array(sh)*100, w, color=OI["grey"], label="Sham")
ax.bar(xs+w/2, np.array(sa)*100, w, color=OI["verm"], label="SAH")
ax.set_xticks(xs); ax.set_xticklabels(["homeo","iron-\nhandling","ferropto-\nsis-prone"], fontsize=9)
ax.set_ylabel("% of microglia/macrophages"); ax.set_title(f"Macrophage state (argmax)\nchi2 p={mg_arg['p']:.1e}", fontsize=10.5, weight="bold")
ax.legend(fontsize=8, frameon=False)

# panel 2: KMeans cross-check
ax=axes[1]
sh=[mg_km["sham_share"][s] for s in states]; sa=[mg_km["sah_share"][s] for s in states]
ax.bar(xs-w/2, np.array(sh)*100, w, color=OI["grey"], label="Sham")
ax.bar(xs+w/2, np.array(sa)*100, w, color=OI["verm"], label="SAH")
ax.set_xticks(xs); ax.set_xticklabels(["homeo","iron-\nhandling","ferropto-\nsis-prone"], fontsize=9)
ax.set_title(f"Macrophage state (KMeans check)\nchi2 p={mg_km['p']:.1e}", fontsize=10.5, weight="bold")
ax.legend(fontsize=8, frameon=False)

# panel 3: neutrophil heme-stress score
ax=axes[2]
ax.boxplot([nu[~nu.is_sah].s_hemeneut.values, nu[nu.is_sah].s_hemeneut.values],
           patch_artist=True, showfliers=False, widths=.5)
for patch,col in zip(ax.patches,[OI["grey"],OI["orange"]]): patch.set_facecolor(col); patch.set_alpha(.8)
ax.set_xticklabels(["Sham","SAH"]); ax.set_ylabel("Heme-stress module score")
ax.set_title(f"Neutrophil heme-stressed state\nSAH {neut['sah_heme_stressed']*100:.0f}% vs Sham {neut['sham_heme_stressed']*100:.0f}% (p={neut['mwu_p_greater']:.1e})", fontsize=10.5, weight="bold")

fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig20_cellstate.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG,"Fig20_cellstate.pdf"), bbox_inches="tight", facecolor="white")
plt.close(fig)
print(json.dumps(summary, indent=2, default=str))
print("wrote Fig20 + cellstate tables")
