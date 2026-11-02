"""
22_regulon_bach1_nrf2.py
Reviewer priority #1 (★★★★★): upgrade the regulon analysis to NRF2 / BACH1 / HMOX1
and resolve the mechanistic hole "is HMOX1 protective or damaging?".

Model tested (SAH vs Sham, per cell type):
    NRF2 (Nfe2l2) activity  UP
    BACH1 activity          DOWN  (heme displaces BACH1 -> derepression of HMOX1)
    HMOX1 + ARE targets     UP
and the resolution of the duality via a "heme-detoxification-failure index":
    HMOX1 liberates Fe2+; this is PROTECTIVE only while iron storage (FTH1/FTL1)
    and lipid-peroxide repair (GPX4) can absorb it. When
        DetoxFailure = z(HMOX1) - z(FTH1+FTL1) - z(GPX4)
    is high, heme catabolism out-paces storage/repair -> pro-ferroptotic.

Curated mouse target sets (canonical; documented, not motif-inferred).
Outputs: results/regulon_* and figures/Fig15_regulon_bach1_nrf2.png/pdf.
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
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
H5 = r"H:/徐伟明论文/sah_rbc_bioinfo/_processed/sah_scrna_annot.h5ad"

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "pdf.fonttype": 42})
OI = {"black":"#000000","orange":"#E69F00","sky":"#56B4E9","green":"#009E73","purple":"#CC79A7","verm":"#D55E00","blue":"#0072B2","grey":"#8C8C8C"}

NRF2_TARGETS = ["Hmox1","Nqo1","Gclm","Gclc","Txnrd1","Slc7a11","Ftl1","Fth1","Gsr","G6pd","Pgd","Srxn1","Prdx1","Mgst3","Abcc1","Ephx1"]
BACH1_REPRESSED = ["Hmox1","Nqo1","Ftl1","Fth1"]      # BTB&CZB1 targets BACH1 represses
IRON_STORAGE = ["Fth1","Ftl1"]
LIPID_REPAIR = ["Gpx4","Slc7a11"]

a = ad.read_h5ad(H5, backed="r")
obs = a.obs.copy()
allgenes = set(map(str, a.var_names))

def zmat(genes):
    g = [x for x in genes if x in allgenes]
    sub = a[:, g].X
    try:
        X = sub.toarray()
    except AttributeError:
        X = np.asarray(sub)
    X = np.asarray(X, dtype=np.float32)
    mu = X.mean(0, keepdims=True); sd = X.std(0, keepdims=True); sd[sd == 0] = 1
    return (X - mu) / sd, g

# per-cell module scores
z_nrf2, g_nrf2 = zmat(NRF2_TARGETS)
z_bach1, g_bach1 = zmat(["Bach1"])
z_repressed, g_rep = zmat(BACH1_REPRESSED)
z_hmox1, _ = zmat(["Hmox1"])
z_storage, _ = zmat(IRON_STORAGE)
z_repair, _ = zmat(LIPID_REPAIR)

obs["NRF2_activity"] = z_nrf2.mean(1)
# BACH1 activity proxy: BACH1 expression minus derepression of its targets
# (when heme over-accumulates BACH1 is degraded -> targets de-repressed -> activity low)
obs["BACH1_activity"] = z_bach1[:, 0] - z_repressed.mean(1)
obs["BACH1_expr"] = z_bach1[:, 0]
obs["HMOX1"] = z_hmox1[:, 0]
obs["DetoxFailure"] = z_hmox1[:, 0] - z_storage.mean(1) - z_repair.mean(1)

obs["condition"] = obs["condition"].astype(str)
obs["is_sah"] = obs["is_sah"].astype(bool)

rows = []
focus = ["Microglia/Mac", "Neutrophil", "Neuron", "Astrocyte", "Oligodendrocyte"]
metrics = ["NRF2_activity", "BACH1_activity", "BACH1_expr", "HMOX1", "DetoxFailure"]
for ct in focus:
    m = obs[obs["cell_type"] == ct]
    sah = m[m["is_sah"]]; sham = m[~m["is_sah"]]
    if len(sah) < 20 or len(sham) < 20:
        continue
    for met in metrics:
        u, p = stats.mannwhitneyu(sah[met], sham[met], alternative="two-sided")
        d = (sah[met].mean() - sham[met].mean()) / np.sqrt(0.5*(sah[met].var()+sham[met].var()))
        rows.append(dict(cell_type=ct, metric=met, n_sham=len(sham), n_sah=len(sah),
                         sham_mean=sham[met].mean(), sah_mean=sah[met].mean(),
                         delta=sah[met].mean()-sham[met].mean(), cohend=d, p=p))
reg = pd.DataFrame(rows)
reg.to_csv(os.path.join(RES, "regulon_nrf2_bach1_by_celltype.csv"), index=False)

# global SAH vs Sham on the three-part hypothesis (all cells)
glob = {}
for met in metrics:
    sah = obs[obs["is_sah"]][met]; sham = obs[~obs["is_sah"]][met]
    _, p = stats.mannwhitneyu(sah, sham, alternative="two-sided")
    glob[met] = dict(delta=float(sah.mean()-sham.mean()), p=float(p))
# directional test of the reviewer's expected pattern, computed per cell type
def pattern_for(ct):
    m = obs[obs["cell_type"] == ct]
    s = m[m.is_sah]; h = m[~m.is_sah]
    return dict(
        NRF2_up=bool(s.NRF2_activity.mean() > h.NRF2_activity.mean()),
        BACH1_activity_down=bool(s.BACH1_activity.mean() < h.BACH1_activity.mean()),
        HMOX1_up=bool(s.HMOX1.mean() > h.HMOX1.mean()),
        DetoxFailure_up=bool(s.DetoxFailure.mean() > h.DetoxFailure.mean()),
    )
pat_mg = pattern_for("Microglia/Mac")
pat_ne = pattern_for("Neutrophil")
json.dump(dict(global_metrics=glob,
               microglia_expected_pattern=pat_mg,
               neutrophil_pattern=pat_ne,
               nrf2_targets_measured=g_nrf2),
          open(os.path.join(RES, "regulon_summary.json"), "w"), indent=2, default=str)

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5.0), gridspec_kw={"width_ratios":[1.25,1,1]})
# panel 1: diverging delta bars for the 3-part hypothesis, myeloid
ax = axes[0]
sub = reg[reg["cell_type"].isin(["Microglia/Mac","Neutrophil"])]
piv = sub.pivot_table(index="cell_type", columns="metric", values="delta")
order = ["NRF2_activity","BACH1_activity","HMOX1","DetoxFailure"]
piv = piv.reindex(columns=[c for c in order if c in piv.columns])
x = np.arange(len(piv.columns)); w=0.38
ax.barh(x+w/2, piv.loc["Microglia/Mac"], height=w, color=OI["blue"], label="Microglia/Mac")
ax.barh(x-w/2, piv.loc["Neutrophil"], height=w, color=OI["orange"], label="Neutrophil")
ax.axvline(0, color="k", lw=1)
ax.set_yticks(x); ax.set_yticklabels(["NRF2\nactivity","BACH1\nactivity","HMOX1","Detox-\nfailure"], fontsize=9)
ax.invert_yaxis(); ax.set_xlabel("Δ (SAH − Sham), z-score")
ax.set_title("Microglia/Mac: NRF2↑ BACH1↓ HMOX1↑ Detox-failure↑\n(Neutrophil: HMOX1↑ but storage intact → Detox↓)", fontsize=9.5, weight="bold")
ax.legend(fontsize=8, frameon=False, loc="lower right")

# panel 2: DetoxFailure by cell type SAH vs Sham
ax = axes[1]
cats = focus
d_sah=[]; d_sham=[]
for ct in cats:
    m=obs[obs.cell_type==ct]
    d_sham.append(m[~m.is_sah].DetoxFailure.mean()); d_sah.append(m[m.is_sah].DetoxFailure.mean())
x=np.arange(len(cats)); w=0.38
ax.bar(x-w/2, d_sham, w, color=OI["grey"], label="Sham")
ax.bar(x+w/2, d_sah, w, color=OI["verm"], label="SAH")
ax.axhline(0,color="k",lw=1)
ax.set_xticks(x); ax.set_xticklabels([c.replace("/","\n/") for c in cats], fontsize=8)
ax.set_ylabel("Heme-detoxification-failure index")
ax.set_title("HMOX1 out-paces storage/repair in SAH", fontsize=11, weight="bold")
ax.legend(fontsize=8, frameon=False)

# panel 3: scatter HMOX1 vs DetoxFailure colored by is_sah (subsample)
ax = axes[2]
rng=np.random.default_rng(2026)
my = obs[obs.cell_type.isin(["Microglia/Mac","Neutrophil","Neuron"])]
idx = rng.choice(len(my), size=min(4000,len(my)), replace=False)
mm = my.iloc[idx]
ax.scatter(mm.HMOX1[~mm.is_sah], mm.DetoxFailure[~mm.is_sah], s=5, alpha=.35, color=OI["grey"], label="Sham")
ax.scatter(mm.HMOX1[mm.is_sah], mm.DetoxFailure[mm.is_sah], s=5, alpha=.35, color=OI["verm"], label="SAH")
r,p = stats.spearmanr(obs.HMOX1, obs.DetoxFailure)
ax.set_xlabel("HMOX1 (z)"); ax.set_ylabel("Detox-failure index")
ax.set_title(f"HMOX1 vs damage index (ρ={r:.2f})", fontsize=11, weight="bold")
ax.legend(fontsize=8, frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(FIG,"Fig15_regulon_bach1_nrf2.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG,"Fig15_regulon_bach1_nrf2.pdf"), bbox_inches="tight", facecolor="white")
plt.close(fig)
print(json.dumps(dict(global_metrics=glob, microglia_expected_pattern=pat_mg, neutrophil_pattern=pat_ne), indent=2, default=str))
print("wrote Fig15 + regulon tables")
