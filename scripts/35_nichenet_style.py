"""
35_nichenet_style.py

Reviewer priority (★★★★): go ONE STEP BEYOND the CellChat-style receptor routes
(23_cellchat_style.py, sec.4.8) to ligand -> TARGET-GENE causality.

    "Which neutrophil-derived ligand drives microglial HMOX1/FTH1/FTL1/SLC40A1
     induction?"

This is a documented, minimal NicheNet-STYLE approximation implemented in pure
Python/scanpy. It is NOT the full Hadush/NicheNet package: there is NO
Connectome protein-signalling network and NO Seurat integration. Every output
is labelled "NicheNet-style" for the paper's integrity.

Pipeline (mirrors NicheNet's three-part logic: ligand activity from
(a) expression + (b) prior knowledge + (c) empirical target correlation):
  Step 1  sender (Neutrophil) ligand expression screen  -- over-expression
          filter identical to sec.4.8 (detected in >10% SAH neutrophils AND
          mean log-expr > 0.05) + SAH-vs-Sham upregulation (|log2FC|>0.3, MWU).
  Step 2  cognate receptor detection on receiver (Microglia/Mac = merged
          microglia+macrophage myeloid receiver), via an embedded curated
          ligand-receptor table (mouse symbols).
  Step 3  receiver response set = the four reviewer-named heme-handling genes
          HMOX1/FTH1/FTL1/SLC40A1 (mean z per cell); GPX4 reported in the
          OPPOSITE direction only.
  Step 4  curated prior ligand->target regulatory-potential weight (0-3).
          This is an honest hand-curated prior, NOT a genome-wide network.
  Step 5  empirical ligand-target activity = Spearman rho across SAH receiver
          cells between (a) local ligand exposure proxy = mean ligand log-expr
          in the K=25 nearest neutrophil sender cells in UMAP space and
          (b) the receiver target-gene score.
  Step 6  combined ranking = normalised prior * normalised sender log2FC z *
          normalised empirical-activity z.

Outputs:
  results/nichenet_style_scores.csv   per-ligand full component table
  figures/Fig25_nichenet_style.png/.pdf  ranked ligand dot/bar plot
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures")
H5 = r"C:/Users/Admin/.qoderworkcn/workspace/muycr5n2rnyoywup/sah_rds/sah_scrna_annot.h5ad"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "pdf.fonttype": 42})
# Okabe-Ito colour-blind-safe palette
OI = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
      "purple": "#CC79A7", "verm": "#D55E00", "blue": "#0072B2", "grey": "#8C8C8C"}

SRC = "Neutrophil"
Rcv = "Microglia/Mac"                       # merged microglia + macrophage receiver
CONDS = ["Sham", "SAH_24h", "SAH_72h"]
SAH = ["SAH_24h", "SAH_72h"]
NPERM = 0                                   # (no permutation here; empirical = Spearman)
K_NN = 25                                   # nearest sender neighbours in UMAP space

# ---- Step 2: curated ligand -> receptor(s) table (mouse symbols present in atlas) ----
LR = {
    "Tnf":    ["Tnfrsf1a", "Tnfrsf1b"],
    "Il1b":   ["Il1r1"],
    "S100a8": ["Tlr4", "Ager"],
    "S100a9": ["Tlr4", "Ager"],
    "Hmgb1":  ["Tlr4", "Ager", "Cxcr4"],
    "Lcn2":   ["Slc22a6"],
    "Spp1":   ["Cd44", "Itgav", "Itgb3"],
    "Il6":    ["Il6ra", "Il6st"],
    "Ccl2":   ["Ccr2"],
    "Tgfb1":  ["Tgfbr1", "Tgfbr2"],
    "Hp":     ["Cd163", "Mrc1"],
    "Mif":    ["Cd74", "Cxcr4"],
    "Csf1":   ["Csf1r"],
    "Ereg":   ["Egfr"],
    "Areg":   ["Egfr"],
    "Hbegf":  ["Egfr"],
    "Vegfa":  ["Kdr", "Flt1"],
    "Fgf2":   ["Fgfr1"],
    "Pdgfb":  ["Pdgfrb"],
    "Ccl5":   ["Ccr5"],
    "Anxa1":  ["Fpr1", "Fpr2"],
}
# ---- Step 4: curated prior ligand->target regulatory potential (0-3) ----
# 3 = strong NF-kB/MAPK drive of HMOX1/FTH1/FTL/SLC40A1; 2 = intermediate; 1 = weak/indirect
PRIOR = {
    "Tnf": 3, "Il1b": 3, "S100a8": 3, "S100a9": 3, "Hmgb1": 3,
    "Hp": 2, "Il6": 2, "Mif": 2,
    "Lcn2": 1, "Spp1": 1, "Ccl2": 1, "Csf1": 1, "Anxa1": 1, "Ccl5": 1,
    "Ereg": 1, "Areg": 1, "Hbegf": 1, "Vegfa": 1, "Fgf2": 1, "Tgfb1": 1,
    "Pdgfb": 0,
}
# ---- Step 3: heme-handling response set (reviewer-named) + opposite-direction GPX4 ----
RESPONSE = ["Hmox1", "Fth1", "Ftl1", "Slc40a1"]
OPPOSITE = ["Gpx4"]

# ---------------- load ----------------
a = ad.read_h5ad(H5, backed="r")
ligands = list(LR.keys())
all_recs = sorted({r for rs in LR.values() for r in rs})
need = sorted(set(ligands) | set(all_recs) | set(RESPONSE) | set(OPPOSITE))
need_present = [g for g in need if g in set(map(str, a.var_names))]
missing = [g for g in need if g not in need_present]
sub = a[:, need_present].X
try: Xs = sub.toarray()
except AttributeError: Xs = np.asarray(sub)
Xs = np.asarray(Xs, dtype=np.float32)
gi = {g: i for i, g in enumerate(need_present)}

obs = a.obs.copy()
obs["condition"] = obs["condition"].astype(str)
umap = np.asarray(a.obsm["X_umap"], dtype=np.float64)

# ---------------- receiver response score (per cell) ----------------
ri = (obs["cell_type"].values == Rcv)
resp_cols = [gi[g] for g in RESPONSE if g in gi]
recv_expr = Xs[np.ix_(ri, resp_cols)]
recv_z = (recv_expr - recv_expr.mean(0)) / (recv_expr.std(0) + 1e-9)
resp_score_recv = recv_z.mean(1)             # length = n receiver cells
opp_cols = [gi[g] for g in OPPOSITE if g in gi]
gpx4_recv = Xs[np.ix_(ri, opp_cols)][:, 0] if opp_cols else np.full(ri.sum(), np.nan)

recv_idx_all = np.where(ri)[0]               # absolute cell indices for receiver cells
cond_recv = obs["condition"].values[ri]      # condition label per receiver cell

# ---------------- sender pool (all neutrophils) for KNN exposure ----------------
neut_idx_all = np.where(obs["cell_type"].values == SRC)[0]
neut_umap = umap[neut_idx_all]
neut_cond = obs["condition"].values[neut_idx_all]
tree = cKDTree(neut_umap)

# receiver UMAP coords (query points)
recv_umap = umap[recv_idx_all]
# precompute the K=25 nearest neutrophil rows for every receiver cell ONCE
k = min(K_NN, len(neut_idx_all))
_, nn_local = tree.query(recv_umap, k=k)     # nn_local: index into neut_idx_all
if nn_local.ndim == 1:
    nn_local = nn_local[:, None]
nn_abs = neut_idx_all[nn_local]              # absolute cell indices of neighbours

# ---------------- Step 1: sender expression screen ----------------
neut_is_sah = np.isin(neut_cond, SAH)
neut_is_sham = neut_cond == "Sham"
rows = []
for L in ligands:
    if L not in gi:
        continue
    vec = Xs[:, gi[L]]
    v_sah = vec[neut_idx_all][neut_is_sah]
    v_sham = vec[neut_idx_all][neut_is_sham]
    frac_pos = float((v_sah > 0.05).mean())
    mean_sah = float(v_sah.mean())
    oe_pass = (frac_pos >= 0.10) and (mean_sah > 0.05)
    if len(v_sham) >= 5 and len(v_sah) >= 5:
        lfc = float(np.mean(v_sah) - np.mean(v_sham))
        mwu_p = float(stats.mannwhitneyu(v_sah, v_sham, alternative="two-sided")[1])
    else:
        lfc, mwu_p = np.nan, np.nan
    up_pass = (abs(lfc) > 0.3) and (mwu_p < 0.05) if not np.isnan(lfc) else False
    # receptor detection on receiver
    recs = [r for r in LR[L] if r in gi]
    rec_detect = {}
    for r in recs:
        rv = Xs[ri, gi[r]]
        rec_detect[r] = float((rv > 0.05).mean() * 100)
    recv_frac_pos = float(np.mean(list(rec_detect.values()))) if rec_detect else 0.0
    # empirical ligand-target activity (SAH receiver cells only)
    lig_neut = vec[neut_idx_all]                       # per-neutrophil ligand log-expr
    nbr_vals = lig_neut[nn_local]                       # (n_recv, k)
    exposure_recv = nbr_vals.mean(1)                   # local ligand exposure proxy per receiver cell
    sah_mask = np.isin(cond_recv, SAH)
    if sah_mask.sum() > 10:
        rho, pval = stats.spearmanr(exposure_recv[sah_mask], resp_score_recv[sah_mask])
    else:
        rho, pval = np.nan, np.nan
    rows.append(dict(
        ligand=L, sender_pct_pos_sah=round(frac_pos * 100, 2), sender_mean_logexpr_sah=round(mean_sah, 4),
        oe_pass=bool(oe_pass), sender_log2FC_sah_vs_sham=round(lfc, 4) if not np.isnan(lfc) else np.nan,
        sender_mwu_p=round(mwu_p, 6) if not np.isnan(mwu_p) else np.nan, up_pass=bool(up_pass),
        receptors=";".join(recs), receptor_detect_pct_mean=round(recv_frac_pos, 2),
        receptor_detection=rec_detect,
        prior_weight=PRIOR.get(L, 0),
        emp_activity_rho=round(float(rho), 4) if not np.isnan(rho) else np.nan,
        emp_activity_p=round(float(pval), 6) if not np.isnan(pval) else np.nan,
        n_sah_receiver=int(sah_mask.sum()),
    ))
df = pd.DataFrame(rows)

# ---------------- candidate gate + Step 6 combined ranking ----------------
cand = df[df.oe_pass].copy()                            # must pass over-expression in SAH neutrophils

def minmax(s):
    s = pd.to_numeric(s, errors="coerce")
    lo, hi = s.min(), s.max()
    if not np.isfinite(lo) or hi == lo:
        return pd.Series(np.full(len(s), 0.5), index=s.index)
    return (s - lo) / (hi - lo)

max_prior = max(PRIOR.values()) if PRIOR else 1
cand["n_prior"] = pd.to_numeric(cand["prior_weight"]) / max_prior
cand["n_sender"] = minmax(cand["sender_log2FC_sah_vs_sham"])
cand["n_emp"] = minmax(cand["emp_activity_rho"])
# honest: ligands with prior weight 0 cannot rank (regulatory-potential gate)
cand["combined_score"] = cand["n_prior"] * cand["n_sender"] * cand["n_emp"]
cand = cand.sort_values("combined_score", ascending=False).reset_index(drop=True)

# export (flatten receptor_detection dict to readable string)
out = cand.copy()
out["receptor_detection"] = out["receptor_detection"].apply(
    lambda d: "; ".join(f"{k}:{v:.1f}%" for k, v in d.items()))
out.to_csv(os.path.join(RES, "nichenet_style_scores.csv"), index=False)

# ---------------- printed summary (sanity check) ----------------
print("MISSING genes handled (dropped):", missing if missing else "none")
print("Candidate ligands passing over-expression gate:", len(cand), "of", len(df))
print("\n=== NicheNet-STYLE ranking: neutrophil ligands -> microglial/macro heme-handling response ===")
top = cand.head(10)
for _, r in top.iterrows():
    sig = "yes" if (pd.notna(r.emp_activity_p) and r.emp_activity_p < 0.05) else "no"
    print(f"{r.ligand:<8} score={r.combined_score:.3f} | prior={int(r.prior_weight)} "
          f"| sender_log2FC={r.sender_log2FC_sah_vs_sham:+.3f}(mwu p={r.sender_mwu_p:.3g}) "
          f"| receiver_rec%={r.receptor_detect_pct_mean:.1f} "
          f"| emp_rho={r.emp_activity_rho:+.3f}(p={r.emp_activity_p:.3g}, sig={sig})")
print("\nTOP-5 component breakdown:")
print(cand.head(5)[["ligand", "sender_pct_pos_sah", "sender_log2FC_sah_vs_sham", "sender_mwu_p",
                    "receptor_detect_pct_mean", "prior_weight", "emp_activity_rho", "emp_activity_p",
                    "n_prior", "n_sender", "n_emp", "combined_score"]].to_string(index=False))

# ---------------- figure: ranked dot/bar plot ----------------
show = cand.head(9).copy()
n = len(show)
fig, ax = plt.subplots(figsize=(8.4, 0.52 * n + 1.9))
y = np.arange(n)[::-1]
# empirical-activity significance drives bar colour
cols = [OI["verm"] if (pd.notna(p) and p < 0.05) else OI["blue"] for p in show["emp_activity_p"]]
ax.barh(y, show["combined_score"].values, color=cols, height=0.62,
        edgecolor=OI["black"], linewidth=0.5)
# annotate the three components to the right of each bar
maxs = show["combined_score"].max()
for yi, (_, r) in zip(y, show.iterrows()):
    ax.text(r.combined_score + maxs * 0.02, yi,
            f"prior={int(r.prior_weight)}  sender \u0394={r.sender_log2FC_sah_vs_sham:+.2f}  "
            f"\u03c1={r.emp_activity_rho:+.2f}"
            + ("*" if (pd.notna(r.emp_activity_p) and r.emp_activity_p < 0.05) else ""),
            va="center", ha="left", fontsize=7.6, color=OI["black"])
ax.set_yticks(y)
ax.set_yticklabels(show["ligand"].values, fontsize=9.5)
ax.set_xlim(0, maxs * 1.95)
ax.set_xlabel("NicheNet-STYLE combined ligand-target score\n(normalised prior \u00d7 sender log2FC \u00d7 empirical activity)", fontsize=8.6)
ax.set_title("Neutrophil ligand \u2192 microglial/macrophage HMOX1\u00b7FTH1\u00b7FTL1\u00b7SLC40A1\n"
             "(NicheNet-STYLE, not NicheNet; curated prior; K=25 UMAP exposure)",
             fontsize=9.8, weight="bold")
ax.axvline(0, color=OI["black"], lw=0.8)
# legend
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=OI["verm"], edgecolor=OI["black"], label="empirical activity p<0.05"),
                   Patch(facecolor=OI["blue"], edgecolor=OI["black"], label="empirical activity n.s.")],
          fontsize=7.8, frameon=False, loc="lower right")
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "Fig25_nichenet_style.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG, "Fig25_nichenet_style.pdf"), bbox_inches="tight", facecolor="white")
plt.close(fig)
print("\nwrote results/nichenet_style_scores.csv + Fig25_nichenet_style.png/pdf")
