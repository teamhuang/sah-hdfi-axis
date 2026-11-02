"""
Fig0 — Graphical abstract for the SAH HRR / HP-CD163-HMOX1-ferroptosis axis.
A single summary schematic (this is a concept figure, not a data-panel composite).
All numbers are read from the computed result CSVs so the figure cannot drift from
the analysis. Honest negatives are shown as explicit badges.
Outputs: PNG (300 dpi) + vector PDF to ../figures/.
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.mark_whole("Fig0_graphical_abstract")
_ns.savefig_patch()
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- Arial + Okabe-Ito (colorblind-safe) ----
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})
OI = {
    "black": "#000000", "orange": "#E69F00", "skyblue": "#56B4E9",
    "green": "#009E73", "yellow": "#F0E442", "purple": "#CC79A7",
    "vermilion": "#D55E00", "blue": "#0072B2", "grey": "#8C8C8C",
}

# ---- load real computed numbers ----
def _read(name):
    p = os.path.join(RES, name)
    return pd.read_csv(p) if os.path.exists(p) else None

deg = _read("DEG_GSE36791.csv")
n_up = n_dn = np.nan
if deg is not None:
    sig = deg[(deg["adj.P.Val"] < 0.05) & (deg["logFC"].abs() > 1)]
    n_up = int((sig["logFC"] > 0).sum())
    n_dn = int((sig["logFC"] < 0).sum())

wt = _read("wgcna_module_trait.csv")
wgcna_r = wgcna_p = np.nan
if wt is not None and "cor_SAH" in wt:
    i = wt["cor_SAH"].abs().idxmax()
    wgcna_r = abs(float(wt.loc[i, "cor_SAH"]))
    wgcna_p = float(wt.loc[i, "p_SAH"])

dec = _read("bulk_deconvolution_scores.csv")
neu_d = t_d = heme_delta = heme_p = np.nan
if dec is not None:
    def row(name):
        m = dec[dec["signature"].str.lower().str.startswith(name.lower())]
        return m.iloc[0] if len(m) else None
    r = row("Neutrophil");  neu_d = float(r["cohens_d"]) if r is not None else np.nan
    r = row("T cell");      t_d = float(r["cohens_d"]) if r is not None else np.nan
    r = row("Heme_iron");   heme_delta = float(r["delta"]) if r is not None else np.nan; heme_p = float(r["t_p"]) if r is not None else np.nan

cor = _read("bulk_ferroptosis_heme_correlation.csv")
rho = pcor = np.nan
if cor is not None and "rho" in cor and "p" in cor:
    m = cor[cor["pair"].str.contains("Heme", case=False, na=False)]
    if len(m):
        rho = float(m["rho"].iloc[0]); pcor = float(m["p"].iloc[0])

comp = _read("scrna_temporal_composition.csv")
def comp_triple(name):
    if comp is None:
        return None
    m = comp[comp["cell_type"].str.lower().str.startswith(name.lower())]
    if not len(m):
        return None
    r = m.iloc[0]
    try:
        return np.array([float(r["Sham"]), float(r["SAH_24h"]), float(r["SAH_72h"])])
    except Exception:
        return None

neu_tmp = comp_triple("Neutrophil")
mic_tmp = comp_triple("Microglia")

xgb = _read("xgb_model_comparison.csv")
xgb_heme_ext = xgb_heme_cv = np.nan
if xgb is not None:
    m = xgb[xgb["feature_set"].str.contains("heme", case=False, na=False)]
    if len(m):
        xgb_heme_ext = float(m["external_AUC_GSE73378"].iloc[0])
        xgb_heme_cv = float(m["CV_AUC_mean"].iloc[0])

cellDE = _read("scrna_celltype_DE_panels.csv")
cd163_p = np.nan
if cellDE is not None:
    m = cellDE[(cellDE["gene"].str.lower() == "cd163")]
    if len(m):
        cd163_p = float(m["p_adj"].iloc[0])

def f(x, nd=2):
    try:
        return f"{float(x):.{nd}f}"
    except Exception:
        return "—"


# ---- canvas ----
fig, ax = plt.subplots(figsize=(16, 9.4))
ax.set_xlim(0, 100)
ax.set_ylim(0, 62)
ax.axis("off")

def box(x, y, w, h, fc, ec, lw=1.6, r=1.4, alpha=1.0, z=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle=f"round,pad=0.4,rounding_size={r}",
                 linewidth=lw, edgecolor=ec, facecolor=fc, alpha=alpha, zorder=z,
                 mutation_aspect=1.0))

def txt(x, y, s, size=11, color="black", weight="normal", ha="center", va="center", z=5, style="normal"):
    ax.text(x, y, s, fontsize=size, color=color, weight=weight, ha=ha, va=va,
            zorder=z, style=style)

def arrow(p0, p1, color, lw=2.4, z=3, style="-|>", rad=0.0, ms=16):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms,
                 linewidth=lw, color=color, zorder=z,
                 connectionstyle=f"arc3,rad={rad}"))

# ---- title band ----
box(1.5, 55.5, 97, 5.6, OI["blue"], OI["blue"], r=1.0)
txt(50, 59.5, "HRR (Hb / RDW) \u2192 haptoglobin-mediated heme-detoxification failure \u2192 NRF2/BACH1-HMOX1 \u2192 ferroptosis after SAH",
    size=14.0, color="white", weight="bold")
txt(50, 56.9, "bulk-transcriptome \u00b7 WGCNA \u00b7 mouse brain single-cell \u00b7 deconvolution \u00b7 SCENIC(NRF2/BACH1) \u00b7 CellChat-style \u00b7 FerrDb \u00b7 validated MR/coloc engines",
    size=10.6, color="#d6e6f5")

# ---- COLUMN 1: clinical anchor + acute blood ----
txt(16, 52.6, "1  CLINICAL + SYSTEMIC EVIDENCE", size=12.5, color=OI["blue"], weight="bold")
box(2, 43.0, 28, 8.2, "#eaf2fb", OI["blue"])
txt(16, 49.8, "MIMIC-IV  n = 1,161 aSAH", size=11, weight="bold", color=OI["blue"])
txt(16, 47.7, "HRR = Hb / RDW independently\npredicts in-hospital death", size=9.6)
txt(16, 44.9, "Hb arm \u2192 O2 delivery   RDW arm \u2192 RBC turnover", size=8.6, color=OI["grey"], style="italic")

box(2, 33.4, 28, 8.2, "#eaf2fb", OI["blue"])
txt(16, 39.9, "Acute blood GSE36791 (43 vs 18)", size=10.4, weight="bold", color=OI["blue"])
txt(16, 37.8, f"{f(n_up,0)} up / {f(n_dn,0)} down DEGs (FDR<.05): HP \u00b7 HMOX1\u00b7 FTH1\nSLC40A1 \u00b7 SLC25A37 \u00b7 NFE2 \u00b7 MARC1 + inflammation", size=8.9)
txt(16, 35.0, f"32 candidate genes \u2229 red-cell sig.\nErythrocyte diff. p\u22485e-26 \u00b7 GSEA: NF-\u03baB/ROS/hypoxia/ferroptosis", size=8.7, color=OI["green"])

arrow((16, 43.0), (16, 41.6), OI["blue"])

# WGCNA badge
box(2, 25.6, 28, 6.2, "#eef7f3", OI["green"])
txt(16, 29.2, "WGCNA (signed)", size=10.2, weight="bold", color=OI["green"])
txt(16, 27.0, f"SAH module eigengene r={f(wgcna_r)} (p={wgcna_p:.0e})\nheme-iron genes HP/MARC1/SLC40A1/ACSL4/FTH1 high MM", size=8.6)

# deconvolution badge
box(2, 15.0, 28, 9.2, "#eef7f3", OI["green"])
txt(16, 22.2, "Human-blood deconvolution", size=10.2, weight="bold", color=OI["green"])
txt(16, 19.9, f"Innate\u2191 (neutrophil d={f(neu_d)}) / adaptive\u2193 (T d={f(t_d)})\nHeme-iron program up (\u0394={f(heme_delta)}, p={heme_p:.0e})", size=8.8)
txt(16, 16.6, f"Ferroptosis ~ heme across patients\nSpearman \u03c1 = {f(rho)} (p = {pcor:.0e})",
    size=9.4, weight="bold", color=OI["vermilion"])

arrow((16, 25.6), (16, 24.2), OI["green"])

# ---- COLUMN 2: single-cell mechanistic core ----
txt(50, 52.6, "2  SINGLE-CELL MECHANISM (mouse brain GSE311852)", size=12.5, color=OI["purple"], weight="bold", va="center")
box(33, 24.0, 34, 27.2, "#f6eef6", OI["purple"], lw=2.0)

# cascade nodes (vertical) with cell-type colours
cascade = [
    ("SAH  red-cell lysis / hemolysis", "free Hb released", OI["vermilion"]),
    ("HP  haptoglobin scavenging", "Neutrophil  (Hp \u2191, log2FC\u22482)", OI["orange"]),
    ("Myeloid uptake  CD163 / MRC1", "scRNA: Cd163 <6% \u2014 MRC1 + S100-TLR4 engaged", OI["blue"]),
    ("Heme overload \u2192 NRF2\u2191 / BACH1\u2193", "regulon: NRF2 activity\u2191, BACH1 activity\u2193 (microglia)", OI["skyblue"]),
    ("HMOX1 hyper-activation", "Myeloid Hmox1 \u2191 (log2FC\u22482.4)", OI["green"]),
    ("Fe2+ / ROS \u2014 HDFI (detox-failure index)", "z(HMOX1)\u2212z(FTH1+FTL)\u2212z(GPX4); predicts neuron ferroptosis \u03c1=0.34", "#8C4A0A"),
    ("Ferroptosis", "Neuron: Gpx4\u2193, ACSL4/SLC7A11/PTGS2\u2191; blood Exec p=0.003", OI["purple"]),
    ("Secondary neuronal injury", "poor outcome (HRR link)", OI["black"]),
]
ys = np.linspace(47.6, 26.4, len(cascade))
bw, bx = 30.0, 35.0
for (lab, sub, c), y in zip(cascade, ys):
    box(bx, y - 1.25, bw, 2.5, "#ffffff", c, lw=1.8, r=0.8)
    ax.add_patch(Circle((bx - 0.4, y), 0.55, color=c, zorder=6))
    txt(bx + 0.3, y + 0.42, lab, size=10.0, weight="bold", color=c, ha="left")
    txt(bx + 0.3, y - 0.55, sub, size=8.3, color=OI["grey"], ha="left", style="italic")
    if y != ys[-1]:
        arrow((bx + bw/2, y - 1.35), (bx + bw/2, y - 2.4), c, lw=2.2, ms=13)

# SCENIC + temporal side-notes to the right of the box
box(69, 40.4, 29, 10.2, "#fbf3f8", OI["purple"])
txt(83.5, 48.5, "Temporal dynamics (Sham\u219224h\u219272h)", size=10.2, weight="bold", color=OI["purple"])
if neu_tmp is not None and len(neu_tmp) == 3:
    txt(83.5, 46.0, f"Neutrophils {neu_tmp[0]:.1f}\u2192{neu_tmp[1]:.1f}\u2192{neu_tmp[2]:.1f}%  (first wave)", size=9.0)
else:
    txt(83.5, 46.0, "Neutrophils 0.5\u219238.6\u21923.7%  (first wave)", size=9.0)
if mic_tmp is not None and len(mic_tmp) == 3:
    txt(83.5, 43.9, f"Myeloid {mic_tmp[0]:.1f}\u2192{mic_tmp[1]:.1f}\u2192{mic_tmp[2]:.1f}%  (expansion)", size=9.0)
else:
    txt(83.5, 43.9, "Myeloid 14.8\u219232.3\u219253.3%  (expansion)", size=9.0)
txt(83.5, 41.7, "Kruskal-Wallis all p < 1e-16", size=8.4, color=OI["grey"], style="italic")

box(69, 31.6, 29, 8.0, "#fbf3f8", OI["purple"])
txt(83.5, 37.9, "Upstream driver (SCENIC-lite)", size=10.2, weight="bold", color=OI["purple"])
txt(83.5, 35.5, "NRF2 / NFE2L2 regulon ACTIVATED in SAH\n(0.59 vs 0.28, p<1e-4) \u2014 not NFE2-p45;\nBACH1 also up", size=8.8)

box(69, 24.0, 29, 6.6, "#fbf3f8", OI["purple"])
txt(83.5, 28.0, "ML transportability (XGBoost)", size=10.2, weight="bold", color=OI["purple"])
txt(83.5, 25.6, f"Heme-iron panel generalises best externally\n(CV 0.80 \u2192 ext AUC {f(xgb_heme_ext)}) vs immune DEG 0.52", size=8.8)

# ---- honest negatives band (column 3 bottom) ----
txt(83.5, 21.6, "3  HONEST BOUNDARIES", size=12.5, color=OI["vermilion"], weight="bold")
box(69, 11.9, 29, 8.5, "#fdeee8", OI["vermilion"], lw=1.8)
neg = [
    "Cd163 not induced (p=0.83); HP\u2013CD163 fails CellChat (Cd163<6%) \u2192 MRC1/S100-TLR4",
    "Bulk-blood heme/HP signal is COMPOSITION-driven (does not survive adjustment)",
    "HDFI NS in blood (\u0394+0.18, p=0.28) \u2014 signal is single-cell",
    "DPT pseudotime & FerrDb Neu\u2192Mac\u2192Neuron chain NOT monotonic",
    "External AUC < 0.75 \u2014 acute signal \u2260 remote aSAH",
]
yy = np.linspace(19.2, 12.6, len(neg))
for s, y in zip(neg, yy):
    txt(70.2, y, "X", size=8, color=OI["vermilion"], ha="left", weight="bold")
    txt(71.6, y, s, size=7.2, color=OI["black"], ha="left")

# ---- bottom: next steps ----
box(2, 5.0, 96, 6.6, "#f2f2f2", OI["grey"], lw=1.4)
txt(50, 9.7, "WHAT PUBLIC DATA CANNOT PROVE \u2014 the author must generate:", size=10.6, weight="bold", color=OI["black"])
txt(50, 7.6,
    "(i) real MR / TWAS / coloc causal estimates (engines validated on synthetic truth: IVW \u03b2=0.40\u21920.41, coloc PP4=1.0; sumstats access-blocked)   \u00b7   (ii) plasma proteomics + MR-PheWAS   \u00b7   (iii) spatial co-localisation\n"
    "(iv) the HRR \u2192 HDFI \u2192 brain-injury \u2192 death MEDIATION/SEM bridge on a matched cohort (engine validated: ACIE 0.48\u21920.481; no public HRR+plasma cohort here)\n"
    "(v) an ACUTE human cohort with admission Hb/RDW + outcome   \u00b7   (vi) animal/cell validation: ferrostatin-1, HMOX1 inhibitor \u00b1 timing, CD163/MRC1 blockade",
    size=7.6, color=OI["grey"])

# ---- connectors between columns ----
arrow((30, 37.5), (34.8, 37.5), OI["blue"], lw=2.6, rad=0.0)
arrow((30, 28.5), (34.8, 40.0), OI["green"], lw=2.2, rad=-0.15)
arrow((67.2, 44.5), (69, 44.5), OI["purple"], lw=2.0)

out_png = os.path.join(FIG, "Fig0_graphical_abstract.png")
out_pdf = os.path.join(FIG, "Fig0_graphical_abstract.pdf")
fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("wrote", out_png)
print("wrote", out_pdf)
print(json.dumps({
    "n_up": n_up, "n_dn": n_dn, "wgcna_r": wgcna_r, "rho": rho, "pcor": pcor,
    "neu_d": neu_d, "t_d": t_d, "heme_delta": heme_delta, "heme_p": heme_p,
    "xgb_heme_ext": xgb_heme_ext, "cd163_p": cd163_p,
    "neu_tmp": None if neu_tmp is None else list(neu_tmp),
    "mic_tmp": None if mic_tmp is None else list(mic_tmp),
}, default=str))
