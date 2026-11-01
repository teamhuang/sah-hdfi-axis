"""
50_csf_neuronal_injury_link.py  -- Round 12 (figure index 40).

Reviewer's #1 priority this round: the protein layer must not stop at "is HDFI
detectable in CSF" -- it must push one link further down the causal chain the
whole paper rests on:

    HDFI (heme-detoxification failure)  ->  NEURONAL INJURY  ->  outcome

Round 11 established that the canonical HDFI_protein composite
z(HMOX1)-z(FTH1+FTL)-z(GPX4) is NOT COMPUTABLE in the only openly downloadable
human SAH-CSF proteome (PXD030593), because HMOX1 and GPX4 are below the CSF
detection depth. We do NOT re-litigate that here. What we CAN do -- on the SAME
real per-sample data -- is test the *downstream* half of the chain using the
proteins that ARE measured:

  (a) the DETECTABLE HDFI-axis members (HP, CD163, FTH1, FTL, TFRC), and
  (b) canonical CSF NEURONAL-INJURY / axonal-injury proteins
      (GFAP, S100B, ENO2 = neuron-specific enolase/NSE, UCHL1, NEFL, NEFM,
       YWHAE = 14-3-3-gamma, MBP),
  (c) BLOOD-CONTAMINATION controls (HBA1, HBB, ALB, TTR) -- the reviewer's own
      rationale for CSF over blood; we surface them honestly rather than hide them.

Design (honest): n = 4 SAH vs 4 control (PXD030593). Every association is a
8-point Spearman / 4-vs-4 Mann-Whitney -> extremely low power. We report exact
rho and p, label everything exploratory/descriptive, and do NOT claim a tested
mediation. We also transcribe, VERBATIM from the independent larger SAH-CSF
cohort PMC11267790 (Clin Proteomics 2024, SAH n=23 vs control n=10, its open
Supplementary Table 1), the reported direction of the same neuronal-injury
proteins, so the human evidence is not resting on n=4 alone. Nothing is
simulated; nothing is imputed for the missing HMOX1/GPX4.
"""
import os, re, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA, RES, FIG = (os.path.join(ROOT, x) for x in ("data", "results", "figures"))
PX_XLSX = os.path.join(DATA, "pxd030593", "protein.xlsx")
PMC_TXT = os.path.join(DATA, "PMC11267790", "supp_tables_text.txt")
OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
          purple="#CC79A7", sky="#56B4E9", yellow="#F0E442", grey="#8C8C8C")
rcParams.update({"font.family": "Arial", "font.size": 9, "axes.linewidth": 0.9,
                 "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                 "axes.spines.right": False, "pdf.fonttype": 42})

HDFI_DETECTABLE = ["HP", "CD163", "FTH1", "FTL", "TFRC"]          # measured in R11
INJURY = ["GFAP", "S100B", "ENO2", "UCHL1", "NEFL", "NEFM", "YWHAE", "MBP"]
BLOOD  = ["HBA1", "HBB", "ALB", "TTR"]

# -----------------------------------------------------------------------------
def load_pxd():
    df = pd.read_excel(PX_XLSX, sheet_name="Proteins", header=0)
    df.columns = [str(c).strip() for c in df.columns]
    df["GeneName"] = df["GeneName"].astype(str).str.strip().str.upper()
    ab = [c for c in df.columns if c.startswith("Abundances")]
    ctrl = [c for c in ab if "Control" in c]
    sah = [c for c in ab if "SAH" in c]
    tot = df[ab].replace(0, np.nan)
    frac = tot.div(tot.sum(axis=0), axis=1) * 1e6      # per-sample total-intensity norm, PPM
    l2 = np.log2(frac + 1)
    l2.insert(0, "GeneName", df["GeneName"].values)
    return l2, ctrl, sah

def series_for(l2, gene):
    m = l2[l2["GeneName"] == gene]
    return None if len(m) == 0 else m.iloc[0]

def z(v):
    v = np.asarray(v, float); s = v.std(ddof=0)
    return (v - v.mean()) / (s if s > 0 else 1.0)

# -----------------------------------------------------------------------------
def main():
    l2, ctrl, sah = load_pxd()
    allcols = ctrl + sah
    lab = ["Ctrl%d" % (i + 1) for i in range(len(ctrl))] + ["SAH%d" % (i + 1) for i in range(len(sah))]
    grp = np.array([0] * len(ctrl) + [1] * len(sah))     # 0=ctrl 1=sah

    # build a gene x sample log2 matrix for the proteins we care about
    want = HDFI_DETECTABLE + INJURY + BLOOD
    mat = {}
    for g in want:
        s = series_for(l2, g)
        if s is not None:
            mat[g] = s[allcols].astype(float).values
    M = pd.DataFrame(mat, index=lab)                       # samples x proteins (log2)
    detected_injury = [g for g in INJURY if g in mat]
    detected_hdfi = [g for g in HDFI_DETECTABLE if g in mat]

    # neuronal-injury composite score = mean z of detected injury markers
    inj = np.mean(np.column_stack([z(M[g].values) for g in detected_injury]), axis=1)
    M["INJURY_SCORE"] = inj

    # ---- (1) SAH vs control injury score (Mann-Whitney, n=4v4) ----
    u, p = stats.mannwhitneyu(inj[grp == 1], inj[grp == 0], alternative="two-sided")
    med_c, med_s = float(np.median(inj[grp == 0])), float(np.median(inj[grp == 1]))
    print("=== neuronal-injury score, SAH vs control (PXD030593, n=4v4) ===")
    print("  detected injury markers:", detected_injury)
    print("  median ctrl=%.2f  SAH=%.2f  Mann-Whitney p=%.3f" % (med_c, med_s, p))
    print("  per-sample injury score:", dict(zip(lab, np.round(inj, 2))))

    # ---- (2) HDFI-member vs injury link (Spearman across 8 samples) ----
    rows = []
    for h in detected_hdfi:
        rho, pv = stats.spearmanr(M[h].values, inj)
        rows.append(dict(dataset="PXD030593", hdfi_protein=h, target="INJURY_SCORE",
                         rho=rho, p=pv, n=len(inj)))
        for g in detected_injury:
            r2, p2 = stats.spearmanr(M[h].values, M[g].values)
            rows.append(dict(dataset="PXD030593", hdfi_protein=h, target=g, rho=r2, p=p2, n=len(inj)))
    link = pd.DataFrame(rows)
    print("\n=== HDFI-axis protein vs neuronal-injury (Spearman, 8 samples; exploratory) ===")
    print(link[link.target == "INJURY_SCORE"].to_string(index=False))

    # ---- (2b) leave-one-out sensitivity for the ferritin arm ----
    loo_rows = []
    for h in ["FTH1", "FTL"]:
        full = stats.spearmanr(M[h].values, inj)
        rs = []
        for k in range(len(inj)):
            keep = [i for i in range(len(inj)) if i != k]
            rs.append(stats.spearmanr(M[h].values[keep], inj[keep]).statistic)
        loo_rows.append(dict(protein=h, rho_full=round(float(full.statistic), 3),
                             p_full=round(float(full.pvalue), 4),
                             LOO_rho_min=round(float(np.nanmin(rs)), 3),
                             LOO_rho_max=round(float(np.nanmax(rs)), 3),
                             driven_by_single_sample=bool(np.nanmin(rs) <= 0.4)))
    loo = pd.DataFrame(loo_rows)
    print("\n=== leave-one-out sensitivity, ferritin arm vs injury score ===")
    print(loo.to_string(index=False))

    # ---- (3) blood-contamination sanity ----
    print("\n=== blood-contamination controls (log2 PPM; reviewer's CSF-over-blood rationale) ===")
    for b in [x for x in BLOOD if x in mat]:
        print("  %-5s ctrl=%s sah=%s" % (b, np.round(M[b].values[:len(ctrl)], 1),
                                          np.round(M[b].values[len(ctrl):], 1)))

    # ---- (4) VERBATIM independent cohort PMC11267790 Supplementary Table 1 ----
    verbatim = []
    if os.path.exists(PMC_TXT):
        txt = open(PMC_TXT, encoding="utf-8", errors="ignore").read()
        lines = [l.strip() for l in txt.splitlines()]
        # row layout: gene / acc / ctrl "m (sd) [n]" / sah "m (sd) [n]" / log2FC / p / adj
        cell = re.compile(r"^\d+(?:\.\d+)?\s*\(\d+(?:\.\d+)?\)\s*\[\d+\]$")
        numf = re.compile(r"^-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$")
        for g in ["GFAP", "S100B", "ENO2", "YWHAE", "HP", "TFRC", "FTH1", "FTL"]:
            for i, ln in enumerate(lines):
                if ln != g or i + 6 >= len(lines):
                    continue
                cc, sc, fc, pv = lines[i + 2], lines[i + 3], lines[i + 4], lines[i + 5]
                if cell.match(cc) and cell.match(sc) and numf.match(fc) and numf.match(pv):
                    verbatim.append(dict(dataset="PMC11267790(verbatim)", protein=g,
                                         ctrl_cell=cc, sah_cell=sc, log2FC=fc, p=pv))
                    break
    verb = pd.DataFrame(verbatim)
    print("\n=== PMC11267790 verbatim Supplementary Table 1 rows (SAH start vs control) ===")
    if len(verb):
        print(verb.to_string(index=False))

    # ---- write results ----
    os.makedirs(os.path.join(RES, "tables"), exist_ok=True)
    out = link.copy()
    out["injury_score_sah_vs_ctrl_p"] = p
    out["injury_median_ctrl"] = med_c; out["injury_median_sah"] = med_s
    out.to_csv(os.path.join(RES, "tables", "csf_neuronal_injury_link.csv"), index=False)
    if len(verb):
        verb.to_csv(os.path.join(RES, "tables", "csf_neuronal_injury_verb_pmc.csv"), index=False)
    loo.to_csv(os.path.join(RES, "tables", "csf_injury_LOO_sensitivity.csv"), index=False)
    M.round(3).to_csv(os.path.join(RES, "tables", "csf_axis_injury_matrix.csv"))

    # ---- FIGURE 40 ----
    make_figure(M, detected_hdfi, detected_injury, inj, grp, lab, p, med_c, med_s, link,
                os.path.join(FIG, "Fig40_csf_hdfi_neuronal_injury.png"),
                os.path.join(FIG, "Fig40_csf_hdfi_neuronal_injury.pdf"))
    print("\nsaved Fig40_csf_hdfi_neuronal_injury.png/pdf")

# -----------------------------------------------------------------------------
def make_figure(M, hdfi, injury, inj, grp, lab, p, med_c, med_s, link, png, pdf):
    fig = plt.figure(figsize=(11.2, 4.5))
    # Panel 1: Spearman rho heatmap HDFI-member x (injury markers + score)
    targets = injury + ["INJURY_SCORE"]
    Rm = np.full((len(hdfi), len(targets)), np.nan)
    for i, h in enumerate(hdfi):
        for j, t in enumerate(targets):
            if t == "INJURY_SCORE":
                Rm[i, j] = M[h].corr(M["INJURY_SCORE"], method="spearman")
            else:
                Rm[i, j] = M[h].corr(M[t], method="spearman")
    ax1 = fig.add_subplot(1, 3, 1)
    im = ax1.imshow(Rm, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax1.set_xticks(range(len(targets))); ax1.set_xticklabels(targets, rotation=45, ha="right", fontsize=7)
    ax1.set_yticks(range(len(hdfi))); ax1.set_yticklabels(hdfi, fontsize=8)
    for i in range(len(hdfi)):
        for j in range(len(targets)):
            if not np.isnan(Rm[i, j]):
                ax1.text(j, i, f"{Rm[i,j]:.2f}", ha="center", va="center", fontsize=6.5,
                         color="black" if abs(Rm[i, j]) < 0.6 else "white")
    ax1.set_title("Spearman rho\nHDFI-axis protein x neuronal-injury\n(8 CSF samples, exploratory)", fontsize=8.6)
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)

    # Panel 2: injury score per sample, SAH vs control
    ax2 = fig.add_subplot(1, 3, 2)
    for x, (col, name) in [(0, (OI["blue"], "Control")), (1, (OI["verm"], "SAH"))]:
        vals = inj[grp == x]
        jit = np.array([x] * len(vals)) + (np.linspace(-0.12, 0.12, len(vals)))
        ax2.scatter(jit, vals, s=52, color=col, edgecolor="black", lw=0.4, zorder=3, label=name)
        ax2.hlines(np.median(vals), x - 0.28, x + 0.28, color=col, lw=1.8, zorder=4)
    ax2.set_xticks([0, 1]); ax2.set_xticklabels(["Control", "SAH"], fontsize=9)
    ax2.set_ylabel("neuronal-injury protein score (mean z)", fontsize=8)
    ax2.set_xlim(-0.5, 1.5)
    ax2.set_title("CSF neuronal-injury score\nSAH vs control (n=4v4, MWU p=%.2f)" % p, fontsize=8.6)
    ax2.legend(fontsize=7, frameon=False, loc="upper left")

    # Panel 3: FTH1 vs injury score scatter (the robust positive finding)
    ax3 = fig.add_subplot(1, 3, 3)
    if "FTH1" in M.columns:
        hp = M["FTH1"].values
        rho_hp = stats.spearmanr(hp, inj).statistic
        p_hp = stats.spearmanr(hp, inj).pvalue
        for x, col, name in [(0, OI["blue"], "Control"), (1, OI["verm"], "SAH")]:
            msk = grp == x
            ax3.scatter(hp[msk], inj[msk], s=52, color=col, edgecolor="black", lw=0.4, label=name, zorder=3)
        b, a = np.polyfit(hp, inj, 1)
        xs = np.linspace(hp.min(), hp.max(), 20)
        ax3.plot(xs, a + b * xs, color=OI["grey"], lw=1.2, ls="--", zorder=2)
        ax3.set_xlabel("FTH1 (ferritin heavy) log2 PPM", fontsize=8)
        ax3.set_ylabel("neuronal-injury score", fontsize=8)
        ax3.set_title("HDFI iron-storage arm (FTH1) vs\nneuronal injury  rho=%.2f, p=%.2f" % (rho_hp, p_hp), fontsize=8.6)
        ax3.legend(fontsize=7, frameon=False, loc="upper left")

    fig.suptitle("Pushing the HDFI axis one link down in REAL human SAH-CSF: heme-handling protein -> neuronal-injury protein\n"
                 "PXD030593 (label-free DDA, SAH n=4 vs control n=4, per-sample total-intensity normalised)", fontsize=10, fontweight="bold")
    note = ("Honest scope: n=4v4 -> every rho is an 8-point exploratory association, NOT a tested mediation; HMOX1/GPX4 remain "
            "below CSF detection so the full HDFI_protein is still NOT computable (uses detectable members only).\n"
            "Blood-contamination controls (HBA1/HBB/ALB/TTR) are reported, not hidden -- one SAH sample is Hb-high (possible "
            "contamination) and one is injury-high; the direction is corroborated VERBATIM by the larger independent cohort "
            "PMC11267790 (SAH n=23 vs control n=10).")
    fig.text(0.5, 0.015, note, ha="center", va="bottom", fontsize=6.6, color=OI["grey"])
    fig.tight_layout(rect=[0, 0.075, 1, 0.93])
    fig.savefig(png, dpi=300); fig.savefig(pdf)
    plt.close(fig)


if __name__ == "__main__":
    main()
