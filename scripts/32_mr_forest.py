"""32_mr_forest.py -- publication forest plot for the REAL two-sample MR
(genetically higher plasma Haptoglobin -> subarachnoid haemorrhage risk).
Reads results/mr_real_HP_SAH_results.csv + results/gwas/MR_step2_instruments_final.csv.
Okabe-Ito colour-safe palette, Arial, no panel letters, PNG(300dpi)+vector PDF.
Nothing is re-computed here beyond what 31_mr_real.py produced from real sumstats."""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES  = os.path.join(BASE, "results")
GW   = os.path.join(RES, "gwas")
FIG  = os.path.join(BASE, "figures")
OI   = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
        "verm": "#D55E00", "purple": "#CC79A7", "sky": "#56B4E9", "grey": "#555555"}
rcParams.update({"font.family": "Arial", "font.size": 9, "axes.linewidth": 0.9,
                 "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.spines.top": False,
                 "axes.spines.right": False})

def main():
    summ = pd.read_csv(os.path.join(RES, "mr_real_HP_SAH_results.csv"))
    iv   = pd.read_csv(os.path.join(GW, "MR_step2_instruments_final.csv"))
    s = summ.iloc[0]
    iv = iv.copy()
    iv["ratio"] = iv.bo / iv.be
    iv["ratio_se"] = iv.se_o / iv.be.abs()
    iv = iv.sort_values("ratio")
    n = len(iv)
    fig, ax = plt.subplots(figsize=(5.6, 0.34*n + 1.7))
    y = np.arange(n, 0, -1)
    for yi, (_, r) in zip(y, iv.iterrows()):
        lo, hi = r.ratio - 1.96*r.ratio_se, r.ratio + 1.96*r.ratio_se
        ax.plot([lo, hi], [yi, yi], color=OI["grey"], lw=1.0, zorder=1)
        ax.plot([r.ratio], [yi], marker="s", ms=5, color=OI["blue"], zorder=2)
    ax.axvline(0, color="black", lw=0.8, ls="-", zorder=3)
    # method summary bands
    methods = []
    for name, bcol, scol, pcol, col in [
        ("IVW (fixed)", "IVW_beta", "IVW_se", "IVW_p", OI["verm"]),
        ("IVW (random)", "IVW_beta", "IVW_se_random", "IVW_p_random", OI["orange"]),
        ("MR-Egger", "Egger_beta", "Egger_se", "Egger_p", OI["purple"]),
        ("Weighted median", "WMedian_beta", "WMedian_se", "WMedian_p", OI["green"]),
    ]:
        if bcol in s and not pd.isna(s.get(bcol)) and s.get(scol) and not pd.isna(s.get(scol)):
            methods.append((name, float(s[bcol]), float(s[scol]), col))
    ytop = n + 1.6
    for i, (name, b, se, col) in enumerate(methods):
        yi = ytop + i*1.15
        lo, hi = b - 1.96*se, b + 1.96*se
        ax.plot([lo, hi], [yi, yi], color=col, lw=1.6, zorder=4)
        ax.plot([b], [yi], marker="D", ms=6.5, color=col, zorder=5)
    # axis / labels
    ax.set_yticks(list(y) + [ytop + i*1.15 for i in range(len(methods))])
    ax.set_yticklabels([f"{r['rsid']}" for _, r in iv.iterrows()]
                       + [m[0] for m in methods], fontsize=7.5)
    allx = list(iv.ratio - 1.96*iv.ratio_se) + list(iv.ratio + 1.96*iv.ratio_se)
    for _, b, se, _ in methods:
        allx += [b-1.96*se, b+1.96*se]
    ax.set_xlim(min(allx)-0.05, max(allx)+0.05)
    ax.set_ylim(0.3, ytop + len(methods)*1.15 + 0.6)
    ax.set_xlabel("Per-allele log-odds of SAH per SD genetically higher plasma Haptoglobin", fontsize=8.5)
    ax.set_title("MR: plasma Haptoglobin -> subarachnoid haemorrhage (real GWAS sumstats)",
                 fontsize=9.5, pad=8)
    note = (f"IVW beta={s['IVW_beta']:.3f} (95% CI {s['IVW_beta']-1.96*s['IVW_se']:.3f}, "
            f"{s['IVW_beta']+1.96*s['IVW_se']:.3f}), p={s['IVW_p']:.2g}; "
            f"Egger intercept p={s.get('Egger_intercept_p', float('nan')):.2g}; "
            f"Q={s['IVW_Q']:.1f} (df={int(s['IVW_df'])}); {int(s['n_instruments'])} IVs, mean F={s['F_stat_mean']:.0f}")
    ax.text(0.0, -0.16, note, transform=ax.transAxes, fontsize=6.6, color=OI["grey"], wrap=True)
    ax.grid(axis="x", ls=":", lw=0.5, color="#cccccc", zorder=0)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    png = os.path.join(FIG, "Fig22_mr_forest.png"); pdf = os.path.join(FIG, "Fig22_mr_forest.pdf")
    fig.savefig(png, dpi=300); fig.savefig(pdf)
    print("saved", png); print("saved", pdf)

if __name__ == "__main__":
    main()
