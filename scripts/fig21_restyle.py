"""fig21_restyle.py -- offline re-render of Fig21 (regulatory motif) from cached CSV.

Script 30 re-queries Ensembl live and hangs here; the results are already cached in
results/regulatory_motif_HMOX1.csv, so this redraws the three panels in the shared
Nature style (naturestyle splits them into Fig21_regulatory_motif_a/_b/_c).
"""
import os, ast, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import naturestyle as ns
ns.savefig_patch()

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES, FIG = os.path.join(BASE, "results"), os.path.join(BASE, "figures")
d = pd.read_csv(os.path.join(RES, "regulatory_motif_HMOX1.csv"))
genes = d.gene.tolist()
C = ns.OKABE

fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
# a: ARE consensus sites
ax = axes[0]
cols = [C["vermil"] if g == "HMOX1" else C["blue"] for g in genes]
ax.bar(range(len(genes)), d.are_consensus, color=cols)
ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=60, ha="right", fontsize=7)
ax.set_ylabel("ARE consensus sites (promoter)")
# b: NFE2L2 PWM best score
ax = axes[1]
ax.bar(range(len(genes)), d.nfe2l2_best, color=C["green"])
ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=60, ha="right", fontsize=7)
ax.set_ylabel("NFE2L2 (JASPAR) best PWM score")
# c: ENCODE regulatory features (promoter / enhancer / CTCF)
ax = axes[2]
rt = d.reg_types.apply(lambda s: ast.literal_eval(s) if isinstance(s, str) and s.strip() else {})
prom = rt.apply(lambda x: x.get("promoter", 0))
enh = rt.apply(lambda x: x.get("enhancer", 0))
ctcf = rt.apply(lambda x: x.get("CTCF_binding_site", 0))
x = np.arange(len(genes)); w = 0.27
ax.bar(x - w, prom, w, color=C["blue"], label="promoter")
ax.bar(x, enh, w, color=C["orange"], label="enhancer")
ax.bar(x + w, ctcf, w, color=C["purple"], label="CTCF")
ax.set_xticks(x); ax.set_xticklabels(genes, rotation=60, ha="right", fontsize=7)
ax.set_ylabel("ENCODE regulatory features"); ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "Fig21_regulatory_motif.png"), dpi=300)
fig.savefig(os.path.join(FIG, "Fig21_regulatory_motif.pdf"))
print("re-rendered Fig21 from cache")
