"""52_scrna_reference_panels.py -- Nature-styled single-cell reference panels.

Re-renders the three raw scanpy outputs (UMAP cell type, heme/iron UMAPs,
marker dotplot) without any in-figure titles/notes, as independent panels
Supplementary Figure S5 a/b/c.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import naturestyle as ns

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ns.here = BASE
import anndata as ad

H5 = os.path.join(BASE, "_processed", "sah_scrna_annot.h5ad")
OUT = os.path.join(BASE, "figures")

HEME_GENES = ["Hmox1", "Fth1", "Ftl1", "Slc40a1", "Nfe2", "Hprt"]
MARKERS = {
    "Astrocyte": ["Gja1", "Slc1a3"], "B cell": ["Cd79a", "Ms4a1"],
    "Endothelial": ["Cldn5", "Pecam1"], "Erythroid/RBC": ["Hba-a1", "Hbb-bs"],
    "Mast": ["Cpa3", "Kit"], "Microglia/Mac": ["Cx3cr1", "P2ry12"],
    "NK/GM": ["Klrd1", "Gzma"], "Neuron": ["Snap25", "Syt1"],
    "Neutrophil": ["Elane", "Camp"], "OPC": ["Pdgfra", "Sox10"],
    "Oligodendrocyte": ["Mbp", "Plp1"], "T cell": ["Cd3e", "Cd8a"],
}


def main():
    a = ad.read_h5ad(H5)
    rng = np.random.default_rng(0)
    order = rng.permutation(a.n_obs)

    # --- panel a: UMAP coloured by cell type --------------------------------
    ct = a.obs["cell_type"]
    types = list(ct.cat.categories)
    codes = ct.cat.codes.values
    cmap = matplotlib.colors.ListedColormap([
        ns.OKABE["blue"], ns.OKABE["vermil"], ns.OKABE["green"], ns.OKABE["orange"],
        ns.OKABE["purple"], ns.OKABE["sky"], "#8C564B", "#E79F9F", "#9AA531",
        "#B07AA1", "#BAB0AC", "#4C72B0"][: len(types)])
    codes = a.obs["cell_type"].cat.codes.values
    fig, ax = plt.subplots(figsize=(5 / 2.54, 5 / 2.54))
    sc = ax.scatter(a.obsm["X_umap"][order, 0], a.obsm["X_umap"][order, 1],
                    c=codes[order], cmap=cmap, s=0.6, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=3.5,
                          color=cmap(i / max(1, len(types) - 1))) for i in range(len(types))]
    ax.legend(handles, types, fontsize=5.5, loc="center left",
              bbox_to_anchor=(1.02, 0.5), frameon=False, handletextpad=0.3)
    ns.save(fig, os.path.join(OUT, "umap_celltype"), width_cm=6.5, strip_notes=False)

    # --- panel b: heme/iron module score + component genes ------------------
    import scanpy as sc_mod
    sc_mod.tl.score_genes(a, HEME_GENES, score_name="heme_score", use_raw=False)
    fields = [("heme_score", None), ("Hp", None), ("Cd163", None), ("Hmox1", None)]
    fig, axes = plt.subplots(2, 2, figsize=(6.5 / 2.54, 6.0 / 2.54))
    for (name, _), ax in zip(fields, axes.ravel()):
        v = (a.obs[name].values if name == "heme_score"
             else np.asarray(a[:, name].X.todense()).ravel())
        pts = ax.scatter(a.obsm["X_umap"][order, 0], a.obsm["X_umap"][order, 1],
                         c=v[order], cmap=ns.SEQ, s=0.6, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        cb = fig.colorbar(pts, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=5.5)
        cb.outline.set_linewidth(0.5)
    fig.tight_layout(pad=0.3)
    ns.save(fig, os.path.join(OUT, "umap_heme"), width_cm=5, strip_notes=False)

    # --- panel c: marker dotplot -------------------------------------------
    genes = sum(MARKERS.values(), [])
    labels = list(MARKERS.keys())
    X = a[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    ct = a.obs["cell_type"].astype(str).values
    mean = np.zeros((len(labels), len(genes)))
    frac = np.zeros_like(mean)
    for i, lab in enumerate(labels):
        m = ct == lab
        if m.sum():
            mean[i] = X[m].mean(0)
            frac[i] = (X[m] > 0).mean(0)
    fig, ax = plt.subplots(figsize=(6.5 / 2.54, 4.6 / 2.54))
    xs, ys = np.meshgrid(np.arange(len(genes)), np.arange(len(labels)))
    sz = 260 * frac
    im = ax.scatter(xs.ravel(), ys.ravel(), s=sz.ravel(), c=mean.ravel(),
                    cmap=ns.DIV, vmin=-np.percentile(np.abs(mean), 98),
                    vmax=np.percentile(np.abs(mean), 98), linewidths=0)
    ax.set_xticks(np.arange(len(genes)), genes, rotation=90, fontsize=5.5)
    ax.set_yticks(np.arange(len(labels)), labels, fontsize=5.5)
    ax.set_xlim(-0.6, len(genes) - 0.4); ax.set_ylim(-0.6, len(labels) - 0.4)
    ax.invert_yaxis()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=2)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_linewidth(0.5)
    fig.tight_layout(pad=0.3)
    ns.save(fig, os.path.join(OUT, "dotplot__markers"), width_cm=6.5, strip_notes=False)
    print("re-rendered S5 reference panels (a/b/c), no in-figure titles")


if __name__ == "__main__":
    main()
