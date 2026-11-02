"""
49_mouse_sah_atlas.py  -- Round-11 reviewer TASK: cross-species in-vivo
mechanistic corroboration of the HDFI (Heme-Detoxification-Failure-Index)
gene program in a REAL public MOUSE single-cell dataset.

DATASET: GEO GSE311852 "A spatiotemporal immune atlas of subarachnoid
hemorrhage from single-cell and spatial transcriptomics [scRNA-Seq]".
  * Mus musculus, CD45+ sorted immune cells from the injured cerebral hemisphere.
  * 5 samples, per-sample CellRanger "filtered" matrices (genes x cells):
      GSM9333532 = Sham (day0)
      GSM9333533 = SAH_24h_1
      GSM9333534 = SAH_24h_2
      GSM9333535 = SAH_72h_1
      GSM9333536 = SAH_72h_2

HARD-HONESTY RULE (non-negotiable): only real downloaded data is used. Nothing
is fabricated / imputed / invented.  Genes not present in the feature list are
reported as NOT MEASURED.  Comparisons that are under-powered (Sham has n=1
sample) are reported descriptively and explicitly caveated.  Cell-level p-values
are PSEUDO-REPLICATED (cells are not independent samples) and are labelled as
such; a sample-level descriptive summary is also given.  This is MOUSE CD45+
immune tissue -- NOT human, NOT SAH-CSF.

HDFI composite follows the paper's convention (see scripts/47_proteomics_hdfi.py):
    HDFI_mouse = z(Hmox1) - z(Fth1) - z(Ftl1) - z(Gpx4)
  (the "z(FTH1+FTL)" term is implemented as z(FTH1)+z(FTL), each gene z-scored
   ACROSS ALL CELLS with ddof=0, exactly as in the bulk/protein pipeline).
  Input values = log1p(CPM) expression per cell.

OUTPUTS (only these files are written; no other project file is touched):
  results/tables/mouse_sah_atlas_hdfi.csv
  results/mouse_sah_atlas_summary.json
  figures/Fig39_mouse_sah_atlas.png   +   .pdf
"""
import os, re, json, gzip, glob
import numpy as np
import pandas as pd
from scipy import io as sio, stats, sparse
import matplotlib; matplotlib.use("Agg")

# --- naturestyle: strip in-figure notes; split composites to independent panels; 5 cm footprint ---
import naturestyle as _ns
_ns.savefig_patch()
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

rng_seed = 0
RNG = np.random.default_rng(rng_seed)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUPP = os.path.join(ROOT, "data", "external", "GSE311852")
TAB  = os.path.join(ROOT, "results", "tables")
FIG  = os.path.join(ROOT, "figures")
for d in (TAB, FIG): os.makedirs(d, exist_ok=True)

# ---- sample registry ---------------------------------------------------------
SAMPLES = {
    "GSM9333532": dict(tag="Sham",      cond="Sham"),
    "GSM9333533": dict(tag="SAH_24h1",  cond="SAH24"),
    "GSM9333534": dict(tag="SAH_24h2",  cond="SAH24"),
    "GSM9333535": dict(tag="SAH_72h1",  cond="SAH72"),
    "GSM9333536": dict(tag="SAH_72h2",  cond="SAH72"),
}
CONDS = ["Sham", "SAH24", "SAH72"]
COND_LABEL = {"Sham": "Sham (d0)", "SAH24": "SAH 24h", "SAH72": "SAH 72h"}

# ---- HDFI gene program (mouse symbols) --------------------------------------
HDFI_GENES = ["Hmox1", "Fth1", "Ftl1", "Gpx4", "Acsl4", "Hp", "Cd163", "Tfrc", "Slc40a1"]

# ---- cell-type marker panels (mouse symbols) --------------------------------
MARKERS = {
    "Myeloid(microglia/macrophage)": ["P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Adgre1", "Cd68", "Sall1"],
    "Neutrophil":                    ["Elane", "Mpo", "CSF3R", "S100a8", "S100a9"],
    "T cell":                        ["Cd3e", "Cd3d", "Lck", "Il7r", "Cd247"],
    "B cell":                        ["Cd79a", "Cd79b", "Ms4a4a", "Pax5"],
    "NK":                            ["Ncr1", "Klrd1", "Gzma"],
    "Dendritic":                     ["Itgax", "Flt3", "Cd83", "Il12b"],
}
CT_ORDER = list(MARKERS.keys())

# ---- Okabe-Ito palette + rcParams (match scripts/47) -------------------------
OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
          purple="#CC79A7", sky="#56B4E9", grey="#555555", yellow="#F0E442")
COND_COLOR = {"Sham": OI["blue"], "SAH24": OI["orange"], "SAH72": OI["verm"]}
rcParams.update({"font.family": "Arial", "font.size": 9, "axes.linewidth": 0.9,
                 "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                 "axes.spines.right": False, "pdf.fonttype": 42})


def _path(gsm, kind):
    tag = SAMPLES[gsm]["tag"]
    ext = {"bar": "tsv", "fea": "tsv", "mtx": "mtx"}[kind]
    fn = f"{gsm}_filtered_{tag}_{kind if kind!='mtx' else 'matrix'}.{ext}.gz"
    if kind == "bar":
        fn = f"{gsm}_filtered_{tag}_barcodes.tsv.gz"
    elif kind == "fea":
        fn = f"{gsm}_filtered_{tag}_features.tsv.gz"
    else:
        fn = f"{gsm}_filtered_{tag}_matrix.mtx.gz"
    return os.path.join(SUPP, fn)


def load_features(gsm):
    """Return dict symbol -> list of row indices (handles duplicate symbols)."""
    sym2idx = {}
    with gzip.open(_path(gsm, "fea"), "rt") as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            name = parts[1]
            sym2idx.setdefault(name, []).append(i)
    return sym2idx


def load_barcodes(gsm):
    with gzip.open(_path(gsm, "bar"), "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def gzsize(path):
    return os.path.getsize(path) if os.path.exists(path) else 0


def aggregate_to_panel(M, sym2idx, panel):
    """Given sparse genes x cells matrix M and symbol map, return dense
    (n_cells x n_panel) raw counts, summing duplicate symbol rows.
    Missing gene -> all-zero column."""
    nc = M.shape[1]
    out = np.zeros((nc, len(panel)), dtype=np.float32)
    for j, g in enumerate(panel):
        idxs = sym2idx.get(g)
        if not idxs:
            continue
        block = M[idxs, :]
        col = np.asarray(block.sum(axis=0)).ravel()
        out[:, j] = col
    return out


def main():
    log = {}
    print("[1] loading 5 per-sample matrices ...")

    files_downloaded = []
    for gsm in SAMPLES:
        for kind in ("bar", "fea", "mtx"):
            p = _path(gsm, kind)
            files_downloaded.append({"file": os.path.basename(p), "bytes": gzsize(p)})

    PANEL = HDFI_GENES + sorted({g for v in MARKERS.values() for g in v})
    PANEL = list(dict.fromkeys(PANEL))          # preserve order, unique

    all_counts, all_meta = [], []
    panel_measured = {g: [] for g in PANEL}      # per-sample presence
    for gsm in SAMPLES:
        meta = SAMPLES[gsm]
        sym2idx = load_features(gsm)
        barcodes = load_barcodes(gsm)
        M = sio.mmread(_path(gsm, "mtx"))         # genes x cells
        M = sparse.csr_matrix(M)
        assert M.shape[1] == len(barcodes), f"{gsm}: matrix cols {M.shape[1]} != barcodes {len(barcodes)}"
        present = set(sym2idx.keys())
        for g in PANEL:
            panel_measured[g].append(g in present)
        sub = aggregate_to_panel(M, sym2idx, PANEL)          # n_cells x n_panel raw
        lib  = np.asarray(M.sum(axis=0)).ravel()             # total counts per cell (all genes)
        nnz  = np.asarray((M != 0).sum(axis=0)).ravel()      # detected genes per cell
        all_counts.append(sub)
        all_meta.append(dict(
            sample=gsm, tag=meta["tag"], cond=meta["cond"],
            barcodes=np.array(barcodes, dtype=object), lib=lib, detected_genes=nnz,
            n_cells=len(barcodes)))
        print(f"    {gsm} ({meta['cond']:>5}): {len(barcodes)} cells, "
              f"{M.shape[0]} genes, panel present={sum(1 for g in PANEL if g in present)}/{len(PANEL)}")

    counts = np.vstack(all_counts)                            # n_cells_all x n_panel
    lib    = np.concatenate([m["lib"] for m in all_meta])
    dg     = np.concatenate([m["detected_genes"] for m in all_meta])
    cond   = np.repeat([m["cond"] for m in all_meta], [m["n_cells"] for m in all_meta])
    samp   = np.repeat([m["sample"] for m in all_meta], [m["n_cells"] for m in all_meta])
    gsm_tag= np.repeat([m["tag"] for m in all_meta], [m["n_cells"] for m in all_meta])
    panel_idx = {g: j for j, g in enumerate(PANEL)}

    n_total = counts.shape[0]
    log["n_cells_raw"] = int(n_total)

    # ---- measured panel genes ----------------------------------------------
    measured = {g: bool(np.any(panel_measured[g])) for g in PANEL}
    hdfi_measured = [g for g in HDFI_GENES if measured[g]]
    hdfi_missing  = [g for g in HDFI_GENES if not measured[g]]
    hdfi_computable = all(measured[g] for g in ("Hmox1", "Fth1", "Ftl1", "Gpx4"))
    log["hdfi_measured"] = hdfi_measured
    log["hdfi_missing"] = hdfi_missing
    log["hdfi_computable"] = bool(hdfi_computable)
    print(f"    HDFI genes measured: {hdfi_measured}")
    print(f"    HDFI genes NOT measured: {hdfi_missing} ; composite computable={hdfi_computable}")

    # ---- QC: drop cells with <200 detected genes ---------------------------
    keep = dg >= 200
    log["n_cells_after_qc"] = int(keep.sum())
    print(f"[2] QC: kept {int(keep.sum())}/{n_total} cells (>=200 detected genes)")
    counts, lib, dg, cond, samp, gsm_tag = (counts[keep], lib[keep], dg[keep],
                                            cond[keep], samp[keep], gsm_tag[keep])

    # ---- subsample to <=15000 cells (seed 0) -------------------------------
    n_after = counts.shape[0]
    MAXCELLS = 15000
    if n_after > MAXCELLS:
        sel = RNG.choice(n_after, size=MAXCELLS, replace=False)
        sel.sort()
        counts, lib, dg, cond, samp, gsm_tag = (counts[sel], lib[sel], dg[sel],
                                               cond[sel], samp[sel], gsm_tag[sel])
        subsampled = True
    else:
        subsampled = False
    n_used = counts.shape[0]
    log["n_used"] = int(n_used)
    log["subsampled"] = bool(subsampled)
    log["subsample_max"] = MAXCELLS
    print(f"    analysis set n_used = {n_used} (subsampled={subsampled}, seed={rng_seed})")

    # ---- normalize to log1p(CPM) -------------------------------------------
    cpm = counts / lib[:, None] * 1e6
    logx = np.log1p(cpm)                                       # n_used x n_panel
    L = {g: logx[:, panel_idx[g]] for g in PANEL}              # gene -> per-cell log expr

    # ---- cell-type assignment (argmax of mean marker z-score) --------------
    def z(v):
        return (v - v.mean()) / (v.std(ddof=0) + 1e-12)
    scores = np.zeros((n_used, len(CT_ORDER)), dtype=float)
    panel_marker_used = {}
    for k, ct in enumerate(CT_ORDER):
        genes = [g for g in MARKERS[ct] if measured[g] and L[g].std() > 0]
        if len(genes) == 0:
            scores[:, k] = -np.inf
            panel_marker_used[ct] = []
            continue
        panel_marker_used[ct] = genes
        scores[:, k] = np.mean([z(L[g]) for g in genes], axis=0)
    celltype = np.array([CT_ORDER[i] for i in scores.argmax(axis=1)], dtype=object)
    log["markers_used_per_panel"] = panel_marker_used

    # ---- HDFI_mouse composite ---------------------------------------------
    hdfi = None
    if hdfi_computable:
        hdfi = z(L["Hmox1"]) - z(L["Fth1"]) - z(L["Ftl1"]) - z(L["Gpx4"])
        hdfi = hdfi.astype(float)
    else:
        log["hdfi_not_computable_reason"] = "one or more of Hmox1/Fth1/Ftl1/Gpx4 not measured"
    log["markers_missing"] = {ct: [g for g in MARKERS[ct] if not measured[g]] for ct in CT_ORDER}

    # ---- detection fraction per HDFI gene (cells with count>0) --------------
    detfrac = {g: float((counts[:, panel_idx[g]] > 0).mean()) for g in HDFI_GENES}
    log["detection_fraction"] = detfrac

    # =========================================================================
    # 5. STATISTICS
    # =========================================================================
    print("[5] statistics ...")
    stat_rows = []

    def bootstrap_median_diff(a, b, n=1000):
        """95% CI for median(a)-median(b) via non-parametric bootstrap."""
        rng = np.random.default_rng(rng_seed)
        diffs = []
        for _ in range(n):
            aa = rng.choice(a, size=len(a), replace=True)
            bb = rng.choice(b, size=len(b), replace=True)
            diffs.append(np.median(aa) - np.median(bb))
        return float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))

    results_stats = {}
    if hdfi_computable:
        # per-cell-type tests + myeloid focus
        for ct in CT_ORDER:
            mask = celltype == ct
            if mask.sum() == 0:
                continue
            group = {c: hdfi[mask & (cond == c)] for c in CONDS}
            nz = {c: group[c].size for c in CONDS}
            present = [c for c in CONDS if nz[c] >= 2]
            entry = dict(n_cells_by_cond={c: int(nz[c]) for c in CONDS},
                         medians={c: (float(np.median(group[c])) if nz[c] > 0 else None) for c in CONDS})
            # Kruskal-Wallis (cell-level, PSEUDO-REPLICATED)
            if len(present) >= 2:
                kw = stats.kruskal(*[group[c] for c in present])
                entry["kruskal_H"] = float(kw.statistic)
                entry["kruskal_p_celllevel"] = float(kw.pvalue)
                entry["kw_conds_tested"] = present
            # pairwise Mann-Whitney + effect size
            pw = {}
            for i, a in enumerate(CONDS):
                for b in CONDS[i+1:]:
                    if nz[a] >= 2 and nz[b] >= 2:
                        u, p = stats.mannwhitneyu(group[a], group[b], alternative="two-sided")
                        lo, hi = bootstrap_median_diff(group[b], group[a])  # b - a
                        pw[f"{a}_vs_{b}"] = dict(
                            p_celllevel=float(p),
                            median_diff_b_minus_a=float(np.median(group[b]) - np.median(group[a])),
                            boot95CI_median_diff=[lo, hi],
                            n_a=int(nz[a]), n_b=int(nz[b]))
            entry["pairwise"] = pw
            results_stats[ct] = entry
            stat_rows.append(dict(section="HDFI_stats", celltype=ct,
                                  metric="KW_p_celllevel", value=entry.get("kruskal_p_celllevel"),
                                  n_cells=int(mask.sum())))

        # sample-level descriptive (Sham n=1 -> underpowered)
        sample_level = {}
        for ct in ["Myeloid(microglia/macrophage)"]:
            mask = (celltype == ct)
            per_sample = {}
            for s in np.unique(samp[mask]):
                v = hdfi[mask & (samp == s)]
                per_sample[str(s)] = dict(tag=SAMPLES[str(s)]["tag"], cond=SAMPLES[str(s)]["cond"],
                                          n_cells=int(v.size),
                                          mean_hdfi=float(v.mean()), median_hdfi=float(np.median(v)))
            sample_level[ct] = per_sample
        results_stats["_sample_level_myeloid"] = sample_level
        log["sample_level_myeloid"] = sample_level
        log["pseudo_replication_caveat"] = (
            "Cell-level Kruskal-Wallis / Mann-Whitney p-values treat each cell as an "
            "independent observation. Cells within a sample are NOT independent "
            "(pseudo-replication); these p-values are descriptive of cell-population "
            "separation, NOT inferential evidence at the biological-sample level. "
            "Only 5 samples exist (Sham n=1, SAH24 n=2, SAH72 n=2), so sample-level "
            "hypothesis testing is under-powered and is reported descriptively only.")

        # correlation across cells (Spearman)
        corr = {}
        proxy = (z(L["Acsl4"]) + z(L["Tfrc"]) - z(L["Gpx4"])) if (measured["Acsl4"] and measured["Tfrc"] and measured["Gpx4"]) else (L["Acsl4"] + L["Tfrc"] - L["Gpx4"])
        r_h, p_h = stats.spearmanr(hdfi, L["Hmox1"])
        r_p, p_p = stats.spearmanr(hdfi, proxy)
        corr["HDFI_vs_Hmox1_logexpr"] = dict(rho=float(r_h), p=float(p_h), n=int(n_used))
        corr["HDFI_vs_ferroptosis_proxy"] = dict(rho=float(r_p), p=float(p_p), n=int(n_used),
                                                 proxy_def="(Acsl4+Tfrc-Gpx4)")
        results_stats["_correlations"] = corr
        log["correlations"] = corr
        print(f"    Spearman HDFI~Hmox1 rho={r_h:.3f} p={p_h:.2e}")
        print(f"    Spearman HDFI~ferroptosis_proxy rho={r_p:.3f} p={p_p:.2e}")
        my = results_stats.get("Myeloid(microglia/macrophage)", {})
        log["myeloid_summary"] = dict(medians=my.get("medians"),
                                      kw_p=my.get("kruskal_p_celllevel"),
                                      pairwise=my.get("pairwise"))

    # =========================================================================
    # gene mean log-expression per celltype x condition -> CSV
    # =========================================================================
    print("[6] writing CSV ...")
    rows = []
    for ct in CT_ORDER:
        for c in CONDS:
            mask = (celltype == ct) & (cond == c)
            nc = int(mask.sum())
            if nc == 0:
                continue
            r = dict(section="gene_mean_logexpr", celltype=ct, condition=c, n_cells=nc)
            for g in HDFI_GENES:
                r[f"{g}_mean_logCPM"] = round(float(L[g][mask].mean()), 4)
                r[f"{g}_pct_detected"] = round(float((counts[mask, panel_idx[g]] > 0).mean()) * 100, 2)
            if hdfi_computable:
                r["HDFI_mouse_median"] = round(float(np.median(hdfi[mask])), 4)
                r["HDFI_mouse_mean"] = round(float(hdfi[mask].mean()), 4)
            rows.append(r)
    # HDFI test-stat rows (cell-level, myeloid focus + per-type KW)
    if hdfi_computable:
        for ct, e in results_stats.items():
            if ct.startswith("_"):
                continue
            rows.append(dict(section="HDFI_KW_test", celltype=ct, condition="all",
                             n_cells=sum(e["n_cells_by_cond"].values()),
                             kruskal_p_celllevel=e.get("kruskal_p_celllevel"),
                             med_Sham=e["medians"]["Sham"], med_SAH24=e["medians"]["SAH24"],
                             med_SAH72=e["medians"]["SAH72"],
                             n_Sham=e["n_cells_by_cond"]["Sham"],
                             n_SAH24=e["n_cells_by_cond"]["SAH24"],
                             n_SAH72=e["n_cells_by_cond"]["SAH72"],
                             notes="CELL-LEVEL p = pseudo-replicated; descriptive only"))
            for pair, d in e.get("pairwise", {}).items():
                rows.append(dict(section="HDFI_pairwise", celltype=ct, condition=pair,
                                 n_cells=d["n_a"] + d["n_b"],
                                 mannwhitney_p_celllevel=round(d["p_celllevel"], 6),
                                 median_diff=round(d["median_diff_b_minus_a"], 4),
                                 boot95_lo=round(d["boot95CI_median_diff"][0], 4),
                                 boot95_hi=round(d["boot95CI_median_diff"][1], 4),
                                 notes="CELL-LEVEL p = pseudo-replicated; CI = bootstrap of median diff"))
    # detection rows
    for g in HDFI_GENES:
        rows.append(dict(section="gene_detection", gene=g,
                         measured_in_feature_space=bool(measured[g]),
                         mean_pct_detected_across_used_cells=round(detfrac[g]*100, 3)))
    csv = pd.DataFrame(rows)
    out_csv = os.path.join(TAB, "mouse_sah_atlas_hdfi.csv")
    csv.to_csv(out_csv, index=False)
    print("    wrote", out_csv, f"({os.path.getsize(out_csv)} bytes)")

    # =========================================================================
    # JSON summary
    # =========================================================================
    comp_frac = {}
    for c in CONDS:
        sub = cond == c
        tot = sub.sum()
        comp_frac[c] = {ct: float(((celltype == ct) & sub).sum() / tot) if tot else 0.0 for ct in CT_ORDER}
    log["composition_fraction_by_cond"] = comp_frac
    log["n_cells_by_cond"] = {c: int((cond == c).sum()) for c in CONDS}
    log["files_downloaded"] = files_downloaded
    log["dataset"] = "GSE311852 (Mus musculus, CD45+ immune, injured hemisphere; in-vivo cross-species corroboration -- NOT human, NOT SAH-CSF)"
    log["hdfi_formula"] = "z(Hmox1) - z(Fth1) - z(Ftl1) - z(Gpx4), z across all cells ddof=0, input log1p(CPM) [matches scripts/47]"
    log["results_stats"] = results_stats

    out_json = os.path.join(FIG, "..", "results", "mouse_sah_atlas_summary.json")
    out_json = os.path.normpath(os.path.join(ROOT, "results", "mouse_sah_atlas_summary.json"))
    with open(out_json, "w") as fh:
        json.dump(log, fh, indent=2, default=str)
    print("    wrote", out_json, f"({os.path.getsize(out_json)} bytes)")

    # =========================================================================
    # FIGURE
    # =========================================================================
    print("[7] drawing figure ...")
    fig = plt.figure(figsize=(11.0, 8.2))
    gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.34,
                  left=0.10, right=0.97, top=0.92, bottom=0.10)

    # ---- Panel 1: composition ----
    ax1 = fig.add_subplot(gs[0, 0])
    w = 0.8
    bottoms = np.zeros(len(CONDS))
    ct_colors = [OI["blue"], OI["orange"], OI["green"], OI["purple"], OI["sky"], OI["yellow"]]
    for k, ct in enumerate(CT_ORDER):
        vals = np.array([comp_frac[c][ct] * 100 for c in CONDS])
        ax1.bar(np.arange(len(CONDS)), vals, w, bottom=bottoms, color=ct_colors[k], label=ct, edgecolor="white", linewidth=0.6)
        bottoms += vals
    ax1.set_xticks(np.arange(len(CONDS)))
    ax1.set_xticklabels([COND_LABEL[c] for c in CONDS])
    ax1.set_ylabel("Cell composition (%)")
    ax1.set_title("Cell-type composition of mouse CD45+ atlas", fontweight="bold")
    ax1.set_ylim(0, 100)
    ax1.legend(fontsize=6.0, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, frameon=False)

    # ---- Panel 2: HDFI by cell type x condition (box) ----
    ax2 = fig.add_subplot(gs[0, 1])
    if hdfi_computable:
        types_with_data = [ct for ct in CT_ORDER if ((celltype == ct)).sum() > 0]
        pos = np.arange(len(types_with_data))
        bw = 0.26
        data = []
        for j, c in enumerate(CONDS):
            dd = [hdfi[(celltype == ct) & (cond == c)] for ct in types_with_data]
            bp = ax2.boxplot(dd, positions=pos + (j-1)*bw, widths=bw*0.9, patch_artist=True,
                             boxprops=dict(linewidth=0.7), whiskerprops=dict(linewidth=0.7),
                             capprops=dict(linewidth=0.7), medianprops=dict(color="black", lw=1.1),
                             flierprops=dict(marker=".", markersize=1.2, alpha=0.3))
            for b in bp["boxes"]:
                b.set_facecolor(COND_COLOR[c]); b.set_alpha(0.85)
            data.append(dd)
        ax2.set_xticks(pos)
        ax2.set_xticklabels([ct.replace("Myeloid(microglia/macrophage)", "Myeloid") for ct in types_with_data],
                            rotation=30, ha="right", fontsize=7)
        ax2.set_ylabel("HDFI_mouse score")
        ax2.axhline(0, color="grey", lw=0.6, ls="--")
        ax2.set_title("HDFI_mouse by cell type x condition", fontweight="bold")
        ax2.legend(handles=[Patch(facecolor=COND_COLOR[c], label=COND_LABEL[c]) for c in CONDS],
                   fontsize=7, frameon=False, ncol=3,
                   loc="upper center", bbox_to_anchor=(0.5, -0.28))
    else:
        ax2.text(0.5, 0.5, "HDFI composite\nNOT COMPUTABLE", ha="center", va="center", transform=ax2.transAxes, color=OI["verm"])
        ax2.set_title("HDFI_mouse by cell type x condition", fontweight="bold")

    # ---- Panel 3: HDFI gene mean log expr heatmap (celltype x cond rows), spans bottom row ----
    ax3 = fig.add_subplot(gs[1, :])
    genes_for_hm = hdfi_measured if hdfi_measured else HDFI_GENES
    row_labels, row_keys = [], []
    for ct in CT_ORDER:
        for c in CONDS:
            mask = (celltype == ct) & (cond == c)
            if mask.sum() >= 5:
                row_labels.append(f"{ct.replace('Myeloid(microglia/macrophage)','Myeloid')} | {COND_LABEL[c]}")
                row_keys.append(mask)
    mat = np.zeros((len(row_keys), len(genes_for_hm)))
    for i, mask in enumerate(row_keys):
        for j, g in enumerate(genes_for_hm):
            mat[i, j] = L[g][mask].mean()
    # z across rows for visual contrast
    show = mat.copy()
    if show.shape[0] > 1:
        show = (show - show.mean(axis=0)) / (show.std(axis=0) + 1e-9)
    im = ax3.imshow(show, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax3.set_yticks(np.arange(len(row_labels)))
    ax3.set_yticklabels(row_labels, fontsize=6)
    ax3.set_xticks(np.arange(len(genes_for_hm)))
    ax3.set_xticklabels(genes_for_hm, rotation=45, ha="right", fontsize=7)
    ax3.set_title("HDFI gene expression (mean logCPM, row-z)", fontweight="bold")
    ax3.set_ylabel("Cell type | condition")
    fig.colorbar(im, ax=ax3, fraction=0.022, pad=0.02, label="row-z")

    fig.suptitle("Mouse SAH immune atlas (GSE311852): HDFI programme across cell types over time",
                 fontweight="bold", fontsize=11)
    png = os.path.join(FIG, "Fig39_mouse_sah_atlas.png")
    pdf = os.path.join(FIG, "Fig39_mouse_sah_atlas.pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    _written = [p for p in (png, pdf) if os.path.exists(p)] or sorted(glob.glob(png[:-4] + "_*.png"))
    print("    wrote", ", ".join(os.path.basename(p) for p in _written) or "(panels split)")

    # =========================================================================
    # printed SUMMARY
    # =========================================================================
    print("\n=== SUMMARY ===")
    print("Dataset:", log["dataset"])
    print("HDFI genes detected in mouse:", hdfi_measured)
    print("HDFI genes NOT detected:", hdfi_missing if hdfi_missing else "none")
    print("n_used cells:", n_used, "| by cond:", log["n_cells_by_cond"])
    myl = results_stats.get("Myeloid(microglia/macrophage)", {}) if hdfi_computable else {}
    med = myl.get("medians", {})
    def gfmt(v):
        if v is None: return "NA"
        return ("<2e-16" if v < 2e-16 else f"{v:.2g}")
    if hdfi_computable and myl:
        print("Myeloid HDFI_mouse medians:",
              {COND_LABEL[c]: (round(med[c],3) if med.get(c) is not None else None) for c in CONDS})
        print("Myeloid KW cell-level p:", gfmt(myl.get("kruskal_p_celllevel")))
        pw = myl.get("pairwise", {})
        for pair in ["Sham_vs_SAH24", "Sham_vs_SAH72", "SAH24_vs_SAH72"]:
            d = pw.get(pair)
            if d:
                print(f"  {pair}: p={gfmt(d['p_celllevel'])} medDiff(b-a)={d['median_diff_b_minus_a']:+.3f} "
                      f"CI95=[{d['boot95CI_median_diff'][0]:+.3f},{d['boot95CI_median_diff'][1]:+.3f}] "
                      f"(n_a={d['n_a']},n_b={d['n_b']})")
        print("Spearman HDFI~Hmox1:", log["correlations"]["HDFI_vs_Hmox1_logexpr"])
        print("Spearman HDFI~ferroptosis_proxy:", log["correlations"]["HDFI_vs_ferroptosis_proxy"])
        print("sample-level myeloid (underpowered):")
        for ct, ps in log["sample_level_myeloid"].items():
            for s, d in ps.items():
                print(f"   {d['tag']} ({d['cond']}): n={d['n_cells']} mean={d['mean_hdfi']:.3f} median={d['median_hdfi']:.3f}")
    else:
        print("HDFI composite NOT COMPUTABLE ->", log.get("hdfi_not_computable_reason"))
    print("NOTE:", log["pseudo_replication_caveat"])
    return log


if __name__ == "__main__":
    main()
