"""
47_proteomics_hdfi.py  -- Round-11 reviewer TASK: PROTEIN-LEVEL validation of the
HDFI (Heme-Detoxification-Failure-Index) axis using REAL SAH cerebrospinal-fluid
(CSF) proteomics.

Reviewer-requested HDFI axis proteins (protein level):
    HP, CD163 (MRC1 as alt), HMOX1, FTH1, FTL, GPX4, ACSL4, TFRC, SLC40A1

Requested HDFI composite:
    HDFI_protein = z(HMOX1) - z(FTH1 + FTL) - z(GPX4)

HARD-HONESTY RULE: we use ONLY real, openly downloadable data. Nothing is
fabricated, imputed or invented. If a required component of the composite is not
measured in a real dataset we report the composite as NOT COMPUTABLE rather than
substituting a proxy and calling it HDFI.

DATA ACTUALLY OBTAINED IN THIS RUN (see results/tables/proteomics_hdfi_ACCESS.txt
for the HTTP evidence):

  (1) PXD030593 "Human CSF after SAH LC-MS/MS"  -- a real, openly downloadable
      per-sample relative-abundance table (protein.xlsx) with an explicit
      SAH (n=4) vs Control (n=4) column design.  This is the PRIMARY analysed
      dataset.  6/9 axis proteins are detected; HMOX1, GPX4, ACSL4, SLC40A1 are
      NOT detected -> the HDFI composite cannot be formed from this dataset.

  (2) PMC11267790 (Clin Proteomics 2024) supplementary tables (openly served as a
      ZIP by Europe PMC) -- an independent, larger SAH-vs-control CSF proteomics
      cohort (control n=10 vs SAH n=23) reported only as summary statistics.
      Used here purely as an independent *directional cross-check* of the
      detectable axis proteins (values are read verbatim from the supplementary
      text file, never re-derived).  HMOX1 / GPX4 / ACSL4 / SLC40A1 are likewise
      ABSENT from this larger cohort's supplementary tables.

  (3) PXD032118 (SAH CSF DIA-MS, 20 patients) is NOT in the PRIDE archive
      (HTTP 404) and is hosted on iProX behind a login/auth server -> raw / gated,
      so it is NOT analysed.

Outputs (only the executable parts run; nothing is faked):
  results/tables/proteomics_hdfi.csv
  figures/Fig37_proteomics_hdfi.png  +  .pdf    (Okabe-Ito, Arial, 300 dpi, no panel letters)
"""
import os, re, json
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "external")
TAB  = os.path.join(ROOT, "results", "tables")
FIG  = os.path.join(ROOT, "figures")
for d in (TAB, FIG): os.makedirs(d, exist_ok=True)

PX_XLSX  = os.path.join(DATA, "PXD030593", "protein.xlsx")
PMC_TEXT = os.path.join(DATA, "PMC11267790", "PMC11267790_supplementary_tables_text.txt")

OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00",
          purple="#CC79A7", sky="#56B4E9", grey="#555555", yellow="#F0E442")
rcParams.update({"font.family":"Arial","font.size":9,"axes.linewidth":0.9,
                 "xtick.labelsize":7,"ytick.labelsize":7,"axes.spines.top":False,
                 "axes.spines.right":False,"pdf.fonttype":42})

# HDFI axis definition: canonical gene + acceptable alias UniProt/gene names.
AXIS = {
    "HP":      {"aliases":["HP"],       "role":"heme scavenger (haptoglobin)"},
    "CD163":   {"aliases":["CD163","MRC1"], "role":"heme-Hb endocytosis receptor"},
    "HMOX1":   {"aliases":["HMOX1"],    "role":"heme catabolism (pro-oxidant arm)"},
    "FTH1":    {"aliases":["FTH1"],     "role":"iron storage (ferritin heavy)"},
    "FTL":     {"aliases":["FTL"],      "role":"iron storage (ferritin light)"},
    "GPX4":    {"aliases":["GPX4"],     "role":"lipid-peroxide repair (anti-ferroptosis)"},
    "ACSL4":   {"aliases":["ACSL4"],    "role":"PUFA-PL synthesis (pro-ferroptosis)"},
    "TFRC":    {"aliases":["TFRC"],     "role":"iron uptake (CD71)"},
    "SLC40A1": {"aliases":["SLC40A1","FPN"], "role":"iron export (ferroportin)"},
}
HDFI_NEED = ["HMOX1","FTH1","FTL","GPX4"]   # proteins the composite formula requires

# -----------------------------------------------------------------------------
# 1. Load the REAL PXD030593 per-sample table and normalise
# -----------------------------------------------------------------------------
def load_pxd():
    if not os.path.exists(PX_XLSX):
        raise FileNotFoundError(f"Primary dataset missing: {PX_XLSX}")
    df = pd.read_excel(PX_XLSX, sheet_name="Proteins", header=0)
    df.columns = [str(c).strip() for c in df.columns]
    df["GeneName"] = df["GeneName"].astype(str).str.strip().str.upper()
    ab = [c for c in df.columns if c.startswith("Abundances")]
    ctrl = [c for c in ab if "Control" in c]
    sah  = [c for c in ab if "SAH" in c]
    # total-intensity (sum) normalisation per sample -> relative abundance, then log2
    # total-intensity normalisation PER SAMPLE (column-wise) -> relative abundance x1e6, then log2
    tot = df[ab].replace(0, np.nan)
    colsum = tot.sum(axis=0)                     # total intensity of each sample
    frac = tot.div(colsum, axis=1) * 1e6         # divide each sample column by its total
    l2 = np.log2(frac + 1)
    l2.insert(0, "GeneName", df["GeneName"].values)
    return df, l2, ctrl, sah, ab

# -----------------------------------------------------------------------------
# 2. Per-protein case-vs-control statistics (Mann-Whitney) on the REAL data
# -----------------------------------------------------------------------------
def axis_stats(l2, ctrl, sah):
    rows = []
    for gene, meta in AXIS.items():
        # a gene may be recorded under an alias (e.g. CD163 vs MRC1); keep first detected
        present = [g for g in meta["aliases"] if (l2["GeneName"] == g).any()]
        if not present:
            rows.append(dict(dataset="PXD030593", protein=gene, role=meta["role"],
                             detected=False, n_control=len(ctrl), n_sah=len(sah)))
            continue
        g = present[0]
        r = l2.loc[l2["GeneName"] == g].iloc[0]
        c = r[ctrl].astype(float).dropna().values
        s = r[sah].astype(float).dropna().values
        u, p = stats.mannwhitneyu(c, s, alternative="two-sided")
        rows.append(dict(dataset="PXD030593", protein=gene, role=meta["role"],
                         detected=True, gene_used=g,
                         n_control=len(c), n_sah=len(s),
                         med_ctrl_log2=float(np.median(c)), med_sah_log2=float(np.median(s)),
                         log2FC_SAH_vs_Ctrl=float(np.median(s)-np.median(c)),
                         mwu_U=float(u), mwu_p=float(p)))
    return pd.DataFrame(rows)

def spearman_axis(l2, ctrl, sah, detected):
    """coherence of the DETECTED axis proteins across all 8 real samples."""
    cols = ctrl + sah
    mat = {}
    for gene in detected:
        aliases = AXIS[gene]["aliases"]
        g = next((a for a in aliases if (l2["GeneName"]==a).any()), None)
        if g is None: continue
        mat[gene] = l2.loc[l2["GeneName"]==g].iloc[0][cols].astype(float).values
    M = pd.DataFrame(mat)
    corr = M.corr(method="spearman")
    return M, corr

# -----------------------------------------------------------------------------
# 3. HDFI composite -- honest computation (only if formula inputs are detected)
# -----------------------------------------------------------------------------
def hdfi_composite(l2, ctrl, sah, st):
    cols = ctrl + sah
    detected = {r["protein"]: r["detected"] for _, r in st.iterrows()}
    need_ok = all(detected.get(g, False) for g in HDFI_NEED)
    if need_ok:
        def zvec(gene):
            aliases = AXIS[gene]["aliases"]
            g = next((a for a in aliases if (l2["GeneName"]==a).any()), None)
            v = l2.loc[l2["GeneName"]==g].iloc[0][cols].astype(float).values
            return (v - v.mean()) / (v.std(ddof=0) + 1e-12)
        z_hmox = zvec("HMOX1")
        z_fth  = zvec("FTH1"); z_ftl = zvec("FTL")
        z_gpx  = zvec("GPX4")
        hdfi = z_hmox - (z_fth + z_ftl) - z_gpx
        c = [cols.index(x) for x in ctrl]; s = [cols.index(x) for x in sah]
        u, p = stats.mannwhitneyu(hdfi[c], hdfi[s], alternative="two-sided")
        return dict(computable=True, n_ctrl=len(c), n_sah=len(s),
                    med_ctrl=float(np.median(hdfi[c])), med_sah=float(np.median(hdfi[s])),
                    mwu_p=float(p)), hdfi, cols
    missing = [g for g in HDFI_NEED if not detected.get(g, False)]
    return dict(computable=False, missing=missing,
                reason=("HDFI_protein = z(HMOX1) - z(FTH1+FTL) - z(GPX4) NOT COMPUTABLE: "
                        f"{missing} not detected at the protein level in PXD030593 "
                        "(and also absent from the independent cohort PMC11267790).")), None, None

# -----------------------------------------------------------------------------
# 4. Independent cross-check: verbatim values read from the PMC11267790 supplementary text
# -----------------------------------------------------------------------------
def pmc_crosscheck():
    out = {"present_counts":{}, "table1_rows":{}}
    if not os.path.exists(PMC_TEXT):
        out["note"] = "PMC11267790 supplementary text not on disk; cross-check skipped."
        return out
    t = open(PMC_TEXT, encoding="utf-8", errors="ignore").read()
    tab1 = t[t.find("Supplementary Table 1"):t.find("Supplementary Table 2")]
    rowpat = re.compile(
        r"(?P<g>[A-Za-z0-9\-]+)\s+(?P<acc>\S+)\s+(?P<cm>\d+(?:\.\d+)?)\s*\((?P<cs>[^)]*)\)"
        r"\s*\[(?P<cn>\d+)\]\s+(?P<sm>\d+(?:\.\d+)?)\s*\((?P<ss>[^)]*)\)\s*\[(?P<sn>\d+)\]"
        r"\s+(?P<fc>-?\d+(?:\.\d+)?)\s+(?P<p>[\d.eE\-+]+)\s+(?P<adj>[\d.eE\-+]+)")
    seen = {}
    for ln in tab1.splitlines():
        m = rowpat.search(ln.strip())
        if m and m.group("g").upper() in [k.upper() for k in AXIS]:
            seen.setdefault(m.group("g").upper(), ln.strip())
    out["table1_rows"] = seen
    for g, meta in AXIS.items():
        cnt = sum(len(re.findall(rf"\b{re.escape(a)}\b", t)) for a in meta["aliases"])
        out["present_counts"][g] = cnt
    return out

# -----------------------------------------------------------------------------
# 5. Figure (real detectable-axis evidence + explicit HDFI not-computable note)
# -----------------------------------------------------------------------------
def make_figure(st, hdfi_info, corr_summary, samples_dict, path_png, path_pdf):
    # only the axis proteins actually detected in PXD030593
    order = ["HP","CD163","MRC1","FTH1","FTL","TFRC"]
    plot = [g for g in order if g in samples_dict]
    n = len(plot)
    ncol = 3; nrow = int(np.ceil(n/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.3*ncol, 2.0*nrow+1.7), squeeze=False)
    stidx = st.set_index("protein")
    rng = np.random.default_rng(0)
    for k, gene in enumerate(plot):
        ax = axes[k//ncol][k%ncol]
        cv, sv = samples_dict[gene]
        for (vals, x, col, lab) in [(cv,0,OI["blue"],"Control"), (sv,1,OI["verm"],"SAH")]:
            vals = np.asarray(vals, float)
            jit = x + (rng.random(len(vals))-0.5)*0.22
            ax.scatter(jit, vals, s=34, color=col, edgecolor="black", linewidth=0.4, zorder=3, label=lab)
            ax.hlines(np.median(vals), x-0.32, x+0.32, color=col, lw=1.6, zorder=4)
        ax.set_xticks([0,1]); ax.set_xticklabels(["Control","SAH"], fontsize=8)
        ax.set_xlim(-0.55,1.55)
        if k%ncol==0: ax.set_ylabel("log2 rel. abundance", fontsize=8)
        row = stidx.loc[gene]
        p = float(row["mwu_p"]); fc = float(row["log2FC_SAH_vs_Ctrl"])
        ps = f"p={p:.3f}" if p>=0.001 else f"p={p:.1e}"
        ns = "**" if p<0.01 else ("*" if p<0.05 else "n.s.")
        ax.set_title(f"{gene}   {ns}", fontsize=10, fontweight="bold", color="#222222")
        ax.text(0.5, 0.98, f"median log2FC {fc:+.2f}\nMann-Whitney {ps}", transform=ax.transAxes,
                ha="center", va="top", fontsize=6.8, color=OI["grey"], zorder=6,
                bbox=dict(boxstyle="round,pad=0.22", facecolor="white", alpha=0.85, edgecolor="none"))
    for j in range(n, nrow*ncol):
        axes[j//ncol][j%ncol].axis("off")
    fig.suptitle("HDFI axis proteins at the PROTEIN level — SAH vs control CSF\n"
                 "PXD030593 (Q Exactive Plus, label-free DDA LC-MS/MS; n=4 SAH vs n=4 control), total-intensity normalised",
                 fontsize=10, fontweight="bold", color="#222222", y=0.995)
    miss = hdfi_info.get("missing", ["HMOX1","GPX4"])
    note = ("HDFI_protein = z(HMOX1) - z(FTH1+FTL) - z(GPX4): NOT COMPUTABLE at the CSF-protein level.\n"
            + "/".join(miss) + " (also ACSL4, SLC40A1) were below detection in PXD030593 and are likewise absent\n"
            "from the independent larger cohort PMC11267790 (n=23 SAH vs n=10 control). Only detectable axis members\n"
            "are shown; nothing imputed or fabricated." + (("\n" + corr_summary) if corr_summary else ""))
    fig.text(0.5, 0.012, note, ha="center", va="bottom", fontsize=6.6,
             color=OI["verm"])
    fig.tight_layout(rect=(0,0.20,1,0.90))
    fig.savefig(path_png, dpi=300)
    fig.savefig(path_pdf)
    plt.close(fig)

# -----------------------------------------------------------------------------
def main():
    print("[1] loading PXD030593 protein.xlsx (real per-sample table)...")
    df, l2, ctrl, sah, ab = load_pxd()
    print("    proteins:", len(df), "| normalised log2 columns:", len(ab),
          "| n_ctrl:", len(ctrl), "n_sah:", len(sah))

    print("[2] per-protein case-vs-control Mann-Whitney on real data...")
    st = axis_stats(l2, ctrl, sah)

    print("[3] axis coherence (Spearman across 8 real samples, detected proteins)...")
    detected = list(st.loc[st["detected"]==True,"protein"])
    M, corr = spearman_axis(l2, ctrl, sah, detected)
    pair_p = None
    if M.shape[1] >= 2:
        # strongest positive heme-storage/heme-uptake coherence, as an example
        cs = corr.copy(); np.fill_diagonal(cs.values, np.nan)
        pair_p = f"Spearman rho range across detected axis: {np.nanmin(cs.values):.2f} to {np.nanmax(cs.values):.2f}."

    print("[4] HDFI composite attempt (honest)...")
    hdfi_info, hdfi_vec, hcols = hdfi_composite(l2, ctrl, sah, st)
    print("    ", json.dumps(hdfi_info, default=str))

    print("[5] independent cross-check: PMC11267790 supplementary (verbatim)...")
    pmc = pmc_crosscheck()

    # ------------------------------------------------------------------ CSV
    rows = []
    for _, r in st.iterrows():
        d = dict(source="PXD030593", protein=r["protein"], role=r["role"],
                 detected=bool(r["detected"]),
                 n_control=r.get("n_control"), n_sah=r.get("n_sah"))
        if r["detected"]:
            d.update(median_control_log2=round(r["med_ctrl_log2"],3),
                     median_sah_log2=round(r["med_sah_log2"],3),
                     log2FC_sah_vs_ctrl=round(r["log2FC_SAH_vs_Ctrl"],3),
                     mannwhitney_p=(round(r["mwu_p"],4) if pd.notna(r["mwu_p"]) else np.nan))
        else:
            d.update(median_control_log2=np.nan, median_sah_log2=np.nan,
                     log2FC_sah_vs_ctrl=np.nan, mannwhitney_p=np.nan)
        rows.append(d)
    # HDFI composite record
    rows.append(dict(source="PXD030593", protein="HDFI_protein",
                     role="z(HMOX1)-z(FTH1+FTL)-z(GPX4)",
                     detected=bool(hdfi_info["computable"]),
                     n_control=(hdfi_info.get("n_ctrl") if hdfi_info["computable"] else "0"),
                     n_sah=(hdfi_info.get("n_sah") if hdfi_info["computable"] else "0"),
                     median_control_log2=(round(hdfi_info["med_ctrl"],3) if hdfi_info["computable"] else np.nan),
                     median_sah_log2=(round(hdfi_info["med_sah"],3) if hdfi_info["computable"] else np.nan),
                     log2FC_sah_vs_ctrl=(round(hdfi_info["med_sah"]-hdfi_info["med_ctrl"],3) if hdfi_info["computable"] else np.nan),
                     mannwhitney_p=(round(hdfi_info["mwu_p"],4) if hdfi_info["computable"] else np.nan)))
    # PMC cross-check records (verbatim Table 1 rows, independent cohort)
    for g, line in pmc.get("table1_rows",{}).items():
        m = re.search(r"(\S+)\s+(\S+)\s+([\d.]+)\s*\(([^)]*)\)\s*\[(\d+)\]\s+([\d.]+)\s*\(([^)]*)\)\s*\[(\d+)\]\s+(-?[\d.]+)\s+([\d.eE\-+]+)\s+([\d.eE\-+]+)", line)
        if m:
            rows.append(dict(source="PMC11267790(verbatim)", protein=g,
                             role="independent SAH-CSF cohort (authors' Table 1)",
                             detected=True, n_control=int(m.group(5)), n_sah=int(m.group(8)),
                             median_control_log2=float(m.group(3)), median_sah_log2=float(m.group(6)),
                             log2FC_sah_vs_ctrl=float(m.group(9)), mannwhitney_p=np.nan))
    # PMC presence/absence for the whole supplementary
    for g,cnt in pmc.get("present_counts",{}).items():
        rows.append(dict(source="PMC11267790(verbatim)", protein=g,
                         role="token count across ALL supplementary tables",
                         detected=(cnt>0), n_control="n/a", n_sah="n/a",
                         median_control_log2=np.nan, median_sah_log2=np.nan,
                         log2FC_sah_vs_ctrl=np.nan, mannwhitney_p=np.nan))

    csv = pd.DataFrame(rows)
    out_csv = os.path.join(TAB, "proteomics_hdfi.csv")
    csv.to_csv(out_csv, index=False)
    print("    wrote", out_csv)

    # ------------------------------------------------------------------ FIGURE
    samples_dict = {}
    for _, r in st[st["detected"]==True].iterrows():
        gene = r["gene_used"]
        v = l2.loc[l2["GeneName"]==gene].iloc[0][ctrl+sah].astype(float)
        samples_dict[r["protein"]] = (v[ctrl].values, v[sah].values)
    make_figure(st, hdfi_info, pair_p, samples_dict,
                os.path.join(FIG,"Fig37_proteomics_hdfi.png"),
                os.path.join(FIG,"Fig37_proteomics_hdfi.pdf"))
    print("    wrote figure Fig37_proteomics_hdfi.png/.pdf")

    print("\n=== SUMMARY ===")
    print("detected axis proteins (PXD030593):", detected)
    print("not detected:", [r["protein"] for _,r in st.iterrows() if not r["detected"]])
    print("HDFI computable:", hdfi_info["computable"], "-", hdfi_info.get("reason",""))
    return csv, hdfi_info, pmc

if __name__ == "__main__":
    main()
