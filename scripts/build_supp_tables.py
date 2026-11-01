"""build_supp_tables.py -- emit submission Supplementary Tables from results/*.csv|.json.

Numbers are read directly from the analysis outputs (no hand entry) so every value is
traceable. Writes SUPP_tables_body.md (captions bold; three-line-table content only,
no notes). Tables S1-S14, contiguous, ordered by first appearance in the Results text.
"""
import os, json, math, re
import pandas as pd
from scipy.stats import chi2, spearmanr

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
TAB = os.path.join(RES, "tables")
OUT = os.path.join(BASE, "manuscript", "SUPP_tables_body.md")

def pval(x):
    try:
        x = float(x)
    except Exception:
        return str(x)
    if x != x:
        return "NA"
    if x == 0:
        return "<0.001"
    if x < 1e-3:
        return f"{x:.2g}"
    return f"{x:.3f}"

def num(x, nd=3):
    try:
        v = float(x)
    except Exception:
        return str(x)
    if v != v:
        return "NA"
    s = f"{v:.{nd}f}"
    if nd > 0:
        s = s.rstrip("0").rstrip(".")
    if s in ("-0", ""):
        s = "0"
    return s

def sci(x):
    try:
        return f"{float(x):.3g}"
    except Exception:
        return str(x)

def mdtable(df, headers):
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(v) for v in r.tolist()) + " |")
    return "\n".join(lines)

out = []
A = out.append
def cap(n, t):
    A(f"**Table S{n}.** {t}")
    A("")

# ---- S1: candidate genes ----
d = pd.read_csv(os.path.join(RES, "candidate_RBC_SAH_genes.csv"))
d = d[["gene", "logFC", "P.Value", "adj.P.Val", "source"]].copy()
d["logFC"] = d["logFC"].apply(lambda v: num(v, 2))
d["P.Value"] = d["P.Value"].apply(pval)
d["adj.P.Val"] = d["adj.P.Val"].apply(pval)
d["source"] = d["source"].fillna("-")
d = d.sort_values("adj.P.Val")
cap(1, "The 32 candidate genes at the intersection of SAH differentially expressed genes and the red-cell/HRR signature.")
A(mdtable(d, ["Gene", "log2FC", "P", "adj. P", "Signature source"])); A("")

# ---- S2: hub genes ----
d = pd.read_csv(os.path.join(RES, "hub_genes_annotated.csv"))
d = d[["gene", "logFC", "adj.P.Val", "votes"]].copy()
d["logFC"] = d["logFC"].apply(lambda v: num(v, 2))
d["adj.P.Val"] = d["adj.P.Val"].apply(pval)
cap(2, "Consensus machine-learning hub genes (top-ranked by at least three of LASSO / random-forest / SVM-RFE criteria).")
A(mdtable(d, ["Gene", "log2FC", "adj. P", "Ranking votes"])); A("")

# ---- S3: GSEA ----
d = pd.read_csv(os.path.join(RES, "GSEA_GSE36791_all.csv"))
keep = ["H:heme Metabolism", "H:Reactive Oxygen Species Pathway", "H:TNF-alpha Signaling via NF-kB",
        "H:Hypoxia", "H:Ferroptosis", "H:Complement", "K:Complement and coagulation cascades",
        "R:Neutrophil Degranulation R-HSA-6798695", "H:IL6 JAK STAT3 Signaling", "H:Ribosome"]
sub = d[d.Term.isin(keep)].copy() if any(k in set(d.Term) for k in keep) else d.head(12)
sub["pathway"] = sub.Term.apply(lambda t: re.sub(r"^(R|K|H|C|GOBP):", "", str(t))
                                .replace("R-HSA-", "").strip()).str.replace(r"\s*\d{6,}\s*$", "", regex=True)
sub = sub[["pathway", "NES", "NOM p-val", "FDR q-val"]].copy()
sub["NES"] = sub["NES"].apply(lambda v: num(v, 2))
sub["NOM p-val"] = sub["NOM p-val"].apply(pval)
sub["FDR q-val"] = sub["FDR q-val"].apply(pval)
sub = sub.sort_values("NES", key=lambda s: pd.to_numeric(s, errors="coerce"), ascending=False).head(12)
cap(3, "Gene-set enrichment analysis (preranked) of the discovery cohort, selected pathways.")
A(mdtable(sub, ["Pathway", "NES", "Nominal P", "FDR q"])); A("")

# ---- S4: WGCNA module-trait ----
d = pd.read_csv(os.path.join(RES, "wgcna_module_trait.csv"), index_col=0)
d = d[["n_genes", "cor_SAH", "p_SAH"]].copy()
d["Module"] = [str(i) for i in d.index]
d["cor_SAH"] = d["cor_SAH"].apply(lambda v: num(v, 2))
d["p_SAH"] = d["p_SAH"].apply(pval)
d = d.sort_values("cor_SAH", key=lambda s: pd.to_numeric(s, errors="coerce").abs(), ascending=False).head(10)
d = d.reset_index(drop=True)[["Module", "n_genes", "cor_SAH", "p_SAH"]]
cap(4, "WGCNA module-eigengene correlation with SAH status (top modules).")
A(mdtable(d, ["Module", "Genes", "r with SAH", "P"])); A("")

# ---- S5: real MR HP->SAH ----
r = pd.read_csv(os.path.join(RES, "mr_real_HP_SAH_results.csv")).iloc[0]
rows = [
    ("Instruments retained", num(r["n_instruments"], 0)),
    ("IVW beta (log OR per SD)", num(r["IVW_beta"], 3)),
    ("IVW SE", num(r["IVW_se"], 3)),
    ("IVW P", pval(r["IVW_p"])),
    ("IVW OR (95% CI)", f"{sci(math.exp(r['IVW_beta']))} [{sci(math.exp(r['IVW_beta']-1.96*r['IVW_se']))}, {sci(math.exp(r['IVW_beta']+1.96*r['IVW_se']))}]"),
    ("Cochran Q (df), P", f"{num(r['IVW_Q'],2)} ({num(r['IVW_df'],0)}), {pval(chi2.sf(r['IVW_Q'], r['IVW_df']))}"),
    ("MR-Egger beta, P", f"{num(r['Egger_beta'],3)}, {pval(r['Egger_p'])}"),
    ("MR-Egger intercept P", pval(r["Egger_intercept_p"])),
    ("Weighted-median beta, P", f"{num(r['WMedian_beta'],3)}, {pval(r['WMedian_p'])}"),
    ("Mean / min instrument F", f"{num(r['F_stat_mean'],0)} / {num(r['F_stat_min'],0)}"),
]
cap(5, "Two-sample Mendelian randomization: haptoglobin plasma level (GCST90019428) on SAH (GCST90473583).")
A("| Estimate | Value |"); A("|---|---|")
for k, v in rows:
    A(f"| {k} | {v} |")
A("")

# ---- S6: pQTL-MR panel (executed only) ----
d = pd.read_csv(os.path.join(RES, "mr_panel_HDFI_axis.csv"))
d = d[d["status"].astype(str).str.lower() == "executed"].copy()
d = d[["exposure", "outcome", "instruments_retained", "IVW_OR", "IVW_OR_CI_low", "IVW_OR_CI_high", "IVW_p"]].copy()
d["instruments_retained"] = d["instruments_retained"].apply(lambda v: num(v, 0))
d["IVW_OR"] = d["IVW_OR"].apply(sci)
d["ci"] = d.apply(lambda r: f"[{sci(r['IVW_OR_CI_low'])}, {sci(r['IVW_OR_CI_high'])}]", axis=1)
d["IVW_p"] = d["IVW_p"].apply(pval)
d = d[["exposure", "outcome", "instruments_retained", "IVW_OR", "ci", "IVW_p"]]
cap(6, "Mendelian-randomization panel: heme-axis pQTLs on SAH and intracerebral hemorrhage (IVW odds ratios).")
A(mdtable(d, ["Exposure", "Outcome", "IVs", "OR", "95% CI", "P"])); A("")

# ---- S7: cell-type-specific DE ----
d = pd.read_csv(os.path.join(RES, "scrna_celltype_DE_panels.csv"))
sel = {"Neutrophil": ["Hp", "Hmox1", "Fth1", "Ftl1", "Lcn2", "Mmp9", "Cd163", "Mrc1"],
       "Microglia/Macrophage": ["Hmox1", "Fth1", "Cd163", "Mrc1", "P2ry12"],
       "Neuron": ["Gpx4", "Acsl4", "Slc7a11", "Ptgs2", "Hmox1", "Slc40a1"]}
def keep_row(r):
    return r["gene"] in sel.get(r["cell_type"], [])
sub = d[d.apply(keep_row, axis=1)][["cell_type", "gene", "log2FC", "p_adj"]].copy()
sub["log2FC"] = sub["log2FC"].apply(lambda v: num(v, 2))
sub["p_adj"] = sub["p_adj"].apply(pval)
cap(7, "Cell-type-specific differential expression (SAH vs Sham), selected heme-handling and ferroptosis genes.")
A(mdtable(sub, ["Cell type", "Gene", "log2FC", "adj. P"])); A("")

# ---- S8: regulon ----
d = pd.read_csv(os.path.join(RES, "regulon_nrf2_bach1_by_celltype.csv"))
d = d[["cell_type", "metric", "delta", "cohend", "p"]].copy()
d["delta"] = d["delta"].apply(lambda v: num(v, 2))
d["cohend"] = d["cohend"].apply(lambda v: num(v, 2))
d["p"] = d["p"].apply(pval)
cap(8, "Regulon and heme-detoxification-failure activity, SAH vs Sham, by cell type.")
A(mdtable(d, ["Cell type", "Metric", "Delta", "Cohen d", "P"])); A("")

# ---- S9: composition-adjusted ----
d = pd.read_csv(os.path.join(RES, "composition_adjusted_de.csv"))
d = d[["outcome", "p_unadj", "p_adj", "q_adj"]].copy()
for c in ["p_unadj", "p_adj", "q_adj"]:
    d[c] = d[c].apply(pval)
cap(9, "Whole-blood heme outcomes before and after adjustment for neutrophil, monocyte, and T pseudo-fractions.")
A(mdtable(d, ["Outcome", "P (unadjusted)", "P (adjusted)", "q (adjusted)"])); A("")

# ---- S10: CSF proteome ----
d = pd.read_csv(os.path.join(TAB, "proteomics_hdfi.csv"))
d = d[d["detected"].astype(str).str.lower().isin(["true", "1"])].copy()
d = d[~d["role"].astype(str).str.contains("token count", case=False, na=False)]
d = d[["source", "protein", "role", "log2FC_sah_vs_ctrl", "mannwhitney_p"]].copy()
d["log2FC_sah_vs_ctrl"] = d["log2FC_sah_vs_ctrl"].apply(lambda v: num(v, 2))
d["mannwhitney_p"] = d["mannwhitney_p"].apply(pval)
cap(10, "Human SAH cerebrospinal-fluid proteome: axis proteins detected, SAH vs control.")
A(mdtable(d, ["Source", "Protein", "Role", "log2FC", "Mann-Whitney P"])); A("")

# ---- S11: within-CSF injury correlations ----
m = pd.read_csv(os.path.join(TAB, "csf_axis_injury_matrix.csv"), index_col=0)
if "INJURY_SCORE" in m.columns:
    rl = pd.DataFrame({"Protein": ["HP", "CD163", "FTH1", "FTL", "TFRC"]})
    rr, pp = [], []
    for p in rl["Protein"]:
        rho, pv = spearmanr(m[p].values, m["INJURY_SCORE"].values)
        rr.append(num(rho, 2)); pp.append(pval(pv))
    rl["rho vs injury score"] = rr
    rl["P"] = pp
    cap(11, "Within-CSF association of axis proteins with the neuronal-injury protein score (PXD030593, n = 8 samples).")
    A(mdtable(rl, ["Protein", "rho vs injury score", "P"])); A("")

# ---- S12: ferroptosis execution lipids ----
d = pd.read_csv(os.path.join(RES, "lipidomics_ferroptosis.csv"))
cols = [c for c in ["analyte", "dataset", "marker_class", "fold_change", "MWU_p", "cliffs_delta"] if c in d.columns]
d2 = d[cols].copy()
if "fold_change" in d2: d2["fold_change"] = d2["fold_change"].apply(sci)
if "MWU_p" in d2: d2["MWU_p"] = d2["MWU_p"].apply(pval)
if "cliffs_delta" in d2: d2["cliffs_delta"] = d2["cliffs_delta"].apply(lambda v: num(v, 2))
d2 = d2.head(14)
cap(12, "Oxidised-PE ferroptosis execution markers across two public targeted metabolomics datasets.")
A(mdtable(d2, ["Analyte", "Dataset", "Class", "Fold change", "P", "Cliff's delta"])); A("")

# ---- S13: HDFI robustness / leakage ----
s = json.load(open(os.path.join(RES, "hdfi_robustness_summary.json")))
lf = s["leakage_free"]; rn = s["random_null"]; pca = s["pca"]
logo = pd.read_csv(os.path.join(RES, "hdfi_robustness_logo.csv"))
A("**Table S13.** HDFI robustness and GPX4-leakage audit on the neuronal association (raw-count-z pipeline).")
A(""); A("| Test | Result |"); A("|---|---|")
A(f"| Canonical neuronal association | rho {num(lf['rho_canon_vs_sharedGPX4_outcome'],3)}, P {pval(lf['p0'])} |")
A(f"| Leakage-free (no shared gene) | rho {num(lf['rho_GPX4free_HDFI_vs_GPX4free_outcome'],3)}, P {pval(lf['p2'])} |")
A(f"| GPX4 removed from outcome only | rho {num(lf['rho_canon_vs_GPX4free_outcome'],3)}, P {pval(lf['p1'])} |")
A(f"| Sign-reversed index | rho {num(rn['sign_reversed_rho'],3)} |")
A(f"| Random-weighting null (n={int(rn['n_random'])}) | one-sided empirical P {num(rn['empirical_p_one_sided'],3)}; best random rho {num(rn['random_max'],3)} |")
A(f"| PCA-independence | partial rho {num(pca['partial_rho_given_PC1'],3)}, P {pval(pca['partial_p'])} |")
for _, r in logo.iterrows():
    A(f"| Leave-one-out (drop {r['dropped']}) | rho {num(r['rho'],3)}, P {pval(r['p'])} |")
A("")

# ---- S14: accession inventory ----
cap(14, "Public data-accession inventory.")
A("| Type | Accession | Content |"); A("|---|---|---|")
for t, a, c in [
    ("GEO", "GSE36791", "acute SAH whole blood (GPL10558); discovery"),
    ("GEO", "GSE73378", "remote/past aSAH whole blood (GPL10558); test"),
    ("GEO", "GSE311852", "mouse SAH single-cell brain atlas (Sham/24 h/72 h)"),
    ("GEO", "GSE163256", "ICH blood monocytes (external validation)"),
    ("GEO", "GSE296792", "ICH whole blood (external validation)"),
    ("GEO", "GSE269478", "human monocyte RNA+ADT CITE-seq"),
    ("GEO", "GSE172116", "human macrophage bulk ATAC-seq (LPS contrast)"),
    ("GEO", "GSE254161", "mouse cortical-injury 10x Visium (spatial)"),
    ("PRIDE", "PXD030593", "human SAH-CSF label-free DDA proteome"),
    ("GWAS", "GCST90019428", "haptoglobin plasma pQTL (Eur, N=10,708)"),
    ("GWAS", "GCST90473583", "SAH GWAS (Eur, N=458,440)"),
    ("GWAS", "GCST90473584", "ICH GWAS (Non-Finnish Eur, N=458,440)"),
    ("GWAS", "GCST90270865 / GCST90451683 / GCST90469465 / GCST90470840", "ferritin / hepcidin / HMOX1 / TFRC pQTLs"),
    ("ENCODE", "ENCSR584GHV / ENCSR707IUN / ENCSR488EES / ENCSR197WGI", "NFE2L2 ChIP-seq (A549/HeLa/HepG2/IMR-90)"),
    ("JASPAR", "MA0150.3", "NFE2L2 PWM"),
    ("Metabolomics-Workbench", "ST002968 / ST003337 / ST004928", "stroke plasma / MEF / MDSC lipid panels"),
]:
    A(f"| {t} | {a} | {c} |")
A("")

with open(OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(out))
print("wrote", OUT)
