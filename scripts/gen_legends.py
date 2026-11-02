"""gen_legends.py -- build figure legends that match the ACTUAL split panels.

Reads the panel files produced by the render, assigns subpanel letters a,b,c... per
main/supp figure (in manifest order), groups consecutive panels that come from the
same source analysis, and emits legend lines with a descriptor per source. Guarantees
subpanel-level numbering with no gaps.
"""
import os, re, glob

FIG = "figures"
MAIN = {
    1: ["Fig10_axis_schematic"],
    2: ["Fig1_volcano_GSE36791", "Fig3_heatmap_candidates", "Fig4_enrichment",
        "Fig14_bulk_deconvolution", "Fig19_composition_adjusted"],
    3: ["Fig8_scrna_summary", "Fig12_scrna_deep", "Fig13_scrna_temporal", "Fig20_cellstate"],
    4: ["Fig16_cellchat_style", "Fig25_nichenet_style", "Fig34_cite_seq_cd163"],
    5: ["Fig15_regulon_bach1_nrf2", "Fig18_hdfi_centrepiece", "Fig36_hdfi_latent_axis",
        "Fig41_hdfi_robustness"],
    6: ["Fig21_regulatory_motif", "Fig26_nrf2_chip_binding", "Fig32_scatac_hmox1"],
    7: ["Fig33_lipidomics_ferroptosis", "Fig37_proteomics_hdfi",
        "Fig40_csf_hdfi_neuronal_injury", "Fig23_hdfi_external"],
    8: ["Fig0_graphical_abstract"],
}
SUPP = {
    1: ["Fig7_wgcna"], 2: ["Fig5_ml_ranking"], 3: ["Fig6_score_and_validation"],
    4: ["Fig9_xgb"], 5: ["umap_celltype", "umap_heme", "dotplot__markers"],
    6: ["Fig39_mouse_sah_atlas"], 7: ["Fig11_scrna_trajectory"],
    8: ["Fig17_ferroptosis_gsva"], 9: ["Fig29_lipidomics_stroke"],
    10: ["Fig38_sah_csf_lipidomics"], 11: ["Fig22_mr_forest"], 12: ["Fig24_mr_panel"],
    13: ["Fig27_spatial_coloc"], 14: ["Fig28_drugrepurposing_targets"],
    15: ["Fig35_drug_reversal_lincs"],
}
# descriptor per SOURCE figure (one phrase; may cover several split panels)
DESC = {
    "Fig10_axis_schematic": "overview schematic of the HRR→heme-detoxification→ferroptosis axis and the multi-layer study design",
    "Fig1_volcano_GSE36791": "volcano plot of differentially expressed genes in acute SAH blood (GSE36791)",
    "Fig3_heatmap_candidates": "heatmap of the 32 candidate genes at the intersection of SAH differentially expressed genes and the red-cell signature",
    "Fig4_enrichment": "over-representation enrichment of the candidate set and Hallmark/Reactome gene-set enrichment scores",
    "Fig14_bulk_deconvolution": "marker-based immune deconvolution of the neutrophilia/monocytosis/lymphopenia shift and the heme/iron program",
    "Fig19_composition_adjusted": "composition-adjusted differential expression showing heme genes lose significance after leukocyte-fraction adjustment",
    "Fig8_scrna_summary": "heme/iron module score by cell type (SAH versus Sham) and marker localisation in the mouse SAH atlas",
    "Fig12_scrna_deep": "cell-type-specific differential expression, the neuronal ferroptosis program, and curated regulon activity",
    "Fig13_scrna_temporal": "temporal composition and heme/iron and ferroptosis scores across Sham, 24 h and 72 h",
    "Fig20_cellstate": "macrophage and neutrophil cell-state transitions between homeostatic, iron-handling and ferroptosis-prone states",
    "Fig16_cellchat_style": "CellChat-style ligand–receptor test (HP–CD163, HP–MRC1, DAMP routes) per timepoint",
    "Fig25_nichenet_style": "NicheNet-style ligand-to-target prioritisation of the microglial heme-handling program",
    "Fig34_cite_seq_cd163": "monocyte CITE-seq showing CD163 surface protein where its RNA dropped out, and MRC1 levels",
    "Fig15_regulon_bach1_nrf2": "NRF2/BACH1/HMOX1 regulon activity and the heme-detoxification-failure index by cell type",
    "Fig18_hdfi_centrepiece": "HDFI in bulk blood, by cell type, versus neuronal ferroptosis-execution within neurons, and its temporal behaviour",
    "Fig36_hdfi_latent_axis": "pathway-specificity and PCA factor-structure audit of HDFI",
    "Fig41_hdfi_robustness": "HDFI structural-variant, leave-one-gene-out, GPX4-leakage, random-weighting and PCA-independence robustness tests",
    "Fig21_regulatory_motif": "ARE/NFE2L2 promoter motif density and ENCODE regulatory-feature counts across the heme-handling genes",
    "Fig26_nrf2_chip_binding": "ENCODE NFE2L2 ChIP-seq enrichment at the HMOX1 locus across four cell lines",
    "Fig32_scatac_hmox1": "macrophage ATAC-seq showing constitutive HMOX1 openness and NRF2/ARE enrichment on accessible peaks",
    "Fig33_lipidomics_ferroptosis": "oxidised arachidonoyl-PE ferroptosis execution markers across two targeted metabolomics models",
    "Fig37_proteomics_hdfi": "human SAH-CSF proteome axis-protein detection and the independent-cohort cross-check",
    "Fig40_csf_hdfi_neuronal_injury": "within-CSF association of the iron-storage arm with a neuronal-injury protein score, with leave-one-out and blood-contamination controls",
    "Fig23_hdfi_external": "HDFI external validation in independent intracerebral-hemorrhage blood cohorts (null)",
    "Fig0_graphical_abstract": "integrated model of the heme-detoxification-failure state (every arrow labelled associated-with, not causes)",
    "Fig7_wgcna": "signed WGCNA network with SAH-correlated module and gene-significance / module-membership for heme–iron genes",
    "Fig5_ml_ranking": "consensus machine-learning hub selection across LASSO, random-forest and SVM-RFE",
    "Fig6_score_and_validation": "molecular-score internal cross-validation and external transfer to the remote-aSAH cohort",
    "Fig9_xgb": "gradient-boosting comparison of four feature sets (internal versus external AUC)",
    "Fig11_scrna_trajectory": "microglia/macrophage diffusion-pseudotime trajectory (homeostatic-rooted)",
    "Fig17_ferroptosis_gsva": "FerrDb-categorised ferroptosis Driver/Suppressor/Iron and execution scores",
    "Fig22_mr_forest": "two-sample Mendelian-randomization of haptoglobin on SAH (instrument and summary forest)",
    "Fig24_mr_panel": "Mendelian-randomization panel of heme-axis pQTLs on SAH and intracerebral hemorrhage",
    "Fig27_spatial_coloc": "spatial co-localisation test on an independent brain-injury Visium dataset (null)",
    "Fig28_drugrepurposing_targets": "drug-target availability shortlist for the axis proteins (OpenTargets)",
    "Fig29_lipidomics_stroke": "rat post-stroke plasma lipidomics (arachidonate and related metabolites)",
    "Fig35_drug_reversal_lincs": "LINCS L1000 direction-aware signature-reversal ranking",
    "Fig38_sah_csf_lipidomics": "human SAH-CSF lipidome PUFA/ether-PE elevations (reported verbatim from the source study)",
    "Fig39_mouse_sah_atlas": "mouse SAH immune-atlas HDFI program across cell types and timepoints",
    "umap_celltype": "UMAP of the SAH brain atlas coloured by cell type",
    "umap_heme": "UMAP coloured by the heme/iron module score",
    "dotplot__markers": "canonical-marker dotplot across cell types",
}

def panels(b):
    hits = sorted(glob.glob(os.path.join(FIG, b + "*.png")))
    sp = [h for h in hits if re.search(r"_[a-z]\.png$", h)]
    return sp if sp else [h for h in hits if not re.search(r"_[a-z]\.png$", h)]

def build(man, prefix, titles):
    out = []
    for n in sorted(man):
        seq = []  # (base, [letters])
        idx = 0
        for b in man[n]:
            pf = panels(b)
            letters = []
            for _ in pf:
                letters.append(chr(ord("a") + idx)); idx += 1
            if letters:
                seq.append((b, letters))
        rng = lambda L: L[0].upper() if len(L) == 1 else f"{L[0].upper()}–{L[-1].upper()}"
        parts = "; ".join(f"({rng(L)}) {DESC.get(b, b)}" for b, L in seq)
        out.append(f"**Figure {prefix}{n}. {titles.get(n,'')}** {parts}.")
    return out

MAIN_TITLES = {
    1: "Overview and study design.", 2: "The acute whole-blood heme–immune program is largely composition-driven.",
    3: "Single-cell brain atlas localises heme handling to myeloid cells and ferroptotic vulnerability to neurons.",
    4: "The HP→CD163 model refines into a broader heme-handling program.",
    5: "HDFI as a conditional detoxification-failure state, with its dependence disclosed.",
    6: "Coordinated NRF2/ARE control of the heme-handling program.",
    7: "Real execution lipids; protein and external-validation limits.",
    8: "Integrated model.",
}
SUPP_TITLES = {
    1: "Co-expression network.", 2: "Machine-learning hub selection.",
    3: "Molecular-score validation.", 4: "Gradient-boosting comparison.",
    5: "Single-cell reference panels.", 6: "Mouse SAH immune-atlas HDFI program.",
    7: "Pseudotime trajectory (negative).", 8: "FerrDb ferroptosis scoring.",
    9: "Stroke plasma lipidomics.", 10: "Human SAH-CSF lipidome (verbatim).",
    11: "Mendelian-randomization of haptoglobin on SAH.", 12: "pQTL–Mendelian-randomization panel.",
    13: "Spatial co-localisation (null).", 14: "Drug-target availability.",
    15: "LINCS signature reversal.",
}

if __name__ == "__main__":
    os.makedirs("manuscript", exist_ok=True)
    with open(os.path.join("manuscript", "LEGENDS_main.md"), "w", encoding="utf-8") as f:
        f.write("\n\n".join(build(MAIN, "", MAIN_TITLES)) + "\n")
    with open(os.path.join("manuscript", "LEGENDS_supp.md"), "w", encoding="utf-8") as f:
        f.write("\n\n".join(build(SUPP, "S", SUPP_TITLES)) + "\n")
    print("wrote manuscript/LEGENDS_main.md and LEGENDS_supp.md")
