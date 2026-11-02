# sah-hdfi-axis

Reproduction code for the manuscript **"A heme-detoxification-failure state links
red-cell dysfunction to ferroptotic vulnerability after subarachnoid hemorrhage"**
(integrative multi-layer analysis of the Heme-Detoxification-Failure Index, HDFI).

This repository contains the **analysis reproduction code**, the **derived result
tables**, and the **publication figures**. Raw source data are **not** redistributed;
every dataset is a public accession (GEO / PRIDE / GWAS Catalog / Metabolomics
Workbench / ENCODE / JASPAR) listed in `results/` and in the manuscript Data
availability section, and can be re-fetched with `scripts/00_download.py` and
`scripts/10_scrna_download.py`.

## Layout

```
scripts/        numbered analysis pipeline (00 -> 51) + naturestyle.py + build_* helpers
  00_download.py, 01_prepare.py     bulk cohort download / prepare
  04_de.py, 05_signature.py, 06_ml.py, 07_gsea.py, 08_score.py, 09_figures.py
  11_wgcna.py, 12_scrna.py, 13_ml_xgb.py, 14_mr.py, 15_newfigures.py
  16_scrna_deep.py, 17_scrna_deep_figures.py, 18_scrna_temporal.py, 19_bulk_deconvolution.py
  20_graphical_abstract.py, 22_regulon_bach1_nrf2.py, 23_cellchat_style.py, 24_ferroptosis_gsva.py
  26_hdfi_centrepiece.py, 27_composition_corrected_de.py, 28_cellstate.py, 30_regulatory_motif.py
  31_mr_real.py, 32_mr_forest.py, 33_hdfi_external.py, 34_mr_panel.py, 35_nichenet_style.py
  36_scatac_nrf2.py, 37_lipidomics_stroke_panel.py, 38_brain_spatial_coloc.py, 39_drugrepurposing_targets.py
  41_cite_seq_cd163.py, 42_scatac_hmox1.py, 44_lipidomics_ferroptosis.py, 45_drug_reversal_lincs.py
  46_hdfi_latent_axis.py, 47_proteomics_hdfi.py, 48_sah_csf_lipidomics.py, 49_mouse_sah_atlas.py
  50_csf_neuronal_injury_link.py, 51_hdfi_robustness.py
  naturestyle.py                    shared publication-figure style (Arial + Okabe-Ito,
                                    in-figure notes stripped, composites split into
                                    independent panels, 5 cm footprint)
  build_supp_tables.py              renders Supplementary Tables from results/
  build_docx.py                     renders the manuscript .md to a submission .docx
results/        derived result tables (CSV / JSON) and GWAS instrument lists
figures/        publication figures (PNG + editable PDF), one file per panel
manuscript/     the manuscript and supplementary markdown
```

## Reproduce

```bash
# 1. fetch public data
python scripts/00_download.py
python scripts/10_scrna_download.py

# 2. run the pipeline in numeric order
for s in scripts/0*.py scripts/1*.py scripts/2*.py scripts/3*.py scripts/4*.py scripts/5*.py; do
    python "$s"
done

# 3. build supplementary tables and the manuscript docx
python scripts/build_supp_tables.py
python scripts/build_docx.py manuscript/MAIN_submission.md MAIN.docx main
python scripts/build_docx.py manuscript/SUPP_submission.md SUPP.docx supp
```

## Environment

Python 3.10 with pandas, numpy, scipy, statsmodels, scikit-learn, gseapy,
scanpy / anndata, xgboost, matplotlib, openpyxl, requests. Random seeds fixed (2026).
The regulatory layer uses live Ensembl + JASPAR REST endpoints; the MR layer reads
GWAS-Catalog harmonised summary statistics.

## Notes on scope

Negative, composition-confounded and access-limited results are reported as-is; the
only simulations are clearly-labelled engine-correctness checks (MR, coloc, mediation).
See the manuscript Limitations section for the boundaries of each claim.
