# Result Snapshot Manifest

Snapshot date/time: 2026-06-06

This snapshot freezes the current P56 PACS batch-effect filtering materials for a later detailed Chinese PDF report. No new PACS run, MatrixMarket streaming, TF-IDF/LSI/UMAP computation, or metric computation was performed in this packaging step.

## Main Run

- Main run name: `p56_top10000_test10000_fdr005_lsi_saved`
- Main result directory: `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/`
- Main figure set: `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_*`
- P56 cells: 13,526
- Batch counts: P56_batch1 = 7,129; P56_batch2 = 6,397
- Top peaks used: 10,000
- PACS-tested peaks: 10,000
- PACS-significant batch peaks at FDR 0.05: 6,305
- Retained peaks after filtering: 3,695

## Input Result Directories

- `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/`
- `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/`
- `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005/`
- `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap/`
- `results/mouse_kidney_figures/paper_style_batch_mixing_score/`
- `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/`
- `results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space/`
- `results/mouse_kidney_figures/geo_precomputed_umap_quantification/`
- `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/`

## Important Source CSV/Reports

- Main P56 report: `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_pacs_batch_filter_umap_report.md`
- Main batch peak summary: `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_batch_peak_summary.csv`
- Main PACS batch peak results: `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_pacs_batch_peak_results.csv`
- Main retained peaks: `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_retained_peak_indices.csv`
- Main before/after LSI embeddings: `p56_before_lsi_embedding.csv`, `p56_after_lsi_embedding.csv`
- Main before/after LSI UMAP embeddings: `p56_before_lsi_umap_embedding.csv`, `p56_after_lsi_umap_embedding.csv`
- UMAP-space paper-style scores: `results/mouse_kidney_figures/paper_style_batch_mixing_score/paper_style_batch_mixing_scores.csv`
- LSI-space paper-style scores: `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_10000_lsi_space_paper_style_scores.csv`
- LSI dimension sensitivity: `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_dimension_sensitivity_scores.csv`
- LSI-depth QC: `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_depth_correlation.csv`
- PCA-logNorm scores: `results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space/p56_10000_pca_space_scores.csv`
- GEO reference report: `results/mouse_kidney_figures/geo_precomputed_umap_quantification/geo_precomputed_umap_quantification_report.md`
- GEO vs PACS metrics CSV: `results/mouse_kidney_figures/geo_precomputed_umap_quantification/geo_vs_p56_pacs_umap_metrics_comparison.csv`
- All-cell top-peak UMAP report: `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md`

## Scripts Used In The Upstream Workflow

- `scripts/mouse_kidney_figures/07_p56_pacs_batch_filter_umap.R`
- `scripts/mouse_kidney_figures/07b_replot_p56_pacs_batch_filter_umap.R`
- `scripts/mouse_kidney_figures/07c_quantify_p56_batch_mixing.R`
- `scripts/mouse_kidney_figures/09_compute_paper_style_batch_mixing_score.R`
- `scripts/mouse_kidney_figures/09b_compare_paper_style_batch_mixing_scores.R`
- `scripts/mouse_kidney_figures/09c_compare_paper_style_lsi_batch_mixing_scores.R`
- `scripts/mouse_kidney_figures/09d_compute_pca_space_batch_mixing_score.R`
- `scripts/mouse_kidney_figures/09e_check_lsi_depth_correlation.R`
- `scripts/mouse_kidney_figures/09f_lsi_dimension_sensitivity_batch_mixing.R`

## Output Figures Intended For Report Materials

See `materials/copied_figures.csv` and `materials/MISSING_FIGURES.md`.

Most source figures exist in `figures/mouse_kidney/`, but binary copy failed in the current Windows SSHFS/PowerShell session due to permission denial. No source figures were modified.

## Output Tables Created In This Snapshot

- `report/pdf_report/materials/tables/table_p56_main_result_summary.csv`
- `report/pdf_report/materials/tables/table_pca_lognorm_scores.csv`
- `report/pdf_report/materials/tables/table_lsi_scores.csv`
- `report/pdf_report/materials/tables/table_lsi_depth_correlation.csv`
- `report/pdf_report/materials/tables/table_parameter_comparison.csv`
- `report/pdf_report/materials/tables/table_geo_reference_comparison.csv`

## No-New-Analysis Statement

This packaging step only read existing result files and prepared report materials. It did not rerun PACS, UMAP, PCA, LSI, nearest-neighbor metrics, or MatrixMarket streaming.
