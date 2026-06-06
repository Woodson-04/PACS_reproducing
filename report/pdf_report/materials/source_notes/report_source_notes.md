# Report Source Notes

This note records where each report material came from and whether it was directly copied, summarized from an existing CSV/report, or missing.

## Analysis Boundary

- No new analysis was run in this material-pack step.
- No files under `/home/woodson/biostatistic` were modified.
- No installed PACS package files were modified.
- No git operation was performed.

## Figure Sources

Binary figure copying failed in the current Windows SSHFS/PowerShell session with `Access to the path ... is denied`. The source figures below were inspected as existing files, but were not copied into `report/pdf_report/materials/figures/` during this turn.

| target name | source | status |
|---|---|---|
| `fig_umap_four_panel_p56_10000.png` | `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png` | source exists; copy failed due session permission |
| `fig_p56_batch_mixing_metrics_10000.png` | `figures/mouse_kidney/gse157079_p56_batch_mixing_p56_top10000_test10000_fdr005_metrics_barplot.png` | source exists; copy failed due session permission |
| `fig_pca_lsi_umap_score_comparison.png` | `figures/mouse_kidney/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png` | source exists; copy failed due session permission |
| `fig_lsi_depth_correlation_before_after.png` | `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before_after.png` | source exists; copy failed due session permission |
| `fig_lsi_dimension_sensitivity.png` | `figures/mouse_kidney/gse157079_p56_lsi_dimension_sensitivity_batch_mixing.png` | source exists; copy failed due session permission |
| `fig_parameter_comparison_metrics.png` | `figures/mouse_kidney/gse157079_p56_pacs_setting_comparison_metrics.png` | source missing |
| `fig_geo_vs_pacs_reference.png` | `figures/mouse_kidney/gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png` | source exists; copy failed due session permission |
| `fig_all_cells_top_peaks_umap_by_sample.png` | `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png` | source exists; copy failed due session permission |
| `fig_all_cells_top_peaks_umap_by_celltype.png` | `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png` | source exists; copy failed due session permission |

## Table Sources

| output table | source and note |
|---|---|
| `table_p56_main_result_summary.csv` | summarized from `p56_pacs_batch_filter_umap_report.md` and `p56_batch_peak_summary.csv` in the main `lsi_saved` result directory |
| `table_pca_lognorm_scores.csv` | read from `results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space/p56_10000_pca_space_scores.csv`; fold change was calculated as after/before from existing values |
| `table_lsi_scores.csv` | read from `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_dimension_sensitivity_scores.csv`; improvement was calculated as after minus before from existing values |
| `table_lsi_depth_correlation.csv` | read from `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_depth_correlation.csv`, using before/after `LSI_1` rows |
| `table_parameter_comparison.csv` | assembled from existing `p56_batch_peak_summary.csv` and `p56_batch_mixing_metrics_summary.csv` files for 5000/5000, 10000/10000, and 20000/10000 settings |
| `table_geo_reference_comparison.csv` | assembled from `paper_style_batch_mixing_scores.csv` and `geo_vs_p56_pacs_umap_metrics_comparison.csv` |

## Values That Should Be Highlighted In The Report

- Main result: 13,526 P56 cells, 10,000 top peaks, 10,000 PACS-tested peaks.
- PACS found 6,305 significant batch-associated peaks at FDR 0.05.
- Retained peaks after filtering: 3,695.
- PCA-logNorm normalized mixing improved:
  - PC1:20: 0.0398 to 0.2576
  - PC1:30: 0.0510 to 0.3097
  - PC1:50: 0.0685 to 0.4239
- LSI-space normalized mixing improved:
  - LSI_2:30: 0.0520 to 0.3343
  - LSI_2:50: 0.0793 to 0.3730
- UMAP-space normalized mixing improved:
  - 0.0265 to 0.6507
- LSI_1 is strongly depth-associated:
  - before Spearman = -0.9983
  - after Spearman = 0.9976

## Interpretation Boundary

The current results support that PACS filtering removes a substantial portion of batch-associated feature signal while preserving visible biological cell-type structure. They are not an exact reproduction of the full PACS paper figure because this is a P56-only, top-peak, derived-cap-rate reconstruction rather than the author's full analysis with all original preprocessing settings.
