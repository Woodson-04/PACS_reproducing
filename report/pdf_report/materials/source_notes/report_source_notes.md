# Report Source Notes

This file maps the final report materials to their source files.

## Tables

| table | output file | source |
|---|---|---|
| Notebook 1 benchmark | `materials/tables/table1_notebook1_benchmark.md` | `results/20260526_2318_large_baseline/summary.csv` |
| GSE157079 alignment | `materials/tables/table2_gse157079_alignment.md` | `results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md` |
| all-cell UMAP | `materials/tables/table3_all_cell_umap.md` | `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md` |
| P56 PACS filtering | `materials/tables/table4_p56_pacs_filtering.md` | `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_pacs_batch_filter_umap_report.md` |
| normalized batch mixing | `materials/tables/table5_normalized_batch_mixing.md` | PCA, LSI, and UMAP score CSVs under `results/mouse_kidney_figures/` |
| LSI-depth QC | `materials/tables/table6_lsi_depth_qc.md` | `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_depth_correlation_report.md` |
| parameter stability | `materials/tables/table7_parameter_stability.md` | `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_comparison/p56_pacs_filter_setting_comparison.csv` plus current 20000/10000 result report |

## Figures

| figure in final report | intended copied file | original source |
|---|---|---|
| Notebook 1 benchmark | `materials/figures/pacs_benchmark_t1e_power_barplot.png` | `figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png` |
| Permuted-label QQ plot | `materials/figures/pacs_permuted_qq_plot.png` | `figures/mouse_kidney/pacs_permuted_qq_plot.png` |
| all-cell UMAP by sample | `materials/figures/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png` | `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png` |
| all-cell UMAP by cell type | `materials/figures/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png` | `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png` |
| P56 PACS four-panel UMAP | `materials/figures/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png` | `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png` |
| PCA/LSI/UMAP score comparison | `materials/figures/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png` | `figures/mouse_kidney/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png` |
| LSI-depth QC | `materials/figures/gse157079_p56_lsi_depth_correlation_before_after.png` | `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before_after.png` |
| LSI dimension sensitivity | `materials/figures/gse157079_p56_lsi_dimension_sensitivity_batch_mixing.png` | `figures/mouse_kidney/gse157079_p56_lsi_dimension_sensitivity_batch_mixing.png` |
| GEO vs P56 PACS UMAP metrics | `materials/figures/gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png` | `figures/mouse_kidney/gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png` |
