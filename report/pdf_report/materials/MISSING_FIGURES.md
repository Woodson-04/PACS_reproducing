# Missing Or Not-Copied Figures

The current Windows SSHFS/PowerShell session could create text files through `apply_patch`, but binary `Copy-Item` into `report/pdf_report/materials/figures/` failed with `Access to the path ... is denied`.

Therefore, source figures were not copied in this turn. This is a packaging limitation of the current session, not a new analysis result.

## Source Figures Found But Not Copied

- `fig_umap_four_panel_p56_10000.png`
  - source: `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png`
- `fig_p56_batch_mixing_metrics_10000.png`
  - source: `figures/mouse_kidney/gse157079_p56_batch_mixing_p56_top10000_test10000_fdr005_metrics_barplot.png`
- `fig_pca_lsi_umap_score_comparison.png`
  - source: `figures/mouse_kidney/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png`
- `fig_lsi_depth_correlation_before_after.png`
  - source: `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before_after.png`
- `fig_lsi_dimension_sensitivity.png`
  - source: `figures/mouse_kidney/gse157079_p56_lsi_dimension_sensitivity_batch_mixing.png`
- `fig_geo_vs_pacs_reference.png`
  - source: `figures/mouse_kidney/gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png`
- `fig_all_cells_top_peaks_umap_by_sample.png`
  - source: `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png`
- `fig_all_cells_top_peaks_umap_by_celltype.png`
  - source: `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png`

## Requested Source Figure Not Found

- `fig_parameter_comparison_metrics.png`
  - requested source: `figures/mouse_kidney/gse157079_p56_pacs_setting_comparison_metrics.png`
  - status: source figure was not present in the inspected figure directory.

## Linux Copy Commands

Run these from `/home/woodson/PACS_reproducing` if the figure material directory should be populated before PDF generation:

```bash
mkdir -p report/pdf_report/materials/figures
cp figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png report/pdf_report/materials/figures/fig_umap_four_panel_p56_10000.png
cp figures/mouse_kidney/gse157079_p56_batch_mixing_p56_top10000_test10000_fdr005_metrics_barplot.png report/pdf_report/materials/figures/fig_p56_batch_mixing_metrics_10000.png
cp figures/mouse_kidney/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png report/pdf_report/materials/figures/fig_pca_lsi_umap_score_comparison.png
cp figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before_after.png report/pdf_report/materials/figures/fig_lsi_depth_correlation_before_after.png
cp figures/mouse_kidney/gse157079_p56_lsi_dimension_sensitivity_batch_mixing.png report/pdf_report/materials/figures/fig_lsi_dimension_sensitivity.png
cp figures/mouse_kidney/gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png report/pdf_report/materials/figures/fig_geo_vs_pacs_reference.png
cp figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png report/pdf_report/materials/figures/fig_all_cells_top_peaks_umap_by_sample.png
cp figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png report/pdf_report/materials/figures/fig_all_cells_top_peaks_umap_by_celltype.png
```
