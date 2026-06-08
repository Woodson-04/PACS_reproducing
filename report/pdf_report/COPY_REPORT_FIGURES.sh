#!/usr/bin/env bash
set -u

mkdir -p report/pdf_report/materials/figures

copy_one() {
  local src="$1"
  local dst_dir="report/pdf_report/materials/figures"
  if [ -f "$src" ]; then
    cp -v "$src" "$dst_dir/"
  else
    echo "WARNING: missing source figure: $src" >&2
  fi
}

copy_one "figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png"
copy_one "figures/mouse_kidney/pacs_permuted_qq_plot.png"
copy_one "figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png"
copy_one "figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png"
copy_one "figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png"
copy_one "figures/mouse_kidney/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png"
copy_one "figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before_after.png"
copy_one "figures/mouse_kidney/gse157079_p56_lsi_dimension_sensitivity_batch_mixing.png"
copy_one "figures/mouse_kidney/gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png"

echo "Figure copy step finished."
