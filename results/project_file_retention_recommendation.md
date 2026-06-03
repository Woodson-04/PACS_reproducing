# Project File Retention Recommendation

This is a recommendation only. No files have been deleted, moved, renamed, or
overwritten.

## A. Must Keep For Teacher Meeting

| relative path | reason | priority |
|---|---|---|
| `results/project_progress_summary_for_teacher.md` | concise meeting summary and current scientific story | high |
| `results/20260526_2318_large_baseline/summary.csv` | final Notebook 1 large benchmark table | high |
| `results/20260526_2318_large_baseline/pacs_kidney_notebook1_result.rds` | full Notebook 1 result object for verification | high |
| `results/20260526_2318_large_baseline/session_info.txt` | reproducibility record for benchmark run | medium |
| `q.r` | main Notebook 1 reproduction script, including local PACS wrapper fix | high |
| `notebook1_reproduction_report.md` | written conclusion for Notebook 1 stage | high |
| `figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png` | simple visual summary of benchmark performance | high |
| `figures/mouse_kidney/pacs_permuted_qq_plot.png` | shows PACS p-value behavior under permuted labels | high |
| `figures/mouse_kidney/gse157079_umap_by_sample.png` | GEO precomputed UMAP overview by sample | medium |
| `figures/mouse_kidney/gse157079_umap_by_celltype.png` | GEO precomputed UMAP overview by cell type | medium |
| `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png` | strongest current matrix-derived before-filtering UMAP by sample | high |
| `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png` | strongest current matrix-derived before-filtering UMAP by cell type | high |
| `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md` | documents all-cell top-peak matrix-derived UMAP run | high |
| `results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md` | proves GSE matrix, metadata, and peak list alignment | high |

## B. Keep But Not Show Unless Asked

| relative path | reason | priority |
|---|---|---|
| `baseline_methods_notebook1.R` | clean-room baseline implementations used in Notebook 1 | medium |
| `baseline_official_methods_review.md` | explains baseline uncertainty and method choices | medium |
| `mouse_kidney_figure_plan.md` | larger figure reproduction roadmap | medium |
| `mouse_kidney_figure_source_report.md` | earlier source/data discovery notes | medium |
| `results/mouse_kidney_figures/official_author_code_search_report.md` | author-code route investigation, useful if asked about dev-kidney-snATAC | medium |
| `results/mouse_kidney_figures/gse157079_metadata_merged.csv` | standardized metadata/UMAP table for scripts | high |
| `results/mouse_kidney_figures/gse157079_metadata_summary.csv` | quick metadata sample/cell type counts | medium |
| `results/mouse_kidney_figures/gse157079_peak_list_preview.csv` | lightweight peak-list preview | low |
| `results/mouse_kidney_figures/gse157079_inspection_report.md` | initial GSE file inspection | medium |
| `results/mouse_kidney_figures/gse157079_pilot_matrix_lsi_umap/pilot_lsi_umap_report.md` | sampled pilot run proof before all-cell run | medium |
| `figures/mouse_kidney/gse157079_pilot_matrix_lsi_umap_by_sample.png` | pilot UMAP, useful as intermediate proof | low |
| `figures/mouse_kidney/gse157079_pilot_matrix_lsi_umap_by_celltype.png` | pilot UMAP, useful as intermediate proof | low |
| `scripts/mouse_kidney_figures/00_inspect_gse157079.R` | GSE file inspection script | medium |
| `scripts/mouse_kidney_figures/00_prepare_gse157079_metadata.R` | metadata standardization script | high |
| `scripts/mouse_kidney_figures/01_overview_plots.R` | overview figure script | medium |
| `scripts/mouse_kidney_figures/03_gse157079_umap_plots.R` | GEO precomputed UMAP plotting script | medium |
| `scripts/mouse_kidney_figures/04_smoke_test_gse157079_matrix_alignment.R` | matrix alignment smoke-test script | high |
| `scripts/mouse_kidney_figures/05_pilot_gse157079_lsi_umap_from_matrix.R` | sampled matrix-derived UMAP pilot script | medium |
| `scripts/mouse_kidney_figures/06_all_cells_top_peaks_lsi_umap_from_matrix.R` | all-cell top-peak before-filtering UMAP script | high |

## C. Can Archive Later, But Do Not Delete Now

| relative path | reason | priority |
|---|---|---|
| `results/kidney_notebook1_20260526_172547/` | earlier PACS-only or small/medium run; keep until final audit is complete | low |
| `results/kidney_notebook1_20260526_205633/` | earlier baseline/debug run; not needed for teacher story | low |
| `results/kidney_notebook1_20260526_221904/` | medium baseline run; useful only as development history | low |
| `results/mouse_kidney_figures/gse157079_pilot_matrix_lsi_umap/pilot_counts_sparse.rds` | sampled pilot sparse matrix; reproducible but intermediate | low |
| `results/mouse_kidney_figures/gse157079_pilot_matrix_lsi_umap/pilot_lsi_embedding.csv` | intermediate pilot embedding | low |
| `results/mouse_kidney_figures/gse157079_pilot_matrix_lsi_umap/pilot_lsi_umap_embedding.csv` | intermediate pilot UMAP coordinates | low |
| `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/lsi_embedding.csv` | large intermediate all-cell LSI embedding; keep for now, archive later if storage matters | medium |
| `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/counts_sparse_top_peaks.rds` | large sparse matrix object; important for reproducibility but not meeting-facing | medium |
| `figures/mouse_kidney/*.pdf` | useful for publication-quality export; show PNGs for meeting, keep PDFs for later | medium |
| `cleanup_manifest.md` | previous cleanup planning document | low |

## D. Temporary / Debug Files Likely Safe To Remove Later, After Confirmation

| relative path | reason | priority |
|---|---|---|
| `results/20260526_2318_large_baseline/pacs_function_sources.txt` | useful debug snapshot of PACS functions; not needed for meeting | low |
| `results/20260526_2318_large_baseline/run.log` | benchmark runtime log; keep until final archival decision | low |
| `results/kidney_notebook1_*/pacs_function_sources.txt` | repeated debug snapshots from earlier runs | low |
| `results/kidney_notebook1_*/run.log` | earlier runtime logs | low |
| `figures/mouse_kidney/.gitkeep` | harmless placeholder; keep unless directory policy changes | low |
| `results/mouse_kidney_figures/.gitkeep` | harmless placeholder; keep unless directory policy changes | low |

## Suggested Top 8 Files/Figures To Show Tomorrow

1. `results/project_progress_summary_for_teacher.md`
2. `results/20260526_2318_large_baseline/summary.csv`
3. `figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png`
4. `figures/mouse_kidney/pacs_permuted_qq_plot.png`
5. `results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md`
6. `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png`
7. `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png`
8. `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md`

## Notes

- Do not delete anything before the meeting.
- Keep the current results as evidence that the project has progressed from
  Notebook 1 benchmark reproduction to full matrix-level GSE157079 analysis.
- After the meeting, cleanup can focus on old `kidney_notebook1_*` directories
  and repeated debug logs, but only after explicit confirmation.
