# Current Important Results Summary

This file summarizes the current PACS reproduction results that are most worth
keeping, presenting, and using for the next analysis decisions.

## 1. Main Story

The project now has two completed result layers:

1. Notebook 1 benchmark reproduction:
   - PACS main method was reproduced on real kidney benchmark data.
   - PACS Type I error and power are close to the author Notebook 1 values.
   - Baseline methods remain clean-room reimplementations because the original
     author baseline helper file was not found.

2. GSE157079 mouse kidney figure reconstruction:
   - The public GSE157079 matrix, metadata, and peak list were aligned.
   - Matrix-derived TF-IDF/LSI/UMAP was successfully computed without using the
     precomputed GEO UMAP coordinates.
   - P56-specific PACS batch-effect peak filtering was run.
   - Batch mixing improved after PACS filtering in UMAP space, LSI space, and
     PCA-logNorm space.

The most defensible current conclusion is:

> PACS filtering reduces P56 batch-associated chromatin accessibility structure
> while preserving meaningful cell-type structure. This is a successful
> PACS-paper-style proof of concept, not yet an exact full Figure 4 reproduction.

## 2. Notebook 1 Benchmark

Key output:

- `results/20260526_2318_large_baseline/summary.csv`

Run configuration:

- `n_repeat = 5`
- `n_cell_sample = 500`
- `n_features_sample = 10000`
- `run_baselines = TRUE`

Result table:

| method | Type I error | power | note |
|---|---:|---:|---|
| PACS / our | 0.04008 | 0.83337 | main reproduced method |
| seurat | 0.06342 | 0.82344 | clean-room baseline |
| archR | 0.04096 | 0.67437 | clean-room approximation |
| snapATAC | 0.01810 | 0.76094 | clean-room edgeR-style baseline |
| fisher | 0.02208 | 0.76630 | binary Fisher exact test |

Important interpretation:

- PACS Type I error is close to 0.05 and power is close to the author notebook.
- Baseline numbers should be labeled as clean-room reimplemented, not exact
  author baseline reproduction.
- `q.r` includes a local fixed PACS sparse wrapper because PACS 0.2.2
  `pacs_test_sparse()` has a mixed cumulative/logit branch rownames bug.

Recommended files to show:

- `figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png`
- `figures/mouse_kidney/pacs_permuted_qq_plot.png`
- `notebook1_reproduction_report.md`

## 3. GSE157079 Data Alignment

Key output:

- `results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md`
- `results/mouse_kidney_figures/gse157079_metadata_merged.csv`

Confirmed data dimensions:

| object | count |
|---|---:|
| cells | 28316 |
| peaks | 300755 |
| nonzero matrix entries | 166121193 |
| metadata rows | 28316 |
| peak list rows | 300755 |

Matrix orientation:

- cell x peak

Metadata standard columns:

- `row_index`
- `cell_barcode`
- `sample`
- `cell_type`
- `umap_1`
- `umap_2`

Sample counts:

| sample | cells |
|---|---:|
| P56_batch1 | 7129 |
| P56_batch2 | 6397 |
| P0_batch1 | 5993 |
| P0_batch2 | 5436 |
| P21_batch1 | 3361 |

Cell type counts:

| cell type | cells |
|---|---:|
| PT | 7412 |
| LOH | 3628 |
| stroma2 | 3539 |
| Endo | 2368 |
| NP | 2232 |
| PT2 | 1943 |
| stroma1 | 1479 |
| PC | 1325 |
| DCT | 950 |
| Podo | 912 |
| CNT | 879 |
| IC | 599 |
| PT_out | 446 |
| immune | 386 |
| NP_LOH | 218 |

Important interpretation:

- The four GSE157079 files are sufficient for core matrix-derived UMAP work.
- They are not sufficient for exact author reproduction without the author's
  preprocessing and parameter details.

## 4. All-Cell Matrix-Derived Top-Peak UMAP

Key output:

- `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md`

Run summary:

| item | value |
|---|---:|
| cells | 28316 |
| selected top peaks | 20000 |
| retained nonzeros | 85801336 |
| sparse matrix dimension | 28316 x 20000 |
| LSI dimension | 28316 x 50 |
| UMAP input dimensions | LSI 2:30 |

Important result:

- A new UMAP was computed from the large count matrix.
- The precomputed GEO `umap_1` / `umap_2` columns were not used.
- The UMAP shows biological cell-type structure and sample-associated structure.
- This is the current best before-filtering all-cell reference.

Recommended files to show:

- `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png`
- `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png`

## 5. P56 PACS Batch Filtering

Current preferred setting:

- `n_top_peaks = 10000`
- `max_pacs_peaks = 10000`
- `fdr_cutoff = 0.05`
- result directory:
  `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved`

P56 subset:

| item | value |
|---|---:|
| P56 cells | 13526 |
| P56_batch1 cells | 7129 |
| P56_batch2 cells | 6397 |
| cell types | 15 |

PACS model:

- PACS input orientation: peaks x cells
- full model: `~ cell_type + batch`
- null model: `~ cell_type`
- `cap_rates`: relative depth-derived rates, because public GSE157079 files do
  not provide the author Notebook 1 `q_vec`.

PACS peak results:

| item | value |
|---|---:|
| tested peaks | 10000 |
| significant batch peaks | 6305 |
| retained peaks | 3695 |

Recommended files to show:

- `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png`
- `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_before_by_batch.png`
- `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_after_by_batch.png`
- `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_before_by_celltype.png`
- `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_after_by_celltype.png`

## 6. Batch Mixing Quantification

The main quantitative finding is that PACS filtering improves batch mixing in
three coordinate spaces.

### UMAP-Space Paper-Style Score

Key output:

- `results/mouse_kidney_figures/paper_style_batch_mixing_score/paper_style_batch_mixing_scores.csv`

| setting | normalized batch mixing score |
|---|---:|
| P56 10000/10000 before | 0.02649 |
| P56 10000/10000 after | 0.65073 |
| P56 20000/10000 before | 0.01978 |
| P56 20000/10000 after | 0.60879 |
| GEO precomputed UMAP P56 | 0.93903 |

Interpretation:

- UMAP visualization-space mixing improves strongly after PACS filtering.
- GEO precomputed UMAP is more mixed, but it likely reflects a different atlas
  pipeline and should be treated as a public reference rather than the PACS
  paper filtered UMAP.

### LSI-Space Paper-Style Score

Key output:

- `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_10000_lsi_space_paper_style_scores.csv`

| setting | normalized batch mixing score |
|---|---:|
| before, LSI 2:30 | 0.05201 |
| after, LSI 2:30 | 0.33434 |

Dimension sensitivity:

| setting | normalized batch mixing score |
|---|---:|
| before, LSI 1:30 | 0.05759 |
| after, LSI 1:30 | 0.28025 |
| before, LSI 2:30 | 0.05201 |
| after, LSI 2:30 | 0.33434 |
| before, LSI 2:50 | 0.07925 |
| after, LSI 2:50 | 0.37298 |
| before, LSI 1:50 | 0.07658 |
| after, LSI 1:50 | 0.30876 |

Interpretation:

- LSI-space score improves after filtering across multiple dimension choices.
- LSI 2:30 is the most appropriate default because LSI 1 is depth-associated.

### LSI-Depth Correlation

Key output:

- `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_depth_correlation_report.md`
- `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_depth_correlation.csv`

| stage | component | Spearman correlation with depth | Pearson correlation with depth |
|---|---|---:|---:|
| before | LSI_1 | -0.9983 | -0.9739 |
| after | LSI_1 | 0.9976 | 0.9830 |

Interpretation:

- LSI_1 is essentially a depth axis.
- Excluding LSI_1 and using LSI 2:30 is justified.

### PCA-LogNorm Paper-Style Score

Key output:

- `results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space/p56_10000_pca_space_scores.csv`

| PCA dimensions | before | after |
|---|---:|---:|
| PC1:20 | 0.03977 | 0.25760 |
| PC1:30 | 0.05104 | 0.30972 |
| PC1:50 | 0.06851 | 0.42390 |

Interpretation:

- PCA-logNorm is the closest current analogue to the PACS paper normalized PCA
  mixing metric.
- The PCA-space result supports PACS filtering: after-filtering scores are much
  higher than before-filtering scores.

Recommended file to show:

- `figures/mouse_kidney/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png`

## 7. Parameter Comparison

Key output:

- `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_comparison/p56_pacs_filter_setting_comparison.csv`

| setting | significant peaks | retained peaks | batch entropy | same-batch fraction | cell type silhouette | same-celltype fraction |
|---|---:|---:|---:|---:|---:|---:|
| 5000/5000/FDR0.05 | 3140 | 1860 | 0.69564 | 0.69693 | 0.15933 | 0.66504 |
| 10000/10000/FDR0.05 | 6305 | 3695 | 0.71111 | 0.68101 | 0.30198 | 0.77592 |

Interpretation:

- The 10000/10000/FDR0.05 setting is preferred over 5000/5000/FDR0.05.
- It improves batch mixing and better preserves cell-type structure.

## 8. Most Important Files To Keep And Show

### Core reports

- `results/current_important_results_summary.md`
- `notebook1_reproduction_report.md`
- `results/20260526_2318_large_baseline/summary.csv`
- `results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md`
- `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md`
- `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_pacs_batch_filter_umap_report.md`
- `results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space/p56_10000_pca_space_score_report.md`
- `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_paper_style_metric_final_interpretation.md`

### Core figures

- `figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png`
- `figures/mouse_kidney/pacs_permuted_qq_plot.png`
- `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png`
- `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png`
- `figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png`
- `figures/mouse_kidney/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png`
- `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before_after.png`

## 9. Results That Are Useful But Secondary

- GEO precomputed UMAP plots and quantification:
  useful as a public reference, but not the main PACS paper-style result.
- Pilot sampled matrix-derived UMAP:
  useful as proof that matrix streaming works, but superseded by the all-cell
  top-peak UMAP.
- 5000/5000/FDR0.05 P56 run:
  useful for parameter comparison, but superseded by 10000/10000/FDR0.05.
- 20000/10000/FDR0.05 P56 run:
  useful as a stability check, but current main setting remains 10000/10000.

## 10. Current Limitations

- The P56 PACS analysis is a two-batch, same-age subset analysis.
- It is not yet a full all-age, all-feature Figure 4 reproduction.
- Public GSE157079 data do not include the author Notebook 1 `q_vec`, so
  `cap_rates` use a depth-derived approximation.
- Baselines in Notebook 1 are clean-room reimplementations.
- GEO precomputed UMAP is not the same thing as PACS-filtered UMAP.

## 11. Recommended Presentation Order

1. Notebook 1 PACS benchmark table and QQ plot.
2. GSE157079 matrix alignment proof.
3. All-cell matrix-derived top-peak UMAP by sample and cell type.
4. P56 PACS before/after four-panel UMAP.
5. PCA/LSI/UMAP normalized batch mixing score comparison.
6. LSI-depth correlation showing why LSI 1 was excluded.

## 12. Short Verbal Summary

Notebook 1 PACS reproduction is successful for the main method. The project then
moved to GSE157079 mouse kidney data. The public matrix, metadata, and peak list
are aligned, and a matrix-derived all-cell LSI/UMAP was generated from the count
matrix without using GEO's precomputed UMAP. In P56 cells, PACS identified 6305
batch-associated peaks among the top 10000 tested peaks while controlling for
cell type. Removing those peaks reduced batch separation and improved normalized
batch mixing in UMAP, LSI, and PCA-logNorm spaces. This supports the PACS
paper-style claim that PACS can identify and filter batch-associated chromatin
features, while preserving biological cell-type structure.
