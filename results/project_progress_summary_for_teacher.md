# PACS Reproduction Project Progress Summary

Prepared for teacher meeting. This summary covers the PACS Notebook 1 benchmark
stage and the current mouse kidney / GSE157079 figure reconstruction stage.

## 1. Project Objective

The project aims to reproduce key analyses from the PACS paper, especially:

- Notebook 1 benchmark: real kidney data Type I error and power workflow.
- Mouse kidney figure reconstruction: all-feature / batch-effect-filtered UMAP
  logic from the PACS paper.

The current route is **PACS paper-style reconstruction**. The project is not
trying to fully rebuild the original `dev-kidney-snATAC` atlas; that repository
is used only as background/reference.

## 2. Environment And Data Paths

- Current editable project:
  `/home/woodson/PACS_reproducing`
- Old/reference/raw data:
  `/home/woodson/biostatistic/pacs`
- `/home/woodson/biostatistic` should be treated as read-only.

GSE157079 source files:

- `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_UMAP_coordinates.csv.gz`
- `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_metadata.csv.gz`
- `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_peak_list.csv.gz`
- `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz`

No `external_data/GSE157079` replacement matrix was found in the current
project tree during this summary pass. The latest all-cell top-peak run used
the original matrix path above and successfully streamed all 166,121,193
coordinate entries. If there was a separate `gzip -t` warning in a manual
terminal session, keep that terminal record, but the current successful run
shows the matrix was usable for streaming analysis.

## 3. Notebook 1 Benchmark Reproduction

Main script:

- `q.r`

Large benchmark run:

- `n_repeat = 5`
- `n_cell_sample = 500`
- `n_features_sample = 10000`
- `run_baselines = TRUE`
- output: `results/20260526_2318_large_baseline`

Results:

| method | Type I error | power | note |
|---|---:|---:|---|
| our / PACS | 0.04008 | 0.83337 | close to author notebook |
| seurat | 0.06342 | 0.82344 | clean-room baseline |
| archR | 0.04096 | 0.67437 | clean-room approximation |
| snapATAC | 0.01810 | 0.76094 | clean-room edgeR-style baseline |
| fisher | 0.02208 | 0.76630 | binary Fisher exact test |

Interpretation:

- PACS main method is close to the author Notebook 1 result.
- Baseline methods are **clean-room reimplementations** because the original
  author helper file `other_methods_for_differential_updated.R` was not found.
- Therefore, it is appropriate to claim PACS workflow reproduction, but not
  exact author baseline reproduction.
- PACS 0.2.2 has a `pacs_test_sparse()` mixed cumulative/logit branch
  `rownames` bug. `q.r` includes a local fixed wrapper that directly calls
  `pacs_test_cumu()` and `pacs_test_logit()` and safely merges p values.

## 4. Existing Overview Figure Outputs

Current mouse kidney figure outputs under `figures/mouse_kidney/` include:

- `cell_type_counts_barplot.png`
- `pt_loh_depth_distribution.png`
- `pacs_benchmark_t1e_power_barplot.png`
- `pacs_permuted_qq_plot.png`
- `gse157079_umap_by_sample.png/pdf`
- `gse157079_umap_by_celltype.png/pdf`
- `gse157079_pilot_matrix_lsi_umap_by_sample.png/pdf`
- `gse157079_pilot_matrix_lsi_umap_by_celltype.png/pdf`
- `gse157079_all_cells_top_peaks_lsi_umap_by_sample.png/pdf`
- `gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png/pdf`

The strongest teacher-facing figures are the Notebook 1 benchmark plot, QQ
plot, GEO UMAP overview, and the all-cell top-peak matrix-derived UMAP pair.

## 5. GSE157079 Metadata / UMAP Intake

Metadata and precomputed GEO UMAP coordinates were standardized and merged.

Main output:

- `results/mouse_kidney_figures/gse157079_metadata_merged.csv`

Standard columns:

- `row_index`
- `cell_barcode`
- `sample`
- `cell_type`
- `umap_1`
- `umap_2`

Important fixes:

- A `data.table` indexing issue was fixed by using data-frame-safe column
  access.
- The GEO CSV files have an unnamed first column; scripts now explicitly handle
  that as `row_index`.
- 10x-style `cell_barcode` values are not globally unique across samples.
  Therefore, `row_index` is the merge key for metadata and UMAP coordinates,
  while `sample + cell_barcode` is used only as a sanity check.

## 6. GSE157079 Matrix Alignment And Integrity

Smoke-test report:

- `results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md`

Confirmed:

- Matrix orientation: **cell x peak**
- Matrix dimensions: `28316 cells x 300755 peaks`
- Nonzero entries: `166121193`
- Metadata rows: `28316`
- Peak list rows: `300755`
- Peak list has genomic columns including `seqnames`, `start`, `end`, `name`.

The smoke test did not materialize a dense matrix. It checked the MatrixMarket
header and a prefix of sparse coordinate lines, confirming coordinate indices
were within expected cell/peak bounds.

## 7. Author-Code Route Investigation

The original `dev-kidney-snATAC` repository was inspected as reference.

Relevant script:

- `R/combining all batches__create object.R`

It includes the original atlas-style SnapATAC workflow:

- combine five batches;
- create 5 kb bin matrix;
- remove blacklist / mitochondrial / random regions;
- run diffusion maps;
- KNN and clustering;
- UMAP;
- Harmony batch correction;
- after-correction UMAP.

This is useful background, but it is **not** the PACS paper-style script for:

- all-feature UMAP;
- PACS detection of batch-effect features;
- removing batch-effect features;
- reconstructing filtered-feature UMAP.

Current chosen route:

- Use PACS paper-style reconstruction as the main route.
- Treat SnapATAC atlas-style code only as reference.

## 8. Matrix-Derived UMAP Progress

Pilot script:

- `scripts/mouse_kidney_figures/05_pilot_gse157079_lsi_umap_from_matrix.R`

The pilot sampled cells and peaks, streamed the MatrixMarket file, built a
sparse matrix, and computed TF-IDF / LSI / UMAP from the matrix. It did **not**
use the precomputed GEO `umap_1` / `umap_2` columns.

Pilot result:

- 2500 cells
- 19893 nonempty peaks after filtering
- 964,463 retained nonzeros
- LSI/UMAP completed.

All-cell top-peak script:

- `scripts/mouse_kidney_figures/06_all_cells_top_peaks_lsi_umap_from_matrix.R`

Latest all-cell top-peak run:

- all 28,316 cells;
- top 20,000 detected peaks;
- first pass processed all 166,121,193 coordinate lines;
- second pass retained 85,801,336 nonzeros;
- sparse matrix: `28316 x 20000`;
- LSI dimensions: `28316 x 50`;
- UMAP used LSI dimensions 2:30.

Interpretation:

- Sample-colored UMAP shows strong sample-associated structure.
- Cell-type-colored UMAP preserves meaningful biological cell type structure.
- This is a good **before-PACS-filtering baseline** for the next stage.

## 9. Current Scientific Conclusion

The project has successfully moved from Notebook 1 benchmark reproduction to
matrix-level GSE157079 reconstruction. The current all-cell top-peak LSI/UMAP
demonstrates both biological cell type structure and strong sample-associated
structure. This provides a reasonable baseline for the next PACS paper-style
step: detecting batch-effect features with PACS, removing significant
batch-associated peaks, and reconstructing UMAP after filtering.

## 10. Recommended Next Steps

Before the meeting:

- Do not spend more time tuning overview plots.
- Prepare a short teacher-facing story:
  Notebook 1 passed, GSE157079 matrix alignment passed, matrix-derived UMAP
  pipeline works, and the before-filtering baseline is ready.
- Show only the strongest figures.

After the meeting:

1. Run PACS batch-effect peak detection on GSE157079.
2. Remove significant batch-associated peaks.
3. Recompute filtered-feature LSI/UMAP.
4. Create a before/after four-panel figure:
   - all/top features by sample;
   - filtered features by sample;
   - all/top features by cell type;
   - filtered features by cell type.
