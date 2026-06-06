# P56 PACS Batch-Effect Peak Filtering UMAP Report

## Goal
Generate a first P56-specific PACS paper-style before/after UMAP demonstration.

## Input Files
- matrix_file: `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz`
- metadata_csv: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv`
- peak_file: `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_peak_list.csv.gz`

## P56 Cell Selection
- P56 cells before cell type filtering: 13526
- P56 cells after cell type filtering: 13526

### Batch Table
```text
      level    n
 P56_batch1 7129
 P56_batch2 6397
```

### Cell Type Table
```text
   level    n
      PT 4757
     LOH 1947
     PT2 1482
    Endo 1143
 stroma2  879
     DCT  741
     CNT  560
      PC  478
      IC  409
  PT_out  237
 stroma1  208
    Podo  196
  immune  185
  NP_LOH  161
      NP  143
```

## Cell Type Filtering
- min_cells_per_celltype: 30
- min_cells_per_batch_per_celltype: 5
- removed cell types: none

## Matrix Dimensions
```text
%%MatrixMarket matrix coordinate integer general
28316 300755 166121193
```
- P56 sparse matrix before empty filtering: 13526 x 20000; nnz=40539709
- removed empty cells: 0
- removed empty peaks: 0
- P56 sparse matrix after empty filtering: 13526 x 20000; nnz=40539709

## Top Peak Selection
- requested n_top_peaks: 20000
- selected top peaks: 20000
- first pass processed coordinate lines: 166121193
- selected P56 entries in first pass: 70640242

## Before-Filtering UMAP Settings
- TF-IDF dim: 13526 x 20000
- LSI dim: 13526 x 50
- UMAP LSI dims: 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30

## PACS Function Signatures Inspected
See `pacs_function_signatures.txt`.

## PACS Model Used
- PACS input orientation: peaks x cells.
- Original matrix orientation: cells x peaks.
- full model: `~ cell_type + batch`
- null model: `~ cell_type`
- `cap_rates` are relative cell-depth rates from selected top peaks, scaled to max 0.99. This is a practical first-demonstration approximation because GSE157079 does not provide the Notebook 1 author `q_vec`.
- PACS peaks per round: 1000

## PACS Results
- tested peaks: 10000
- FDR cutoff: 0.05
- significant batch peaks: 6208
- retained tested peaks: 3792

## After-Filtering UMAP Settings
- TF-IDF dim: 13526 x 3792
- LSI dim: 13526 x 50
- UMAP LSI dims: 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30

## Output Figures
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top20000_test10000_fdr005_before_by_batch.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top20000_test10000_fdr005_before_by_celltype.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top20000_test10000_fdr005_after_by_batch.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top20000_test10000_fdr005_after_by_celltype.png/pdf`
- combined: combined figure generated with patchwork

## Interpretation
Compare the before/after batch-colored plots to assess whether P56_batch1 and P56_batch2 separation weakens after removing PACS-significant batch peaks. Compare cell-type-colored plots to check whether biological cell type structure remains.

## Limitations
- P56-only analysis.
- Top detected peaks only.
- First author-style reconstruction, not the final full-dataset PACS paper figure.
- Uses relative depth-derived cap_rates because GSE157079 lacks the author-provided q_vec used in Notebook 1.
