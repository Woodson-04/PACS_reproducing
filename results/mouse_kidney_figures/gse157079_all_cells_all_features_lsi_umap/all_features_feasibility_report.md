# GSE157079 All-Feature LSI UMAP Feasibility Report

This dry run checks whether an all-cell, all-feature matrix-derived UMAP is feasible.
No full matrix was streamed, no LSI/UMAP was run, and no dense matrix was created.

## Conceptual Correction

- Previous all-cell UMAP used top 20000 detected peaks and should be called top-peak UMAP.
- This target all-feature UMAP uses all 300755 peaks from the GSE157079 cell-by-peak matrix.

## Arguments

```text
gse_dir = /home/woodson/biostatistic/pacs/GSE157079
metadata_csv = /home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv
out_dir = /home/woodson/PACS_reproducing/results/mouse_kidney_figures
fig_dir = /home/woodson/PACS_reproducing/figures/mouse_kidney
seed = 1
chunk_lines = 100000
progress_every = 5000000
n_lsi = 50
umap_lsi_start = 2
umap_lsi_end = 30
save_counts = FALSE
dry_run = TRUE
matrix_file = /home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz
```

## Matrix Header

```text
%%MatrixMarket matrix coordinate integer general
28316 300755 166121193
```

## Alignment Checks

- metadata rows: 28316
- peak list rows: 300755
- row_index covers matrix cells: TRUE
- peak_index covers matrix peaks: TRUE
- metadata rows match matrix cells: TRUE
- peak rows match matrix peaks: TRUE

## Estimated dgCMatrix Memory

```text
           component      bytes gigabytes
           x numeric 1328969544     1.238
       i row indices  664484772     0.619
   p column pointers    1203024     0.001
  raw dgCMatrix core 1994657340     1.858
 rough with overhead 2991986010     2.787
```

## Feasibility Interpretation

The raw dgCMatrix core is expected to require roughly the memory shown above.
The full run will require additional memory for coordinate chunks, sparse operations, TF-IDF, irlba workspace, and UMAP.
If available RAM is limited, the full all-feature run may fail during matrix construction, TF-IDF, or LSI.

## Conclusion

- dry-run inputs aligned: TRUE
- full all-feature UMAP should only be attempted on a Linux session with sufficient RAM.
- If memory fails, keep the previous top-20000 UMAP wording and state that all-feature UMAP was attempted but was computationally limited.
