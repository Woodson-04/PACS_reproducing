# GSE157079 Pilot Matrix-Derived LSI UMAP Report

This pilot starts from the GSE157079 MatrixMarket cell-by-peak count matrix.
The precomputed GEO `umap_1`/`umap_2` columns were ignored for embedding construction.

## Command Arguments

```text
gse_dir = /home/woodson/biostatistic/pacs/GSE157079
metadata_csv = /home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv
out_dir = /home/woodson/PACS_reproducing/results/mouse_kidney_figures
fig_dir = /home/woodson/PACS_reproducing/figures/mouse_kidney
n_cells_per_sample = 500
n_peaks = 20000
seed = 1
chunk_lines = 1000000
progress_every = 5000000
```

## Matrix Header

```text
%%MatrixMarket matrix coordinate integer general
28316 300755 166121193
```

## Sampling

- Selected cells before empty filtering: 2500
- Selected peaks before empty filtering: 20000

## Sparse Matrix

- Dimensions before empty filtering: 2500 x 20000
- Nonzeros before empty filtering: 964463
- Removed empty cells: 0
- Removed empty peaks: 107
- Dimensions after empty filtering: 2500 x 19893
- Nonzeros after empty filtering: 964463

## Sample Table

```text
      level   n
  P0_batch1 500
  P0_batch2 500
 P21_batch1 500
 P56_batch1 500
 P56_batch2 500
```

## Cell Type Table

```text
   level   n
      PT 692
     LOH 313
 stroma2 299
    Endo 218
      NP 189
     PT2 153
 stroma1 136
      PC 132
    Podo  84
     CNT  73
     DCT  64
      IC  52
  PT_out  39
  immune  34
  NP_LOH  22
```

## Peak Detection Summary

```text
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1.00    7.00   15.00   48.48   39.00 1596.00 
```

## TF-IDF / LSI

- TF-IDF dimensions: 2500 x 19893
- LSI dimensions: 2500 x 30
- LSI dimensions used for UMAP: 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30
- Spearman correlations between LSI components and pilot depth:

```text
        LSI_1         LSI_2         LSI_3         LSI_4         LSI_5 
 9.618481e-01 -4.113931e-02 -2.609415e-01 -1.946162e-01  5.267373e-02 
        LSI_6         LSI_7         LSI_8         LSI_9        LSI_10 
-4.751128e-03  2.114147e-01 -1.877937e-01 -2.353118e-01 -1.592083e-01 
       LSI_11        LSI_12        LSI_13        LSI_14        LSI_15 
 1.377251e-01 -1.383246e-02 -4.316810e-02  2.669763e-02  8.228378e-03 
       LSI_16        LSI_17        LSI_18        LSI_19        LSI_20 
-2.544350e-02  1.062639e-02  7.964175e-03  2.651415e-02 -8.104287e-03 
       LSI_21        LSI_22        LSI_23        LSI_24        LSI_25 
-2.083025e-02 -1.182966e-02  1.758890e-02 -2.850266e-03 -3.333067e-02 
       LSI_26        LSI_27        LSI_28        LSI_29        LSI_30 
-7.239826e-03  3.729362e-02 -4.344703e-02 -4.506281e-02 -8.344182e-05 
```

## UMAP

- Package: uwot
- n_neighbors: 30
- min_dist: 0.3
- metric: cosine
- UMAP rows: 2500

## Output Files

- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_pilot_matrix_lsi_umap/pilot_counts_sparse.rds`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_pilot_matrix_lsi_umap/pilot_metadata.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_pilot_matrix_lsi_umap/pilot_peak_indices.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_pilot_matrix_lsi_umap/pilot_lsi_embedding.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_pilot_matrix_lsi_umap/pilot_lsi_umap_embedding.csv`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_pilot_matrix_lsi_umap_by_sample.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_pilot_matrix_lsi_umap_by_sample.pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_pilot_matrix_lsi_umap_by_celltype.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_pilot_matrix_lsi_umap_by_celltype.pdf`

## Conclusion

Pilot matrix-derived UMAP succeeded. This is not the final PACS paper figure;
it only verifies that the large cell-by-peak matrix can be streamed into a
sparse pilot subset and used for TF-IDF/LSI/UMAP without relying on the
precomputed GEO UMAP coordinates.
