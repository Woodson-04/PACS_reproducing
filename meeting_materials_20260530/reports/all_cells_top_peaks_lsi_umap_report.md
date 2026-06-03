# GSE157079 All-Cell Top-Peak LSI UMAP Report

This run computes a new matrix-derived TF-IDF/LSI/UMAP embedding from all cells and top detected peaks.
The precomputed GEO `umap_1`/`umap_2` columns were not used.
PACS and batch-feature removal were not run in this step.

## Command Arguments

```text
gse_dir = /home/woodson/biostatistic/pacs/GSE157079
metadata_csv = /home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv
out_dir = /home/woodson/PACS_reproducing/results/mouse_kidney_figures
fig_dir = /home/woodson/PACS_reproducing/figures/mouse_kidney
n_top_peaks = 20000
seed = 1
chunk_lines = 100000
progress_every = 5000000
matrix_file = /home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz
```

## Matrix File

`/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz`

## MatrixMarket Header

```text
%%MatrixMarket matrix coordinate integer general
28316 300755 166121193
```

## Metadata

- Rows: 28316
- Filtered rows used in UMAP: 28316

### Sample Table

```text
      level    n
 P56_batch1 7129
 P56_batch2 6397
  P0_batch1 5993
  P0_batch2 5436
 P21_batch1 3361
```

### Cell Type Table

```text
   level    n
      PT 7412
     LOH 3628
 stroma2 3539
    Endo 2368
      NP 2232
     PT2 1943
 stroma1 1479
      PC 1325
     DCT  950
    Podo  912
     CNT  879
      IC  599
  PT_out  446
  immune  386
  NP_LOH  218
```

## Peak List

- Rows: 300755

## First Pass Summary

- Coordinate lines processed: 166121193
- Peaks with nonzero detection: 300751
- Cell depth summary:

```text
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    236    2954    5105    5867    8054   21085 
```

- Peak detection summary:

```text
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    0.0    79.0   171.0   552.3   449.0 23896.0 
```

## Selected Peak Summary

- n_top_peaks requested: 20000
- n_top_peaks selected: 20000
- Selected peak detection summary:

```text
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1757    2327    3405    4290    5473   23896 
```

## Sparse Matrix

- Second pass coordinate lines processed: 166121193
- Second pass retained nonzeros: 85801336
- Dimensions before empty filtering: 28316 x 20000
- Nonzeros before empty filtering: 85801336
- Removed empty cells: 0
- Removed empty peaks: 0
- Dimensions after empty filtering: 28316 x 20000
- Nonzeros after empty filtering: 85801336

## TF-IDF / LSI

- TF-IDF dimensions: 28316 x 20000
- LSI dimensions: 28316 x 50
- LSI dimensions used for UMAP: 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30
- Spearman correlations between LSI components and top-peak depth:

```text
        LSI_1         LSI_2         LSI_3         LSI_4         LSI_5 
 0.9939147114 -0.1373269834 -0.2687361425  0.2604096696 -0.2978351819 
        LSI_6         LSI_7         LSI_8         LSI_9        LSI_10 
 0.0022207495 -0.5153149537  0.1663186744  0.1147177915 -0.1325766825 
       LSI_11        LSI_12        LSI_13        LSI_14        LSI_15 
 0.0345293857 -0.0475856967  0.1882645791  0.0993782911  0.1312990165 
       LSI_16        LSI_17        LSI_18        LSI_19        LSI_20 
 0.0881437540 -0.0055819501 -0.0615550186  0.0186943032  0.0160021163 
       LSI_21        LSI_22        LSI_23        LSI_24        LSI_25 
 0.0102223144  0.0234940816 -0.0964565846  0.0615987753  0.0392525082 
       LSI_26        LSI_27        LSI_28        LSI_29        LSI_30 
-0.0545891545 -0.0123815806  0.0048071468 -0.0469927837  0.0132110634 
       LSI_31        LSI_32        LSI_33        LSI_34        LSI_35 
 0.0182583833 -0.0737962105  0.0140603171 -0.0312941170  0.0487016142 
       LSI_36        LSI_37        LSI_38        LSI_39        LSI_40 
-0.0452514372 -0.0123246029 -0.0309132553 -0.0324291482 -0.0581517196 
       LSI_41        LSI_42        LSI_43        LSI_44        LSI_45 
 0.0016123289 -0.1040595652 -0.0367602959 -0.0128252858  0.0068780048 
       LSI_46        LSI_47        LSI_48        LSI_49        LSI_50 
-0.0930679778 -0.0234963925  0.0640397790 -0.0754359727 -0.0003636665 
```

## UMAP Parameters

- Package: uwot
- n_neighbors: 30
- min_dist: 0.3
- metric: cosine
- UMAP rows: 28316

## Output Files

- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/counts_sparse_top_peaks.rds`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/metadata_with_depth.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/top_peak_indices.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/lsi_embedding.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/lsi_umap_embedding.csv`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.pdf`

## Conclusion

All-cell top-peak matrix-derived UMAP completed successfully if this report was written.
This is the before-PACS-filtering reference candidate, not the final PACS-filtered UMAP.
