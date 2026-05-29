# GSE157079 Matrix Alignment Smoke Test

This report checks file alignment for the PACS paper-style reconstruction route.
No UMAP, PACS, or dense matrix materialization was performed.

## Inputs

- Matrix: `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz`
- Merged metadata: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv`
- Peak list: `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_peak_list.csv.gz`

## MatrixMarket Header

```text
%%MatrixMarket matrix coordinate integer general
28316 300755 166121193
```

- Header is MatrixMarket-like: `TRUE`
- Parsed dimensions: cells=28316, peaks=300755, nonzero=166121193
- Expected dimensions: cells=28316, peaks=300755, nonzero=166121193
- Dimensions match expected: `TRUE`

## Metadata

- Rows: 28316
- Columns: row_index, cell_barcode, sample, cell_type, umap_1, umap_2
- Required columns present: `TRUE`
- Row count matches matrix cells: `TRUE`
- `row_index` covers 1:28316: `TRUE`

First metadata rows:

```text
 row_index       cell_barcode     sample cell_type      umap_1    umap_2
         1 AAACGAAAGATGTTGA-1 P56_batch2       PT2  -9.0830451 0.9446736
         2 AAACGAAAGGCTAAAT-1 P56_batch2       PT2  -7.9921843 0.8148184
         3 AAACGAAAGGGTCCCT-1 P56_batch2        PT  -0.8994662 4.8089736
         4 AAACGAAAGTCGAGCA-1 P56_batch2       PT2  -8.3691337 1.0375917
         5 AAACGAAAGTGTCACT-1 P56_batch2    immune -10.6504243 7.6219025
```

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
- Columns: peak_index, seqnames, start, end, width, strand, name
- Required columns present: `TRUE`
- Row count matches matrix peaks: `TRUE`
- Peak index check: `TRUE`
- Peak index note: peak_index is non-NA and covers 1:300755.

First peak rows:

```text
 peak_index seqnames   start     end width strand                 name
          1     chr1 3115016 3115356   341      * chr1:3115016-3115356
          2     chr1 3119710 3120787  1078      * chr1:3119710-3120787
          3     chr1 3121548 3121773   226      * chr1:3121548-3121773
          4     chr1 3126118 3126601   484      * chr1:3126118-3126601
          5     chr1 3285439 3285639   201      * chr1:3285439-3285639
```

## Coordinate-Line Sanity Check

- Coordinate lines parsed: 10000
- Row index min/max: 14 / 28315
- Column index min/max: 1 / 52
- Coordinate bounds within expected dimensions: `TRUE`
- Value summary:

```text
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      1       1       1       1       1       1 
```

Coordinate rows range from 14 to 28315 and coordinate columns range from 1 to 52 in the first 10000 entries. Given the header dimensions 28316 x 300755, this is consistent with a cell x peak sparse matrix.

## Conclusion

- Smoke test passed: `TRUE`
- The GSE157079 matrix, merged metadata, and peak list appear sufficient and aligned for the next PACS paper-style pipeline.
- The matrix orientation appears to be cell x peak, matching the file name and header dimensions.

Recommended next step: write a small sparse-loading prototype that creates a sparse Matrix object from the MatrixMarket file in a controlled output directory, verifies row/column names from metadata and peak list, and then plans all-feature versus PACS-filtered UMAP without dense materialization.
