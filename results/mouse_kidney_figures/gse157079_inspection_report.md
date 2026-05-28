# GSE157079 Inspection Report

Source directory: `/home/woodson/biostatistic/pacs/GSE157079`

This report was generated without extracting or modifying source `.gz` files.

## umap

- File: `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_UMAP_coordinates.csv.gz`
- Exists: TRUE
- Size bytes: 540016
- Preview lines read: 6
- Delimited column count from first line: 3
- First-line fields: ``, `umap-1`, `umap-2`

Preview:

```text
,umap-1,umap-2
1,-9.08304511167941,0.944673617394453
2,-7.99218427758161,0.814818378724315
3,-0.899466163192,4.80897364847544
4,-8.36913369759072,1.03759170449244
5,-10.6504242798244,7.62190252398744
```

- Parsed column names: `V1`, `umap-1`, `umap-2`
- Parsed dimensions: 28316 rows x 3 columns
- Header-like first data row: not detected
- Clean dimensions: 28316 rows x 3 columns
- Candidate barcode columns: none
- Candidate UMAP columns: umap-1, umap-2
- Candidate cell type/annotation columns: none
- Candidate peak coordinate columns: none

## metadata

- File: `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_metadata.csv.gz`
- Exists: TRUE
- Size bytes: 221285
- Preview lines read: 6
- Delimited column count from first line: 4
- First-line fields: ``, `barcodes`, `samples`, `clusters`

Preview:

```text
,barcodes,samples,clusters
1,AAACGAAAGATGTTGA-1,P56_batch2,PT2
2,AAACGAAAGGCTAAAT-1,P56_batch2,PT2
3,AAACGAAAGGGTCCCT-1,P56_batch2,PT
4,AAACGAAAGTCGAGCA-1,P56_batch2,PT2
5,AAACGAAAGTGTCACT-1,P56_batch2,immune
```

- Parsed column names: `V1`, `barcodes`, `samples`, `clusters`
- Parsed dimensions: 28316 rows x 4 columns
- Header-like first data row: not detected
- Clean dimensions: 28316 rows x 4 columns
- Candidate barcode columns: barcodes
- Candidate UMAP columns: none
- Candidate cell type/annotation columns: clusters
- Candidate peak coordinate columns: none

## peak_list

- File: `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_peak_list.csv.gz`
- Exists: TRUE
- Size bytes: 5030623
- Preview lines read: 6
- Delimited column count from first line: 7
- First-line fields: ``, `seqnames`, `start`, `end`, `width`, `strand`, `name`

Preview:

```text
,seqnames,start,end,width,strand,name
1,chr1,3115016,3115356,341,*,chr1:3115016-3115356
2,chr1,3119710,3120787,1078,*,chr1:3119710-3120787
3,chr1,3121548,3121773,226,*,chr1:3121548-3121773
4,chr1,3126118,3126601,484,*,chr1:3126118-3126601
5,chr1,3285439,3285639,201,*,chr1:3285439-3285639
```

- Parsed column names: `V1`, `seqnames`, `start`, `end`, `width`, `strand`, `name`
- Parsed dimensions: 300755 rows x 7 columns
- Header-like first data row: not detected
- Clean dimensions: 300755 rows x 7 columns
- Candidate barcode columns: none
- Candidate UMAP columns: none
- Candidate cell type/annotation columns: none
- Candidate peak coordinate columns: seqnames, start, end, name

## cell_by_peak_matrix

- File: `/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz`
- Exists: TRUE
- Size bytes: 469810081
- Preview lines read: 6
- Delimited column count from first line: 1
- First-line fields: `%%MatrixMarket matrix coordinate integer general`

Preview:

```text
%%MatrixMarket matrix coordinate integer general
28316 300755 166121193
7280 1 1
9528 1 1
10133 1 1
10591 1 1
```

- Matrix format guess: MatrixMarket-like
- Full matrix read: skipped intentionally because this file is large.

