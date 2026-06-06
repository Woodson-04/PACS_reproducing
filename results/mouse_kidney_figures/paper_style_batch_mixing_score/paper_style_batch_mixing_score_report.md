# PACS-paper-style normalized batch mixing score

## Goal

This report computes the PACS-paper-style normalized batch mixing score for existing P56 embeddings and the GEO-provided precomputed UMAP reference.

Higher normalized batch mixing score indicates better batch mixing. A value near 1 means the observed different-batch neighbor fraction approaches the expectation from the cell_type-by-batch composition matrix.

## Coordinate-space caveat

No LSI coordinate columns were found in the available P56 embedding CSV files. Therefore, all scores here are UMAP-space approximations, not PCA/LSI-space scores. They should not be directly compared numerically to the PACS paper's PCA-space 0.122 -> 0.358 statement.

## Summary scores

```text
                            space_name coordinate_space_type n_cells  k
1 GEO precomputed UMAP P56, UMAP space            UMAP-space   13526 30
2               P56 10000/10000 before            UMAP-space   13526 30
3                P56 10000/10000 after            UMAP-space   13526 30
4               P56 20000/10000 before            UMAP-space   13526 30
5                P56 20000/10000 after            UMAP-space   13526 30
  n_batches n_cell_types observed_batch_mixing_score
1         2           15                 0.460320864
2         2           15                 0.012987333
3         2           15                 0.318990586
4         2           15                 0.009697373
5         2           15                 0.298432648
  expected_batch_mixing_score normalized_batch_mixing_score
1                    0.490207                    0.93903366
2                    0.490207                    0.02649357
3                    0.490207                    0.65072631
4                    0.490207                    0.01978220
5                    0.490207                    0.60878905
      coordinate_columns_used
1               umap_1;umap_2
2 before_umap_1;before_umap_2
3   after_umap_1;after_umap_2
4 before_umap_1;before_umap_2
5   after_umap_1;after_umap_2
```

## Interpretation

- The GEO row is a public GSE157079 atlas UMAP reference embedding, not PACS-filtered ground truth.
- Compare before vs after rows within the same P56 setting to assess whether PACS filtering improves normalized batch mixing.
- Because the current available P56 CSVs contain UMAP coordinates but no LSI_2:LSI_30 columns, the current comparison is UMAP-space unless future scripts save LSI coordinates.

## Missing inputs

No expected input files were missing.

## Output files

- `results/mouse_kidney_figures/paper_style_batch_mixing_score/paper_style_batch_mixing_scores.csv`
- `figures/mouse_kidney/gse157079_paper_style_batch_mixing_score_comparison.png/pdf`
- Per-cell CSV files were also written for each computed row.
