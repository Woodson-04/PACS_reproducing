# P56 paper-style metric final interpretation

## Metric hierarchy

1. The PACS paper uses normalized PCA mixing.
2. Our closest author-like metric is PCA-logNorm PC1:30 or PC1:50.
3. LSI_2:30 is a scATAC-aware analogue because TF-IDF/LSI is standard for sparse chromatin accessibility matrices.
4. UMAP-space score is visualization-space only.

## Main quantitative conclusion

- PCA-logNorm score improves strongly after PACS filtering.
- LSI-space score improves strongly after PACS filtering.
- UMAP-space score also improves, but should not be used as the main paper-style metric.
- Main LSI_2:30 normalized score: before = 0.05201, after = 0.3343.

## LSI_1 depth correlation

- before strongest Spearman: LSI_1 = -0.9983
- after strongest Spearman: LSI_1 = 0.9976
- LSI_1 correlations: before Spearman=-0.9983, Pearson=-0.9739; after Spearman=0.9976, Pearson=0.983

## Why not directly equate to the paper values

Do not directly equate our scores with the paper's 0.122 -> 0.358 because this analysis uses a P56-only two-batch subset, a top10000 feature universe, depth-derived cap_rates, different normalization, and may differ from the author-specific adult kidney feature set.

## Final current recommendation

Use PCA-logNorm PC1:30 or PC1:50 as the most paper-like quantitative metric, LSI_2:30 as the scATAC-aware supporting metric, and UMAP-space scores only as visualization-level supporting evidence.
