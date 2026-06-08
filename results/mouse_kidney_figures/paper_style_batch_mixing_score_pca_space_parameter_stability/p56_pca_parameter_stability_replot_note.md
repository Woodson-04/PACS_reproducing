# P56 PCA-logNorm Parameter Stability Replot Note

## Input

- CSV: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space_parameter_stability/p56_pca_parameter_stability_scores.csv`

## Outputs

- PNG: `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pca_space_parameter_stability_scores_replot.png`
- PDF: `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pca_space_parameter_stability_scores_replot.pdf`
- Report material PNG: `/home/woodson/PACS_reproducing/report/pdf_report/materials/figures/gse157079_p56_pca_space_parameter_stability_scores_replot.png`

## Scope

- This replot contains only PCA-logNorm parameter stability scores.
- It does not include GEO reference data.
- It does not include UMAP-space auxiliary metrics.
- It reads existing PCA score CSV output and does not rerun PACS, UMAP, MatrixMarket streaming, TF-IDF/LSI, or PCA scoring.

## Axis Label Handling

- The x-axis labels use `angle = 35, hjust = 1, vjust = 1`, so labels tilt toward the lower right.

## PC1:30 Summary

```text
   setting_label  stage normalized_batch_mixing_score
2      5000/5000 before                    0.05103657
5      5000/5000  after                    0.36788860
8    10000/10000 before                    0.05103657
11   10000/10000  after                    0.30972340
14   20000/10000 before                    0.03982079
17   20000/10000  after                    0.29280673
```

## Interpretation

All three settings show after > before, supporting that PACS filtering improves batch mixing in PCA-logNorm space across parameter choices.

The 5000/5000 after score is highest, but it retains the fewest peaks and may be more aggressive. The 10000/10000 setting remains the main display setting because it provides a more balanced tradeoff among batch mixing, retained features, and cell-type preservation.
