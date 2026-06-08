# PCA Parameter Stability Input Audit

This audit checks whether existing P56 PACS result directories contain the files
needed to compute PCA-logNorm paper-style normalized batch mixing scores.

## Candidate Directories

Searched under:

```text
results/mouse_kidney_figures/
```

Candidate directories matching `gse157079_p56_pacs_batch_filter_umap*`:

1. `gse157079_p56_pacs_batch_filter_umap`
2. `gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005`
3. `gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved`
4. `gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005`

## Setting A: 5000/5000/FDR0.05

| item | value |
|---|---|
| result_dir | `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap` |
| n_top_peaks | 10000 |
| tested_peaks | 5000 |
| fdr_cutoff | 0.05 |
| significant_batch_peaks | 3140 |
| retained_peaks | 1860 |
| `p56_counts_top_peaks_sparse.rds` | present |
| `p56_retained_peak_indices.csv` | present |
| `p56_metadata.csv` | present |
| retained column `top_peak_col` | present |
| can compute PCA-logNorm directly | yes |

Note: this default directory used a top-peak matrix of 10000 peaks but tested
5000 peaks with PACS. The setting label is therefore `5000/5000/FDR0.05` for
PACS comparison, while the before matrix contains 10000 top peaks.

## Setting B: 10000/10000/FDR0.05

Preferred directory:

```text
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved
```

| item | value |
|---|---|
| n_top_peaks | 10000 |
| tested_peaks | 10000 |
| fdr_cutoff | 0.05 |
| significant_batch_peaks | 6305 |
| retained_peaks | 3695 |
| `p56_counts_top_peaks_sparse.rds` | present |
| `p56_retained_peak_indices.csv` | present |
| `p56_metadata.csv` | present |
| retained column `top_peak_col` | present |
| can compute PCA-logNorm directly | yes |

The non-`lsi_saved` directory for the same setting also exists, but the
`lsi_saved` directory is preferred because it is the current main setting.

## Setting C: 20000/10000/FDR0.05

| item | value |
|---|---|
| result_dir | `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005` |
| n_top_peaks | 20000 |
| tested_peaks | 10000 |
| fdr_cutoff | 0.05 |
| significant_batch_peaks | 6208 |
| retained_peaks | 3792 |
| `p56_counts_top_peaks_sparse.rds` | present |
| `p56_retained_peak_indices.csv` | present |
| `p56_metadata.csv` | present |
| retained column `top_peak_col` | present |
| can compute PCA-logNorm directly | yes |

## Conclusion

All three settings have the required files to compute PCA-logNorm parameter
stability directly. No 07-script rerun is required at this stage.

Required next step:

```bash
cd /home/woodson/PACS_reproducing
Rscript scripts/mouse_kidney_figures/09g_pca_space_parameter_stability_scores.R
```
