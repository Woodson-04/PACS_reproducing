# P56 PACS Filter Setting Comparison

## Compared Settings

| setting | run name | tested peaks | significant batch peaks | retained peaks | status |
|---|---|---:|---:|---:|---|
| 5000/5000/FDR0.05 | default | 5000 | 3140 | 1860 | completed |
| 10000/10000/FDR0.05 | p56_top10000_test10000_fdr005 | 10000 | 6305 | 3695 | current main |

## Quantitative Metrics

| setting | after batch silhouette | after normalized batch entropy | after same-batch fraction | after batch prediction accuracy | after cell type silhouette | after same-celltype fraction |
|---|---:|---:|---:|---:|---:|---:|
| 5000/5000/FDR0.05 | 0.0266 | 0.6956 | 0.6969 | 0.5899 | 0.1593 | 0.6650 |
| 10000/10000/FDR0.05 | 0.0133 | 0.7111 | 0.6810 | 0.5754 | 0.3020 | 0.7759 |

## Interpretation

The 10000/10000/FDR0.05 setting is currently preferred. It improves batch mixing relative to the 5000 setting and better preserves cell-type structure.

The next planned stability check is:

```text
n_top_peaks = 20000
max_pacs_peaks = 10000
fdr_cutoff = 0.05
run_name = p56_top20000_test10000_fdr005
```

After that run completes, this comparison should be updated to include the third setting and decide whether 10000/10000 remains the main result or whether 20000/10000 should replace it.
