# P56 PACS Batch Filtering Current Main Result Summary

## Why P56-only

The current PACS paper-style reconstruction focuses on P56 cells because P56 has two batches, `P56_batch1` and `P56_batch2`, with enough cells in both batches to support a direct batch-effect test. This makes it a clean setting for asking whether PACS can identify batch-associated accessibility peaks while controlling for biological cell-type composition.

## Batch Contrast

The batch contrast is:

- `P56_batch1`
- `P56_batch2`

The goal is not to remove age-associated biology across P0, P21, and P56. Instead, this P56-only analysis targets technical or batch-associated peak effects within the same age group.

## Why Control Cell Type

PACS was run with cell type as a covariate:

- full model: `~ cell_type + batch`
- null model: `~ cell_type`

This tests whether each peak has a batch-associated effect after accounting for cell-type composition. This is important because P56 cell types differ in chromatin accessibility, and an unadjusted batch test could mistake biological cell-type structure for a batch effect.

## Settings Compared

| setting | n_top_peaks | max_pacs_peaks | FDR cutoff | tested peaks | significant batch peaks | retained tested peaks |
|---|---:|---:|---:|---:|---:|---:|
| 5000/5000/FDR0.05 | 5000 | 5000 | 0.05 | 5000 | 3140 | 1860 |
| 10000/10000/FDR0.05 | 10000 | 10000 | 0.05 | 10000 | 6305 | 3695 |

## Quantitative Comparison

| setting | after batch silhouette | after normalized batch entropy | after same-batch fraction | after batch prediction accuracy | after cell type silhouette | after same-celltype fraction |
|---|---:|---:|---:|---:|---:|---:|
| 5000/5000/FDR0.05 | 0.0266 | 0.6956 | 0.6969 | 0.5899 | 0.1593 | 0.6650 |
| 10000/10000/FDR0.05 | 0.0133 | 0.7111 | 0.6810 | 0.5754 | 0.3020 | 0.7759 |

Lower batch silhouette, lower same-batch fraction, lower batch prediction accuracy, and higher normalized batch entropy all indicate better batch mixing. Higher cell type silhouette and higher same-celltype fraction indicate better preservation of biological cell-type structure.

## Current Preferred Setting

The current preferred P56 result is:

```text
n_top_peaks = 10000
max_pacs_peaks = 10000
fdr_cutoff = 0.05
run_name = p56_top10000_test10000_fdr005
```

This setting improves batch mixing compared with the 5000 setting and also preserves cell-type structure better. The 5000 setting reduces batch separation, but its after-filtering same-celltype fraction and cell type silhouette are lower, suggesting more loss of local biological structure.

## Interpretation

PACS filtering strongly reduces batch-associated UMAP structure in the P56-only analysis. In the 10000 setting, the post-filtering UMAP has weaker batch separability and clearer cell-type organization than the 5000 setting.

This is still a P56-only, top-peak-based reconstruction. It is not yet the final full-dataset PACS paper-level correction. The analysis also uses depth-derived `cap_rates` because the public GSE157079 files do not provide the author's original `q_vec`.

## Recommended Next Step

Run one stability-check setting:

```text
n_top_peaks = 20000
max_pacs_peaks = 10000
fdr_cutoff = 0.05
run_name = p56_top20000_test10000_fdr005
```

This tests whether the 10000 PACS-tested peak result is stable when the top-peak universe is expanded to 20000. Do not jump directly to `max_pacs_peaks = 20000` or exploratory FDR thresholds such as 0.10/0.20 until this stability check is evaluated quantitatively.

Do not judge only by UMAP appearance. Use batch silhouette, normalized batch entropy, same-batch fraction, batch prediction accuracy, cell type silhouette, and same-celltype fraction together.
