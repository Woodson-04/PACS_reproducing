# Baseline Official Methods Review

Date: 2026-05-26

This note records how the current clean-room baseline methods relate to the
official or commonly documented implementations behind Notebook 1. The original
author helper file `other_methods_for_differential_updated.R` has not been
recovered, so these conclusions are about alignment, not exact equivalence.

## Seurat

Reference implementation to compare against:

- `Seurat::FindMarkers(test.use = "LR")`
- Seurat reference: https://satijalab.org/seurat/reference/findmarkers

Seurat's LR test fits a logistic regression framework for differential feature
testing, and `latent.vars` can be used to include covariates such as sequencing
depth. In a full Seurat workflow, the input would usually be a Seurat object
with identities and metadata, and the implementation details are handled inside
`FindMarkers`.

Current clean-room implementation:

```text
full: group ~ log_depth + peak_value
null: group ~ log_depth
```

The current implementation uses raw peak counts as `peak_value`, total read
depth as `log_depth`, and a likelihood-ratio test between the nested models.
This is close to the conceptual Seurat LR baseline used in the paper, while
avoiding a heavy dependency and avoiding conversion of each sampled matrix into
a Seurat object.

Recommendation: keep the clean-room LR for now. A direct
`Seurat::FindMarkers(test.use = "LR")` path could be added later if exact Seurat
object construction and the author's preprocessing choices are recovered.

## snapATAC

References to compare against:

- SnapATAC DAR helper documentation/source, commonly exposed as `findDAR`
- edgeR exact test: https://rdrr.io/bioc/edgeR/man/exactTest.html

Notebook 1 explicitly labels this baseline as:

```text
snapATAC-- edgeR method
```

and calls:

```text
snapATAC_method(data_matrix_pos_2, data_matrix_neg_2, bcv = 0.4)
```

Current clean-room implementation:

- Uses binarized peak x cell matrices.
- Builds an `edgeR::DGEList`.
- Applies `edgeR::calcNormFactors`.
- Runs `edgeR::exactTest` with fixed `dispersion = bcv^2`.
- Defaults to `bcv = 0.4`, matching the notebook call.

This is substantially closer to the notebook comment than the earlier logistic
regression surrogate. It also explains the conservative Type I error seen in
the medium test.

Recommendation: keep the edgeR fixed-dispersion implementation unless the
author's original `snapATAC_method` is recovered. Installing and calling the
full SnapATAC package is probably unnecessary for Notebook 1 if the helper was
just a thin edgeR wrapper.

## ArchR

References to compare against:

- `ArchR::getMarkerFeatures`
- ArchR reference: https://www.archrproject.com/reference/getMarkerFeatures.html
- Typical marker settings include `testMethod = "wilcoxon"` and bias variables
  such as `bias = c("TSSEnrichment", "log10(nFrags)")`.

ArchR marker testing is designed around an `ArchRProject` and Arrow files, with
cell-level quality metrics and bias-aware background matching. The current
Notebook 1 reproduction only has sampled peak x cell matrices and per-cell read
depth estimates available in the baseline function call:

```text
archR_method(data_matrix_pos, data_matrix_neg)
```

Current clean-room implementation:

- Computes cell depth from the two input matrices.
- Performs nearest-neighbor depth matching between groups.
- Tests log1p depth-normalized counts with Wilcoxon rank-sum.

This captures part of ArchR's bias-aware idea, especially depth matching, but it
does not reproduce ArchR's full marker machinery. Missing pieces include an
`ArchRProject`, Arrow-backed matrices, TSSEnrichment, background group
selection, and any package-internal filtering/model details used by the author's
helper.

Recommendation: treat ArchR as the main remaining uncertainty. Keep the current
depth-matched Wilcoxon implementation for matrix-only runs, but do not describe
it as official ArchR output.

## Fisher

The Fisher baseline is a standard two-sided 2x2 exact test on binarized
accessibility:

```text
rows:    accessible, inaccessible
columns: positive group, negative group
```

The orientation of the 2x2 table does not change the two-sided Fisher p-value.
No external official package is needed beyond base R `stats::fisher.test`.

Recommendation: current implementation is appropriate for the clean-room
baseline.

## Summary

Closest current alignments:

- Fisher: standard and stable.
- snapATAC: now close to the notebook's edgeR/`bcv = 0.4` description.
- Seurat: conceptually close to Seurat LR with depth as a latent covariate.

Main uncertainty:

- ArchR: still an approximation because the available data are matrix-only and
  do not include the full ArchR project context.

Recommended next run:

- Run the large baseline once, still labeling results as
  `clean_room_reimplemented`.
- If large baseline results diverge from the paper/notebook, inspect ArchR
  first, then Seurat object-level details. The snapATAC edgeR implementation is
  less likely to be the dominant mismatch after the fixed-BCV change.
