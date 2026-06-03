# Notebook 1 Reproduction Report

## Goal

Notebook 1 tests PACS sensitivity and specificity on real kidney data. The
workflow uses PT-only cells for permuted-label Type I error and PT + LOH cells
for actual-label power, with PACS input matrices oriented as peaks x cells.

## Large Run Configuration

```text
n_repeat          = 5
n_cell_sample     = 500
n_features_sample = 10000
run_baselines     = TRUE
output_dir        = results/20260526_2318_large_baseline
baseline_source   = clean_room_reimplemented
```

## Current Result

```text
method    t1e      t1e_sd    power    power_sd
our       0.04008  0.00187   0.83337  0.02264
seurat    0.06342  0.00172   0.82344  0.02897
archR     0.04096  0.00326   0.67437  0.01474
snapATAC  0.01810  0.00089   0.76094  0.01204
fisher    0.02208  0.00081   0.76630  0.01257
```

## Author Notebook 1 Result

Approximate author notebook values:

```text
method    t1e      power
our       0.04514  0.85319
seurat    0.06506  0.84191
archR     0.04438  0.72194
snapATAC  0.00150  0.42740
fisher    0.02640  0.73908
```

## Interpretation

PACS is reproduced successfully for the purpose of validating the workflow.
The PACS Type I error and power are close to the author notebook result. Seurat
and Fisher are also reasonably close.

ArchR and snapATAC differ more because the original author baseline helper
`other_methods_for_differential_updated.R` could not be recovered. This project
therefore uses clean-room reimplemented baselines. These baselines are useful
for exercising the Notebook 1 benchmark workflow, but they should not be
reported as exact author baseline reproductions.

Because the main PACS method is close and the baseline differences have a clear
source, further tuning of Notebook 1 baselines is no longer the main objective.

## PACS 0.2.2 Wrapper Note

The installed PACS 0.2.2 `pacs_test_sparse()` has a mixed-branch rowname bug
when both cumulative-logit and logit peaks are present. The failure occurs in
the convergence merging code after PACS has computed p-values, where row names
can be assigned with a length that does not match the convergence matrix.

`q.r` therefore uses a local fixed wrapper that follows the package logic:

- split peaks into cumulative-logit and logit branches using the PACS rule;
- convert each branch block to dense matrices before calling
  `PACS::pacs_test_cumu()` or `PACS::pacs_test_logit()`;
- merge p-values by the original input peak order;
- avoid relying on the buggy convergence rowname merge.

The local wrapper does not modify the installed PACS package.

## Stage Conclusion

Notebook 1 benchmark reproduction is sufficient to support moving to mouse
kidney figure reproduction. The next stage should focus on producing clear
figures from the available kidney data and PACS results, while marking any
panels that require external annotations, motif databases, genome tracks, or
the original author plotting scripts.
