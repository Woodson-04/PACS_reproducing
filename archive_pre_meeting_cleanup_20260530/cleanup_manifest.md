# Cleanup Manifest

Date: 2026-05-27

Scope: `/home/woodson/PACS_reproducing` only.

This manifest was written before cleanup. The original data directories under
`/home/woodson/biostatistic`, including `PACS_data` and `GSE157079`, must remain
read-only and untouched.

## Keep

Core code and reports:

```text
q.r
baseline_methods_notebook1.R
readme.md
Notebook_1_Test_For_Sens_Spec_real_kidney_data.ipynb
notebook1_reproduction_report.md
mouse_kidney_figure_plan.md
mouse_kidney_figure_source_report.md
baseline_official_methods_review.md
cleanup_manifest.md
```

Figure workspace:

```text
scripts/mouse_kidney_figures/
figures/mouse_kidney/
results/mouse_kidney_figures/
```

Reference manuscript PDF:

```text
Depth-corrected multi-factor dissection of chromatin accessibility for scATAC-seq data with PACS (1).pdf
```

## Representative Results To Keep

```text
results/20260526_2318_large_baseline
results/kidney_notebook1_20260526_172547
results/kidney_notebook1_20260526_221904
results/kidney_notebook1_20260526_205633
results/mouse_kidney_figures
```

Reasons:

- `20260526_2318_large_baseline`: final large Notebook 1 benchmark with
  clean-room baselines.
- `kidney_notebook1_20260526_172547`: successful full-scale PACS-only run.
- `kidney_notebook1_20260526_221904`: successful medium baseline run after
  snapATAC edgeR alignment.
- `kidney_notebook1_20260526_205633`: successful small baseline smoke test.
- `mouse_kidney_figures`: destination for lightweight GSE157079 inspection and
  metadata outputs.

## Delete

Early failed, duplicate, or superseded Notebook 1 result directories:

```text
results/kidney_notebook1_20260525_165546
results/kidney_notebook1_20260525_183621
results/kidney_notebook1_20260525_211229
results/kidney_notebook1_20260526_091336
results/kidney_notebook1_20260526_092605
results/kidney_notebook1_20260526_171139
results/kidney_notebook1_20260526_210731
```

Standalone debug/search artifacts:

```text
results/pacs_sparse_mixed_error_terminal.log
baseline_source_search_report.md
```

Deletion reason:

- PACS mixed-branch debugging is now summarized in
  `notebook1_reproduction_report.md`.
- Baseline source search conclusions are summarized in README and formal
  reports.
- The final large benchmark and selected small/medium/full PACS-only runs remain
  as representative evidence.

## Do Not Touch

Uncertain or external paths:

```text
/home/woodson/biostatistic
/home/woodson/biostatistic/pacs
/home/woodson/biostatistic/pacs/PACS_data
/home/woodson/biostatistic/pacs/GSE157079
installed PACS package directories
```

No git staging, committing, or pushing should be performed during cleanup.

## Cleanup Execution Status

Deleted in this pass:

```text
results/pacs_sparse_mixed_error_terminal.log
baseline_source_search_report.md
```

Recursive deletion of the superseded result directories was prepared with
path-safety checks, but the local tool approval/review timed out twice before
the command could run. These directories remain pending manual deletion:

```text
results/kidney_notebook1_20260525_165546
results/kidney_notebook1_20260525_183621
results/kidney_notebook1_20260525_211229
results/kidney_notebook1_20260526_091336
results/kidney_notebook1_20260526_092605
results/kidney_notebook1_20260526_171139
results/kidney_notebook1_20260526_210731
```
