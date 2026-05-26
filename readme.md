# PACS Notebook 1 Reproduction

This project reproduces the PACS Notebook 1 workflow:

```text
Notebook_1_Test_For_Sens_Spec_real_kidney_data.ipynb
Test PACS sensitivity and specificity on real kidney data
```

The original data and old reference files live under `/home/woodson/biostatistic` and are treated as read-only. New code, logs, reports, and results belong in `/home/woodson/PACS_reproducing`.

## Current PACS-Only Result

A full PACS-only run completed successfully with:

```text
our t1e    = 0.04008
our t1e_sd = 0.00186735106501161
our power  = 1
```

The `power = 1` value is expected for `--run_baselines FALSE` because the pseudo-true union contains only `our`. It is not comparable to the complete Notebook 1 power comparison, where the pseudo-true union uses `our`, `seurat`, `archR`, and `snapATAC` and excludes `fisher`.

## Baseline Strategy

Do not use `/home/woodson/biostatistic/pacs/my_methods.r`; it is an old simplified substitute and is not part of this reproduction path.

Baseline loading in `q.r` uses this priority when `--run_baselines TRUE`:

1. `./other_methods_for_differential_updated.R`
   - Original author notebook baseline file, if it is recovered and copied into this project.
2. `./baseline_methods_notebook1.R`
   - Clean-room reimplemented baselines.
   - Not the original author baseline file.
   - Results must be reported as reimplemented/approximate baselines.

The search status is documented in:

```text
baseline_source_search_report.md
```

## Clean-Room Baseline Alignment Notes

The original `other_methods_for_differential_updated.R` has not been found. Until that file is recovered, baseline results must be labeled as `clean_room_reimplemented` and should not be described as the author's original baseline output.

Current clean-room baseline alignment:

- `seurat`: logistic regression likelihood-ratio test with total depth adjustment.
  - Full model: `group ~ log_depth + peak_value`
  - Null model: `group ~ log_depth`
  - Uses raw peak counts as `peak_value`.
- `snapATAC`: edgeR-based fixed-dispersion test, matching the notebook comment `snapATAC-- edgeR method`.
  - Uses binarized peak x cell matrices.
  - Calls `edgeR::DGEList`, `edgeR::calcNormFactors`, and `edgeR::exactTest`.
  - Uses `dispersion = bcv^2`, with notebook default `bcv = 0.4`.
  - Requires the `edgeR` package; it fails clearly if `edgeR` is unavailable.
- `fisher`: binary 2x2 two-sided Fisher exact test.
  - Rows are accessible/inaccessible.
  - Columns are positive/negative group.
- `archR`: approximate depth-aware implementation.
  - Performs nearest-neighbor depth matching between groups.
  - Tests log1p depth-normalized counts with Wilcoxon rank-sum.
  - This is not ArchR's official implementation and should be treated as an approximation unless the original helper is recovered.

## Commands

PACS-only full run:

```bash
cd /home/woodson/PACS_reproducing
Rscript q.r \
  --data_dir /home/woodson/biostatistic/pacs/PACS_data \
  --n_repeat 5 \
  --n_cell_sample 500 \
  --n_features_sample 10000 \
  --run_baselines FALSE
```

Baseline smoke test:

```bash
cd /home/woodson/PACS_reproducing
Rscript q.r \
  --data_dir /home/woodson/biostatistic/pacs/PACS_data \
  --n_repeat 1 \
  --n_cell_sample 50 \
  --n_features_sample 100 \
  --run_baselines TRUE
```

Baseline medium test:

```bash
cd /home/woodson/PACS_reproducing
Rscript q.r \
  --data_dir /home/woodson/biostatistic/pacs/PACS_data \
  --n_repeat 1 \
  --n_cell_sample 100 \
  --n_features_sample 1000 \
  --run_baselines TRUE
```

Full baseline runs may be slow because the clean-room Seurat baseline fits per-peak logistic regression models and the ArchR approximation performs per-peak Wilcoxon tests.
