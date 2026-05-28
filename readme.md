# PACS Notebook 1 Reproduction

This project reproduces the PACS Notebook 1 workflow:

```text
Notebook_1_Test_For_Sens_Spec_real_kidney_data.ipynb
Test PACS sensitivity and specificity on real kidney data
```

The original data and old reference files live under `/home/woodson/biostatistic` and are treated as read-only. New code, logs, reports, and results belong in `/home/woodson/PACS_reproducing`.

## Current Focus

Notebook 1 benchmark reproduction is now complete enough for the project goal.
The PACS main result closely matches the author notebook, and remaining
baseline differences are documented as clean-room reimplementation differences.

The project is now moving to mouse kidney figure reproduction. The first figure
stage uses existing PACS_data and the completed Notebook 1 benchmark output to
create overview/QC/benchmark plots before moving to DAR volcano and heatmap
figures.

Key reports:

```text
notebook1_reproduction_report.md
mouse_kidney_figure_source_report.md
mouse_kidney_figure_plan.md
```

## Current PACS-Only Result

A full PACS-only run completed successfully with:

```text
our t1e    = 0.04008
our t1e_sd = 0.00186735106501161
our power  = 1
```

The `power = 1` value is expected for `--run_baselines FALSE` because the pseudo-true union contains only `our`. It is not comparable to the complete Notebook 1 power comparison, where the pseudo-true union uses `our`, `seurat`, `archR`, and `snapATAC` and excludes `fisher`.

## Completed Large Benchmark

The completed large baseline run is:

```text
results/20260526_2318_large_baseline
```

Summary:

```text
our      t1e = 0.04008  power = 0.83337
seurat   t1e = 0.06342  power = 0.82344
archR    t1e = 0.04096  power = 0.67437
snapATAC t1e = 0.01810  power = 0.76094
fisher   t1e = 0.02208  power = 0.76630
```

The PACS main result is close to the author notebook. Seurat and Fisher are also
close. ArchR and snapATAC should be interpreted through the clean-room baseline
notes because the original author baseline helper is still missing.

## Results Directory Naming

New runs are written under:

```text
results/YYYYMMDD_HHMM_<scale>_<mode>
```

Scale is assigned from the run parameters:

- `small`: `n_repeat = 1`, `n_cell_sample = 50`, `n_features_sample = 100`
- `medium`: `n_repeat = 1`, `n_cell_sample = 100`, `n_features_sample = 1000`
- `large`: `n_repeat = 5`, `n_cell_sample = 500`, `n_features_sample = 10000`
- `custom`: any other combination

Mode is:

- `pacs_only`: `run_baselines = FALSE`
- `baseline`: `run_baselines = TRUE`

If a target directory already exists, `q.r` appends `_v2`, `_v3`, and so on.

## Mouse Kidney Figure Workspace

New figure work is organized under:

```text
figures/mouse_kidney/
scripts/mouse_kidney_figures/
results/mouse_kidney_figures/
```

The first script is:

```text
scripts/mouse_kidney_figures/01_overview_plots.R
```

It generates:

```text
figures/mouse_kidney/cell_type_counts_barplot.png
figures/mouse_kidney/pt_loh_depth_distribution.png
figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png
figures/mouse_kidney/pacs_permuted_qq_plot.png
```

## GSE157079 Intake

The downloaded GSE157079 files are treated as read-only source data:

```text
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_UMAP_coordinates.csv.gz
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_metadata.csv.gz
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_peak_list.csv.gz
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz
```

The project does not modify or extract files in that directory. The first
GSE157079 stage saves only lightweight inspection reports, merged metadata,
peak previews, and UMAP figures under this repository. The large
`cell_by_peak_matrix` is not read fully by the inspection or UMAP scripts; it
needs a dedicated matrix-processing script later.

GSE157079 metadata column handling:

- The UMAP CSV has an unnamed first column and no barcode; it is standardized as
  `row_index`, `umap_1`, `umap_2`.
- The metadata CSV is standardized as `row_index`, `cell_barcode`, `sample`,
  `cell_type`.
- UMAP and metadata are merged by `row_index`.
- `cell_barcode` comes from the metadata `barcodes` column.
- `sample` comes from the metadata `samples` column.
- `cell_type` comes from the metadata `clusters` column.

GSE157079 scripts:

```text
scripts/mouse_kidney_figures/00_inspect_gse157079.R
scripts/mouse_kidney_figures/00_prepare_gse157079_metadata.R
scripts/mouse_kidney_figures/03_gse157079_umap_plots.R
```

## Baseline Strategy

Do not use `/home/woodson/biostatistic/pacs/my_methods.r`; it is an old simplified substitute and is not part of this reproduction path.

Baseline loading in `q.r` uses this priority when `--run_baselines TRUE`:

1. `./other_methods_for_differential_updated.R`
   - Original author notebook baseline file, if it is recovered and copied into this project.
2. `./baseline_methods_notebook1.R`
   - Clean-room reimplemented baselines.
   - Not the original author baseline file.
   - Results must be reported as reimplemented/approximate baselines.

The original baseline search conclusion is now summarized in
`notebook1_reproduction_report.md` and `baseline_official_methods_review.md`.

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

More detailed official-method alignment notes are in:

```text
baseline_official_methods_review.md
```

## Current Medium Baseline Result

A medium clean-room baseline run completed with:

```text
n_repeat          = 1
n_cell_sample     = 100
n_features_sample = 1000
run_baselines     = TRUE
baseline_source   = clean_room_reimplemented
```

Summary:

```text
our      t1e = 0.02500  power = 0.7306
seurat   t1e = 0.06636  power = 0.7763
archR    t1e = 0.02774  power = 0.3699
snapATAC t1e = 0.00100  power = 0.6027
fisher   t1e = 0.00400  power = 0.6027
```

Because this run has only one repeat, the standard deviations are not
informative. The snapATAC baseline is now implemented as edgeR with fixed
`bcv = 0.4`, matching the notebook comment more closely; its Type I error is
conservative in this medium run. ArchR remains the main uncertainty because the
clean-room implementation only has matrix data and depth matching, not a full
ArchRProject with Arrow files and ArchR's full bias matching.

The next large baseline run is useful, but results must still be reported as
`clean_room_reimplemented`, not as the original author baseline.

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

Baseline intermediate custom run:

```bash
cd /home/woodson/PACS_reproducing
Rscript q.r \
  --data_dir /home/woodson/biostatistic/pacs/PACS_data \
  --n_repeat 2 \
  --n_cell_sample 300 \
  --n_features_sample 3000 \
  --run_baselines TRUE
```

Baseline large run:

```bash
cd /home/woodson/PACS_reproducing
Rscript q.r \
  --data_dir /home/woodson/biostatistic/pacs/PACS_data \
  --n_repeat 5 \
  --n_cell_sample 500 \
  --n_features_sample 10000 \
  --run_baselines TRUE
```

Full baseline runs may be slow because the clean-room Seurat baseline fits per-peak logistic regression models and the ArchR approximation performs per-peak Wilcoxon tests.

Mouse kidney overview figures:

```bash
cd /home/woodson/PACS_reproducing
Rscript scripts/mouse_kidney_figures/01_overview_plots.R \
  --data_dir /home/woodson/biostatistic/pacs/PACS_data \
  --notebook1_result_dir /home/woodson/PACS_reproducing/results/20260526_2318_large_baseline \
  --output_dir /home/woodson/PACS_reproducing/figures/mouse_kidney
```

GSE157079 inspection:

```bash
cd /home/woodson/PACS_reproducing
Rscript scripts/mouse_kidney_figures/00_inspect_gse157079.R \
  --gse_dir /home/woodson/biostatistic/pacs/GSE157079 \
  --out_dir /home/woodson/PACS_reproducing/results/mouse_kidney_figures
```

GSE157079 metadata preparation:

```bash
cd /home/woodson/PACS_reproducing
Rscript scripts/mouse_kidney_figures/00_prepare_gse157079_metadata.R \
  --gse_dir /home/woodson/biostatistic/pacs/GSE157079 \
  --out_dir /home/woodson/PACS_reproducing/results/mouse_kidney_figures
```

GSE157079 UMAP plots:

```bash
cd /home/woodson/PACS_reproducing
Rscript scripts/mouse_kidney_figures/03_gse157079_umap_plots.R \
  --metadata_csv /home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv \
  --out_dir /home/woodson/PACS_reproducing/figures/mouse_kidney
```
