# P56 PACS Batch-Effect Filtering Parameter Tuning Plan

This is a planning document only. No PACS, MatrixMarket streaming, TF-IDF, LSI,
or UMAP analysis was rerun for this report.

## Current Status

Current P56 run:

- `n_top_peaks = 10000`
- `max_pacs_peaks = 5000`
- `fdr_cutoff = 0.05`
- tested peaks: `5000`
- significant batch peaks: `3140`
- retained tested peaks: `1860`

Quantification from `07c_quantify_p56_batch_mixing.R`:

| metric | before | after | direction |
|---|---:|---:|---|
| batch silhouette | 0.18305 | 0.02661 | decreased, better mixing |
| normalized batch entropy | 0.30363 | 0.69564 | increased, better mixing |
| same-batch fraction | 0.98658 | 0.69693 | decreased, better mixing |
| batch prediction accuracy | 0.76984 | 0.58995 | decreased, better mixing |
| cell type silhouette | 0.17559 | 0.15933 | mildly decreased |
| same-celltype fraction | 0.84885 | 0.66504 | substantially decreased |

Interpretation:

- PACS filtering clearly improves batch mixing by quantitative metrics.
- Correction is **partial, not complete**: residual batch separation remains in
  the after-filtering UMAP.
- The drop in same-celltype fraction is large enough to raise a possible
  overcorrection / biological-structure-loss warning.
- Therefore, parameter scanning is appropriate, but it should be staged and
  judged by quantitative metrics, not by UMAP visual appearance alone.

## Stage A: Use Existing Results

Use the current run as the baseline for tuning.

Actions:

1. Keep the existing `07c` metrics as the reference point.
2. Classify the current correction as:
   - batch mixing improvement: strong;
   - residual batch separation: still present;
   - cell type preservation: weakened and needs monitoring.
3. Do not tune figures visually before deciding the next PACS setting.

Existing output files:

```text
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap/p56_batch_mixing_metrics_summary.csv
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap/p56_batch_mixing_by_celltype.csv
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap/p56_batch_mixing_quantification_report.md
```

## Stage B: Increase Tested Peaks

Current run tested only the top `5000` peaks, and `3140` of them were called
significant batch-associated peaks. Untested high-detection peaks may still
carry batch signal.

Recommended staged runs:

1. `n_top_peaks = 10000`, `max_pacs_peaks = 10000`, `fdr_cutoff = 0.05`
2. `n_top_peaks = 20000`, `max_pacs_peaks = 10000`, `fdr_cutoff = 0.05`
3. `n_top_peaks = 20000`, `max_pacs_peaks = 20000`, `fdr_cutoff = 0.05`

Rationale:

- Run 1 is the cleanest next comparison because it keeps the same top-peak pool
  size but tests all `10000` selected peaks.
- Run 2 asks whether a larger peak universe improves the before/after
  reconstruction while keeping PACS cost controlled.
- Run 3 is the fuller P56 top-peak version, but should only be run if runtime
  and memory are acceptable.

Recommended next actual run:

```bash
cd /home/woodson/PACS_reproducing

Rscript scripts/mouse_kidney_figures/07_p56_pacs_batch_filter_umap.R \
  --gse_dir /home/woodson/biostatistic/pacs/GSE157079 \
  --metadata_csv /home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv \
  --out_dir /home/woodson/PACS_reproducing/results/mouse_kidney_figures \
  --fig_dir /home/woodson/PACS_reproducing/figures/mouse_kidney \
  --n_top_peaks 10000 \
  --max_pacs_peaks 10000 \
  --fdr_cutoff 0.05 \
  --seed 1 \
  --chunk_lines 100000 \
  --progress_every 5000000
```

Then quantify with:

```bash
cd /home/woodson/PACS_reproducing

Rscript scripts/mouse_kidney_figures/07c_quantify_p56_batch_mixing.R \
  --result_dir /home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap \
  --fig_dir /home/woodson/PACS_reproducing/figures/mouse_kidney \
  --k 30 \
  --seed 1
```

Important implementation note:

- To compare settings cleanly, future runs should ideally write to setting-
  specific subdirectories or include setting tags in output names. The current
  `07_p56` script writes to a fixed output directory, so rerunning it will
  overwrite the current P56 result files unless the script is adjusted first.

## Stage C: Adjust FDR Threshold

If residual batch separation remains after increasing tested peaks, compare:

- `fdr_cutoff = 0.05`
- `fdr_cutoff = 0.10`
- `fdr_cutoff = 0.20`

Interpretation:

- `0.05` is the conservative default.
- `0.10` is a moderate exploratory threshold.
- `0.20` should be clearly labeled exploratory and potentially aggressive.

Warning:

The current run already removed many tested peaks and reduced same-celltype
fraction substantially. Therefore, increasing FDR may improve batch mixing but
may also worsen biological-structure preservation. Do not choose a higher FDR
unless quantitative metrics support it.

## Stage D: Cell-Type-Stratified Filtering

If global PACS still leaves strong batch separation, run batch peak testing
within major cell types separately.

Candidate major cell types:

- `PT`
- `PT2`
- `LOH`
- `stroma1` / `stroma2`
- `Endo`
- other sufficiently large types with both P56 batches represented

Plan:

1. Split P56 cells by major cell type.
2. Within each cell type, test batch effect without global cell type covariate.
3. Call batch-associated peaks per cell type.
4. Take the union of cell-type-specific batch-associated peaks.
5. Recompute filtered UMAP using the union-filtered feature set.

Rationale:

Batch-associated peaks may be cell-type-specific. A global model can dilute
cell-type-specific batch effects, especially if some cell types dominate the
dataset.

Risk:

This approach may remove more biologically meaningful cell-type-specific peaks.
It requires careful monitoring of cell type silhouette and same-celltype
fraction.

## Stage E: Feature Filtering Alternatives

Compare several filtering definitions:

1. Remove FDR-significant batch peaks only.
2. Remove top `N` batch peaks by p value:
   - top 3000;
   - top 5000;
   - top 8000.
3. Remove peaks with `FDR <= 0.10`.

Rationale:

- FDR-based filtering is statistically principled.
- Top-N filtering gives predictable retained feature counts and may be useful
  when FDR calls are too aggressive or too unstable.
- FDR `0.10` can be used as an exploratory middle ground.

## Stage F: Quantitative Selection Criterion

For every run, compute the same metrics with `07c`:

- batch silhouette decrease;
- normalized batch entropy increase;
- same-batch fraction decrease;
- batch prediction accuracy decrease;
- cell type silhouette preservation;
- same-celltype fraction preservation.

Preferred setting:

Choose the setting that improves batch mixing while preserving cell type
structure. A good setting should:

- reduce batch silhouette;
- increase normalized batch entropy;
- reduce same-batch fraction;
- reduce batch prediction accuracy;
- avoid a sharp drop in cell type silhouette;
- avoid a sharp drop in same-celltype fraction.

Practical warning:

Do not judge only by UMAP visual appearance. UMAP can exaggerate or hide
structure depending on local geometry. Use both quantitative metrics and visual
inspection, and avoid overcorrection.

## Current Recommendation

The project is suitable for a **small staged parameter scan**, not a large
uncontrolled grid search.

Recommended immediate next setting:

```text
n_top_peaks = 10000
max_pacs_peaks = 10000
fdr_cutoff = 0.05
```

Reason:

- Current run tested only `5000` peaks.
- `3140 / 5000` tested peaks were significant at FDR 0.05.
- Untested peaks may retain batch signal.
- Increasing tested peaks to `10000` is the most direct next step.
- Keep FDR at `0.05` first because current same-celltype fraction already
  dropped substantially after filtering.

Only after this run should the project consider more aggressive FDR thresholds
or cell-type-stratified PACS filtering.
