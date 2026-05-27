# Mouse Kidney Figure Source Report

Date: 2026-05-27

## Search Scope

Local locations searched:

```text
/home/woodson/PACS_reproducing
/home/woodson/biostatistic
/home/woodson/biostatistic/pacs
/home/woodson/biostatistic/pacs/PACS_data
```

Keywords searched included:

```text
kidney, mouse, Figure 4, fig4, adult, PT, LOH, DAR, differential,
PACS_test, FindMarkersPACS, cell_type, celltype, UMAP, TSNE, motif,
enrichment
```

Online search terms attempted through the available web search interface:

```text
Zhen-Miao/PACS Figure 4 kidney PACS
PACS mouse kidney dataset analysis code
Model-based compound hypothesis testing PACS Figure 4 kidney
```

The PACS paper page is:

```text
https://www.nature.com/articles/s41467-024-55580-5
```

Notebook 1 cites the original mouse kidney dataset paper:

```text
Miao, Z., Balzer, M.S., Ma, Z. et al.
Single cell regulatory landscape of the mouse kidney highlights cellular
differentiation programs and disease targets. Nature Communications 12, 2277
(2021). https://doi.org/10.1038/s41467-021-22266-1
```

## Code Found

Current project:

```text
/home/woodson/PACS_reproducing/q.r
/home/woodson/PACS_reproducing/Notebook_1_Test_For_Sens_Spec_real_kidney_data.ipynb
/home/woodson/PACS_reproducing/baseline_methods_notebook1.R
```

Old project/reference directory:

```text
/home/woodson/biostatistic/pacs/run_PACS.r
/home/woodson/biostatistic/pacs/align_pvals.r
/home/woodson/biostatistic/pacs/my_methods.r
/home/woodson/biostatistic/pacs/Notebook_1_Test_For_Sens_Spec_real_kidney_data.ipynb
```

Notes:

- `run_PACS.r` is an older local attempt at the PT/LOH PACS benchmark, not a
  figure-generation script.
- `my_methods.r` is an older simplified baseline substitute and is intentionally
  not used for the current reproduction.
- No direct author script for PACS paper mouse kidney Figure 3/Figure 4 style
  panels was found locally.

## Data Found

The available PACS data are:

```text
/home/woodson/biostatistic/pacs/PACS_data/data_for_test_for_t1e_power.rdata
/home/woodson/biostatistic/pacs/PACS_data/kidney_features_to_keep.rds
/home/woodson/biostatistic/pacs/PACS_data/r_by_ct_est_kidney_adult.rds
```

These contain the matrix and labels used by Notebook 1:

- `pmats`: original cells x peaks count matrix.
- `x.sp_cluster2`: cell-type labels including PT and LOH.
- `kidney_features_to_keep`: feature subset used by Notebook 1.
- `r_by_ct_est$q_vec_new`: per-cell capture-rate/depth parameter estimates.

The completed large Notebook 1 benchmark result is:

```text
/home/woodson/PACS_reproducing/results/20260526_2318_large_baseline/pacs_kidney_notebook1_result.rds
/home/woodson/PACS_reproducing/results/20260526_2318_large_baseline/summary.csv
```

## Direct Figure Script Availability

No ready-to-run author script for the mouse kidney figure panels was found in
the accessible local files.

The accessible PACS_data can support matrix-level overview and PT/LOH
differential accessibility figures, but not every publication-quality biological
annotation panel.

## Figures That Can Be Built Now

Based on current data and Notebook 1 results:

- Cell type composition bar plot.
- PT vs LOH depth/capture-rate distribution.
- PACS permuted-label QQ plot.
- PACS Type I error/power benchmark bar plot.
- PT vs LOH DAR volcano plot using PACS p-values from the actual-label result,
  with an effect-size proxy computed from PT/LOH accessibility.
- Top DAR peak heatmap using sampled PT/LOH cells and normalized accessibility.

## Figures Requiring Extra Inputs

The following require additional annotation files, external packages, or the
original author processing objects:

- DAR peaks linked genes: requires peak-to-gene or nearest-gene annotation.
- Motif enrichment: requires motif databases and genomic ranges for peaks.
- Genome browser style peak tracks: requires genome coordinates, fragment files
  or pseudobulk tracks, and genome annotation.
- Gene activity/regulatory program figures: requires gene activity matrices,
  RNA expression data, or the original mouse kidney analysis object.
- UMAP/tSNE panels: require reduced-dimension coordinates or enough metadata
  and preprocessing choices to reconstruct them.

## Recommended Next Figures

Priority 1:

1. Overview plots from `01_overview_plots.R`.
2. PACS PT vs LOH DAR volcano plot.
3. Top DAR heatmap.

Priority 2:

1. Add peak coordinate parsing/annotation if row names contain genomic
   coordinates.
2. Add motif enrichment after confirming a genome build and motif database.
3. Add browser-style tracks only if fragment files or pseudobulk coverage files
   become available.

This route keeps the next stage honest: first produce polished figures that are
fully supported by local data, then expand toward biological annotation panels
once the needed inputs are identified.
