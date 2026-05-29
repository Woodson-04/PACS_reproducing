# Official Author Code Search Report

Date: 2026-05-28

Scope: mouse kidney UMAP / batch-effect correction code relevant to the PACS
paper. No heavy matrix analysis was run. The original GSE157079 `.gz` files
under `/home/woodson/biostatistic` were not modified, decompressed, moved, or
deleted.

## Repository Retrieval Status

Target repository:

```text
https://github.com/Zhen-Miao/PACS
```

Requested local path:

```text
/home/woodson/PACS_reproducing/external_reference/PACS_official
```

Status:

- Directory was created under the current project, but the repository was not
  successfully cloned.
- `git clone --depth 1 https://github.com/Zhen-Miao/PACS ./external_reference/PACS_official`
  failed with:

```text
fatal: could not create work tree dir '.\external_reference\PACS_official': Permission denied
```

- A fallback `git init/fetch` attempt inside the created directory also failed
  because `.git` metadata could not be created/written on the SSHFS path.
- `git ls-remote https://github.com/Zhen-Miao/PACS ...` failed with:

```text
fatal: unable to access 'https://github.com/Zhen-Miao/PACS/':
Failed to connect to 127.0.0.1 port 9 after 2039 ms: Could not connect to server
```

Target repository:

```text
https://github.com/Zhen-Miao/dev-kidney-snATAC
```

Requested local path:

```text
/home/woodson/PACS_reproducing/external_reference/dev-kidney-snATAC
```

Status:

- Directory was created under the current project, but the repository was not
  successfully cloned.
- `git ls-remote https://github.com/Zhen-Miao/dev-kidney-snATAC ...` failed
  with the same `127.0.0.1 port 9` connection error.

Interpretation:

- Local git clone/search could not be completed in the current Windows SSHFS
  environment.
- The likely causes are SSHFS work-tree write limitations plus a local git proxy
  configuration pointing to `127.0.0.1:9`.
- No official source files were overwritten or merged into this project.

## Online / Local References Inspected

### PACS package website

URL:

```text
https://zhen-miao.github.io/PACS/
```

Relevant points:

- The PACS package website states that PACS supports feature-level batch effect
  correction enabled by the statistical test framework.
- The example PACS call uses `pacs_test_sparse()` with covariates such as
  `cell_type` and `batch`.
- The website links to PACS vignettes, including Notebook 1.

Relevant files/pages found:

```text
PACS package website
PACS vignette: Vignette_1_PACS_cell_type_annotation.html
PACS GitHub Notebook 1:
https://github.com/Zhen-Miao/PACS/blob/main/vignettes/Notebook_1_Test_For_Sens_Spec_real_kidney_data.ipynb
```

Why this matters:

- Confirms the package-level functions and high-level workflow.
- Does not provide an explicit mouse kidney batch-effect UMAP reproduction
  script.

### PACS Nature Communications paper

URL:

```text
https://www.nature.com/articles/s41467-024-55580-5
```

Relevant text found:

- The adult kidney dataset contains strong batch effects.
- PACS detects peaks significantly affected by batch.
- Batch-effect peaks are removed from the feature set.
- Signac is then applied to both the original data and the batch-effect
  corrected data without any other batch correction.
- Figure 3 panels are described as UMAP plots constructed with all features or
  after excluding features significantly affected by batch effects, colored by
  batch labels or cell types.

Key method details from the paper:

```text
Signac workflow
TF-IDF without feature selection (min.cutoff = 'q0')
SVD dimensionality reduction
clustering and UMAP using dimensions 2-30
first LSI dimension usually reflects sequencing depth and is omitted
sample and cell type labels from the original publication annotations
```

Why this matters:

- This is the clearest available author description of how the PACS paper mouse
  kidney UMAP panels were generated.
- It strongly supports the conclusion that the GEO UMAP coordinate table is only
  an overview object and is not sufficient for reproducing the PACS paper's
  all-feature vs batch-effect-filtered UMAP panels.

### Original mouse kidney paper / GSE157079

URL:

```text
https://www.nature.com/articles/s41467-021-22266-1
```

Relevant points:

- The original mouse kidney paper states that code is available at:

```text
https://github.com/Zhen-Miao/dev-kidney-snATAC
```

- GSE157079 is the source of the downloaded snATAC data.

Why this matters:

- This repository likely contains scripts for the original GSE157079 data
  processing and UMAP generation.
- Because cloning failed in the current environment, no explicit script from
  this repository could be inspected locally in this pass.

## Relevant Local Files

Current project files:

```text
scripts/mouse_kidney_figures/00_inspect_gse157079.R
scripts/mouse_kidney_figures/00_prepare_gse157079_metadata.R
scripts/mouse_kidney_figures/03_gse157079_umap_plots.R
results/mouse_kidney_figures/gse157079_inspection_report.md
results/mouse_kidney_figures/gse157079_metadata_merged.csv
results/mouse_kidney_figures/gse157079_peak_list_preview.csv
```

Why these matter:

- They confirm the local GSE157079 inputs and support the current precomputed
  GEO UMAP overview figures.
- They do not yet implement author-style all-feature/batch-effect-filtered
  UMAP reconstruction.

## GSE157079 Matrix Header

From the existing inspection report, the large matrix is MatrixMarket-like:

```text
%%MatrixMarket matrix coordinate integer general
28316 300755 166121193
```

Interpretation:

- 28,316 cells.
- 300,755 peaks.
- 166,121,193 non-zero entries.
- The file is large and should be handled by a dedicated sparse-matrix loading
  script, not by the lightweight metadata/UMAP scripts.

## Author UMAP Code Findings

Explicit PACS paper mouse kidney batch-effect UMAP script:

```text
Not found in this pass.
```

Reason:

- Local clone/search of official repositories failed.
- The PACS package website and Nature paper provide method descriptions but not
  a ready-to-run script for Figure 3 UMAP reconstruction.

Explicit original GSE157079 UMAP script in `dev-kidney-snATAC`:

```text
Unclear / not inspected.
```

Reason:

- The original mouse kidney paper links to `dev-kidney-snATAC`, but the repo
  could not be cloned or searched locally in this environment.

## Method Clues Found

Normalization:

```text
Signac TF-IDF
```

Feature selection:

```text
For the PACS paper kidney UMAP, Signac was run without feature selection:
min.cutoff = 'q0'
```

Dimensionality reduction:

```text
SVD / LSI
```

UMAP:

```text
Run UMAP using dimensions 2-30
```

Batch-effect peak identification:

```text
PACS differential test module
two-sided test
FDR correction for multiple testing
features significantly affected by batch are removed for corrected UMAP
```

Likely PACS model form for batch-effect peak detection:

```text
formula_full includes batch and relevant biological covariates such as cell type
formula_null removes the batch term while retaining covariates that should be adjusted
```

Exact formulas and parameters:

```text
Not recovered from author code in this pass.
```

## Data Sufficiency Note

### A. Is the single large matrix file alone enough?

No. The matrix alone is not enough. It must be used together with:

```text
GSE157079_snATAC_metadata.csv.gz
GSE157079_snATAC_peak_list.csv.gz
```

The metadata are needed for sample/batch labels and cell type labels, and the
peak list is needed for genomic identities of the features.

### B. Are the four GSE157079 files likely sufficient as raw inputs for an author-style reconstruction?

Likely yes for a core UMAP reconstruction. The four files provide:

```text
cell-by-peak counts
sample/batch labels
cell type labels
peak identities
precomputed GEO UMAP coordinates for reference
```

They should support a Signac-style all-feature UMAP and a PACS-filtered UMAP
once the sparse matrix is loaded correctly.

### C. Are they sufficient for exact author reproduction without author code?

No. Exact reproduction also requires:

- exact preprocessing choices;
- Signac/Seurat versions and parameters;
- PACS batch-test formula choices;
- FDR threshold choices;
- whether binary or count/PIC values were used at each step;
- random seeds and UMAP settings;
- any quality-control filtering applied before the deposited matrix.

## Recommended Next Step

Do not run full UMAP recomputation yet. First create a small, explicit
matrix-loading smoke test that:

1. reads only the MatrixMarket header and confirms dimensions;
2. loads a small subset or uses a sparse MatrixMarket reader in a controlled
   output directory;
3. checks that metadata rows align with matrix rows by `row_index`;
4. checks that peak list rows align with matrix columns;
5. defines an author-style reconstruction plan:
   - all-feature Signac TF-IDF/SVD/UMAP using dims 2-30;
   - PACS batch-effect peak test;
   - remove FDR-significant batch-effect peaks;
   - repeat Signac TF-IDF/SVD/UMAP using remaining peaks.

Only after this smoke test should the full 28,316 x 300,755 sparse matrix be
processed.

## Update: Local `dev-kidney-snATAC` Clone Inspected

Date: 2026-05-28

The user cloned the original developmental kidney snATAC repository under the
read-only reference area:

```text
/home/woodson/biostatistic/pacs/dev-kidney-snATAC*/dev-kidney-snATAC
```

This directory was inspected read-only. No files under
`/home/woodson/biostatistic` were modified.

### Repository Contents

The clone contains a small set of R scripts:

```text
R/combining all batches__create object.R
R/all_pairwise DAR analysis.R
R/DAP between P0 and 3 week and 8 week.R
R/distal element prediction using P0 data.R
R/stratify genomic elements.R
R/scRNA-seq__*.R
README.md
```

### Most Relevant File

```text
/home/woodson/biostatistic/pacs/dev-kidney-snATAC*/dev-kidney-snATAC/R/combining all batches__create object.R
```

This appears to be the key original GSE157079-style snATAC processing script.
It is not a PACS paper script, but it explains how the original developmental
kidney snATAC object and GEO-style UMAP were likely produced.

Important steps found:

- Combines five SnapATAC `.snap` batches:
  `90025`, `90026`, `90028`, `90029`, and `batch_1`.
- Filters cells using per-batch barcode list files.
- Adds 5 kb bin matrices with `addBmatToSnap(bin.size = 5000)`.
- Combines shared bins across batches with `snapRbind`.
- Binarizes the bin matrix with `makeBinary`.
- Removes blacklist regions.
- Removes `random` and `chrM` regions.
- Removes highly conserved/high-coverage bins using the 95th percentile of
  `log10(colSums(bmat) + 1)`.
- Runs SnapATAC diffusion maps:
  `runDiffusionMaps(input.mat = "bmat", num.eigs = 50)`.
- Uses dimensions 1:20 for KNN, clustering, and UMAP.
- Runs UMAP with SnapATAC:
  `runViz(method = "umap", eigs.dims = 1:20, seed.use = 10)`.
- Produces batch-colored UMAP before correction.
- Corrects batch effects with Harmony:
  `runHarmony(eigs.dim = 1:20, meta_data = x.sp@sample)`.
- Runs KNN, clustering, and UMAP again after Harmony.
- Produces batch-colored UMAP after correction.
- Saves:
  `all_comb snapATAC before batch correction.RDS` and
  `all_comb snapATAC after batch_correction.RDS`.

This gives a concrete reconstruction target for the deposited precomputed GEO
UMAP: it is likely SnapATAC diffusion-map/UMAP based, with optional Harmony
correction, rather than Signac TF-IDF/LSI.

### DAR / Motif File

```text
/home/woodson/biostatistic/pacs/dev-kidney-snATAC*/dev-kidney-snATAC/R/all_pairwise DAR analysis.R
```

This script starts from:

```text
all_comb snapATAC after batch correction with peak info new clu.RDS
```

and uses SnapATAC `findDAR` on `pmat` for pairwise cell type/cluster DAR
analysis, followed by BH FDR correction and HOMER motif enrichment.

This is useful for later mouse-kidney DAR/motif figure reproduction, but it is
not the PACS batch-effect feature-removal UMAP script.

### Explicit PACS Paper UMAP Code?

Not found in the local `dev-kidney-snATAC` clone.

The clone contains original mouse kidney SnapATAC/Harmony processing code for
the 2021 developmental kidney atlas, not the later PACS paper workflow that
constructs UMAPs with all features and after excluding PACS batch-effect
features.

### Interpretation for the PACS Reproduction Project

There are now two distinct UMAP targets:

1. **GEO/original atlas-style UMAP**
   - Best matched by the `dev-kidney-snATAC` script.
   - Likely SnapATAC bin matrix, diffusion maps, `runViz(method = "umap")`,
     and Harmony for batch-corrected coordinates.
   - Uses raw `.snap` objects and barcode lists that are not present in the
     four GEO `.gz` files.

2. **PACS paper batch-feature-removal UMAP**
   - Best matched by the PACS paper Methods, not by the developmental kidney
     repository.
   - Described as applying a UMAP workflow to all features, then excluding
     features significantly affected by batch effects.
   - The exact author script is still not found locally.
   - The four GSE157079 files are likely enough for an author-style
     reconstruction, but exact reproduction still requires the missing PACS
     figure script or precise parameters.

### Recommended Next Step After This Update

Do not rewrite the current precomputed UMAP plots. They are useful as GEO
overview panels.

For the next reconstruction stage, create a light smoke-test script that only
validates loading and alignment of:

```text
GSE157079_snATAC_cell_by_peak_matrix.txt.gz
GSE157079_snATAC_metadata.csv.gz
GSE157079_snATAC_peak_list.csv.gz
```

Then choose one of two explicit branches:

- **SnapATAC atlas branch:** mimic `combining all batches__create object.R`
  as closely as possible from the deposited matrix, recognizing that the
  original `.snap` objects and barcode-list filtering are unavailable.
- **PACS paper branch:** build the all-feature UMAP and PACS-filtered UMAP from
  the deposited cell-by-peak matrix, using metadata `sample` as the batch
  variable and documenting every parameter.
