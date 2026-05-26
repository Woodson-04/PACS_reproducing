# Baseline Source Search Report

Date: 2026-05-26

## Target

Original author baseline file:

```text
other_methods_for_differential_updated.R
```

Notebook 1 raw URL:

```text
https://raw.githubusercontent.com/Zhen-Miao/PACS/main/analysis_codes/other_methods_for_differential_updated.R
```

## Local Search

Searched current project:

```text
/home/woodson/PACS_reproducing
```

Commands used:

```bash
rg --files | grep -E 'other_methods_for_differential_updated\.R|other_methods|baseline_methods'
rg -n 'other_methods_for_differential_updated|seurat_method2_subsample|snapATAC_method|archR_method|fisher_method' -S
```

Result: the original `other_methods_for_differential_updated.R` was not found.

Searched old project/reference area:

```text
/home/woodson/biostatistic
/home/woodson/biostatistic/pacs
```

Commands used:

```bash
rg --files | grep -E 'other_methods_for_differential_updated\.R|other_methods'
rg -n 'other_methods_for_differential_updated|seurat_method2_subsample|snapATAC_method|archR_method|fisher_method' -S
```

Result: the original file was not found. The old simplified file `my_methods.r` exists in `/home/woodson/biostatistic/pacs`, but it is not the author original and is intentionally not used.

Searched installed PACS package:

```text
/home/woodson/R/x86_64-pc-linux-gnu-library/4.6/PACS
```

Commands used:

```bash
rg --files | grep -E 'other_methods_for_differential_updated\.R|other_methods|seurat_method2_subsample|snapATAC_method|archR_method|fisher_method'
```

Result: the original baseline file and baseline functions were not found in the installed PACS package.

## Git History Search

Current project git history:

```bash
git log --all --name-only -- '*other_methods*'
git log --all --name-only -- '*.R'
git grep -n 'other_methods_for_differential_updated|seurat_method2_subsample|snapATAC_method|archR_method|fisher_method' $(git rev-list --all)
```

Result: git history references the notebook source URL and baseline function calls, but does not contain the original baseline R file.

Old `/home/woodson/biostatistic/pacs` directory:

```bash
git -C /home/woodson/biostatistic/pacs status --short
```

Result: not a git repository, so no git history was available there.

## Notebook Source Reference

`Notebook_1_Test_For_Sens_Spec_real_kidney_data.ipynb` sources:

```text
https://raw.githubusercontent.com/Zhen-Miao/PACS/main/analysis_codes/other_methods_for_differential_updated.R
```

The notebook also calls these baseline functions:

```text
seurat_method2_subsample
snapATAC_method
fisher_method
archR_method
```

## Online Search

Attempted searches:

```text
"other_methods_for_differential_updated.R"
"Zhen-Miao/PACS other_methods_for_differential_updated.R"
"seurat_method2_subsample"
"snapATAC_method" "archR_method" "fisher_method"
"Zhen-Miao/PACS analysis_codes/other_methods_for_differential_updated.R"
```

Result: no retrievable copy of the original `other_methods_for_differential_updated.R` was found from the available search interface. The GitHub repository page is reachable, but the exact raw file contents were not available through the current environment.

Manual URL to try in a browser or Linux terminal:

```text
https://raw.githubusercontent.com/Zhen-Miao/PACS/main/analysis_codes/other_methods_for_differential_updated.R
```

## PACS Package Exports

Installed PACS package:

```text
/home/woodson/R/x86_64-pc-linux-gnu-library/4.6/PACS
Version: 0.2.2
```

Exports include PACS core functions such as:

```text
pacs_test_sparse
pacs_test_logit
pacs_test_cumu
CCT_internal_horizontal
estimate_parameters
estimate_parameters_null
compare_models
```

Exports do not include:

```text
seurat_method2_subsample
snapATAC_method
fisher_method
archR_method
```

## Conclusion

The original author baseline file was not found locally, in available git history, or in the installed PACS package. Because it is missing, this project currently cannot claim a complete author-baseline reproduction.

The project now uses:

```text
/home/woodson/PACS_reproducing/baseline_methods_notebook1.R
```

when `--run_baselines TRUE` and no original `other_methods_for_differential_updated.R` exists in the project directory. These are clean-room reimplemented/approximate baselines and must be reported as such.
