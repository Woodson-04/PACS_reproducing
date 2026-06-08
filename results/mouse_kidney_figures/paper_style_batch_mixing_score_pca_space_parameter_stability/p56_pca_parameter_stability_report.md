# P56 PCA-logNorm Parameter Stability Report

## Input Audit Summary

```text
 setting_label
     5000/5000
   10000/10000
   20000/10000
                                                                                                                               result_dir
                                         /home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap
 /home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved
           /home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005
 n_top_peaks tested_peaks fdr_cutoff significant_batch_peaks retained_peaks
       10000         5000       0.05                    3140           1860
       10000        10000       0.05                    6305           3695
       20000        10000       0.05                    6208           3792
 counts_rds_exists retained_indices_exists metadata_exists
              TRUE                    TRUE            TRUE
              TRUE                    TRUE            TRUE
              TRUE                    TRUE            TRUE
 retained_has_top_peak_col can_compute_directly
                      TRUE                 TRUE
                      TRUE                 TRUE
                      TRUE                 TRUE
```

## Computed Settings

-  5000/5000
-  10000/10000
-  20000/10000

## Score Table

```text
 setting_label  stage pc_dims n_cells n_features normalized_batch_mixing_score
     5000/5000 before      20   13526      10000                    0.03976549
     5000/5000 before      30   13526      10000                    0.05103657
     5000/5000 before      50   13526      10000                    0.06850624
     5000/5000  after      20   13526       1860                    0.30482687
     5000/5000  after      30   13526       1860                    0.36788860
     5000/5000  after      50   13526       1860                    0.50315662
   10000/10000 before      20   13526      10000                    0.03976549
   10000/10000 before      30   13526      10000                    0.05103657
   10000/10000 before      50   13526      10000                    0.06850624
   10000/10000  after      20   13526       3695                    0.25760095
   10000/10000  after      30   13526       3695                    0.30972340
   10000/10000  after      50   13526       3695                    0.42389711
   20000/10000 before      20   13526      20000                    0.02948981
   20000/10000 before      30   13526      20000                    0.03982079
   20000/10000 before      50   13526      20000                    0.05561639
   20000/10000  after      20   13526       3792                    0.23908561
   20000/10000  after      30   13526       3792                    0.29280673
   20000/10000  after      50   13526       3792                    0.40968510
```

## PC1:30 Comparison

```text
 setting_label  stage normalized_batch_mixing_score
     5000/5000 before                    0.05103657
     5000/5000  after                    0.36788860
   10000/10000 before                    0.05103657
   10000/10000  after                    0.30972340
   20000/10000 before                    0.03982079
   20000/10000  after                    0.29280673
```

## PC1:50 Comparison

```text
 setting_label  stage normalized_batch_mixing_score
     5000/5000 before                    0.06850624
     5000/5000  after                    0.50315662
   10000/10000 before                    0.06850624
   10000/10000  after                    0.42389711
   20000/10000 before                    0.05561639
   20000/10000  after                    0.40968510
```

## Interpretation

- PCA-logNorm parameter stability is more consistent with the PACS paper-style normalized PCA mixing metric than UMAP-space auxiliary metrics.
- If before-to-after improvement holds across settings, this supports PACS filtering beyond visualization-space UMAP.
- The preferred setting should balance batch mixing improvement and biological structure preservation.
- Compare this PCA-space analysis with the existing UMAP/LSI summaries before updating Section 10.

## Output Files

- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space_parameter_stability/p56_pca_parameter_stability_scores.csv`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pca_space_parameter_stability_scores.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pca_space_parameter_stability_scores.pdf`
