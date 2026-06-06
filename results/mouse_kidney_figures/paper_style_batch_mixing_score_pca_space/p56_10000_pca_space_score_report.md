# P56 10000/10000 PCA-logNorm paper-style batch mixing score

## Goal

This report computes a PCA-logNorm-space version of the PACS-paper-style normalized batch mixing score for the current P56 main setting.

## Inputs

- result_dir: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved`
- before sparse counts: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_counts_top_peaks_sparse.rds`
- metadata: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_metadata.csv`
- retained peaks: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_retained_peak_indices.csv`

## Matrix loading and retained peak mapping

- before matrix dim: 13526 x 10000
- after matrix dim: 13526 x 3695
- retained peak columns: 3695
- mapping: After matrix was subset using p56_retained_peak_indices.csv column `top_peak_col`.

## PCA-logNorm method

- Counts were binarized.
- Cells were library-size normalized to the median cell depth.
- Values were log1p transformed.
- Sparse truncated PCA was computed with `irlba::prcomp_irlba`.

## PCA-space scores

```text
  setting coordinate_space_type pc_dims n_cells n_features
1  before           PCA-logNorm      20   13526      10000
2   after           PCA-logNorm      20   13526       3695
3  before           PCA-logNorm      30   13526      10000
4   after           PCA-logNorm      30   13526       3695
5  before           PCA-logNorm      50   13526      10000
6   after           PCA-logNorm      50   13526       3695
  observed_batch_mixing_score expected_batch_mixing_score
1                  0.01949332                    0.490207
2                  0.12627779                    0.490207
3                  0.02501848                    0.490207
4                  0.15182858                    0.490207
5                  0.03358224                    0.490207
6                  0.20779733                    0.490207
  normalized_batch_mixing_score  k
1                    0.03976549 30
2                    0.25760095 30
3                    0.05103657 30
4                    0.30972340 30
5                    0.06850624 30
6                    0.42389711 30
                                                                                                                                                                                                                                                                             coordinate_columns_used
1                                                                                                                                                                                     PC_1;PC_2;PC_3;PC_4;PC_5;PC_6;PC_7;PC_8;PC_9;PC_10;PC_11;PC_12;PC_13;PC_14;PC_15;PC_16;PC_17;PC_18;PC_19;PC_20
2                                                                                                                                                                                     PC_1;PC_2;PC_3;PC_4;PC_5;PC_6;PC_7;PC_8;PC_9;PC_10;PC_11;PC_12;PC_13;PC_14;PC_15;PC_16;PC_17;PC_18;PC_19;PC_20
3                                                                                                                         PC_1;PC_2;PC_3;PC_4;PC_5;PC_6;PC_7;PC_8;PC_9;PC_10;PC_11;PC_12;PC_13;PC_14;PC_15;PC_16;PC_17;PC_18;PC_19;PC_20;PC_21;PC_22;PC_23;PC_24;PC_25;PC_26;PC_27;PC_28;PC_29;PC_30
4                                                                                                                         PC_1;PC_2;PC_3;PC_4;PC_5;PC_6;PC_7;PC_8;PC_9;PC_10;PC_11;PC_12;PC_13;PC_14;PC_15;PC_16;PC_17;PC_18;PC_19;PC_20;PC_21;PC_22;PC_23;PC_24;PC_25;PC_26;PC_27;PC_28;PC_29;PC_30
5 PC_1;PC_2;PC_3;PC_4;PC_5;PC_6;PC_7;PC_8;PC_9;PC_10;PC_11;PC_12;PC_13;PC_14;PC_15;PC_16;PC_17;PC_18;PC_19;PC_20;PC_21;PC_22;PC_23;PC_24;PC_25;PC_26;PC_27;PC_28;PC_29;PC_30;PC_31;PC_32;PC_33;PC_34;PC_35;PC_36;PC_37;PC_38;PC_39;PC_40;PC_41;PC_42;PC_43;PC_44;PC_45;PC_46;PC_47;PC_48;PC_49;PC_50
6 PC_1;PC_2;PC_3;PC_4;PC_5;PC_6;PC_7;PC_8;PC_9;PC_10;PC_11;PC_12;PC_13;PC_14;PC_15;PC_16;PC_17;PC_18;PC_19;PC_20;PC_21;PC_22;PC_23;PC_24;PC_25;PC_26;PC_27;PC_28;PC_29;PC_30;PC_31;PC_32;PC_33;PC_34;PC_35;PC_36;PC_37;PC_38;PC_39;PC_40;PC_41;PC_42;PC_43;PC_44;PC_45;PC_46;PC_47;PC_48;PC_49;PC_50
```

## PCA vs LSI vs UMAP comparison

```text
                        setting coordinate_space_type
1                 before PC1:20           PCA-logNorm
2                  after PC1:20           PCA-logNorm
3                 before PC1:30           PCA-logNorm
4                  after PC1:30           PCA-logNorm
5                 before PC1:50           PCA-logNorm
6                  after PC1:50           PCA-logNorm
7        P56 10000/10000 before             LSI-space
8         P56 10000/10000 after             LSI-space
9  P56 10000/10000 before, UMAP            UMAP-space
10  P56 10000/10000 after, UMAP            UMAP-space
   normalized_batch_mixing_score
1                     0.03976549
2                     0.25760095
3                     0.05103657
4                     0.30972340
5                     0.06850624
6                     0.42389711
7                     0.05200683
8                     0.33433678
9                     0.02649357
10                    0.65072631
```

## Interpretation

- PCA-space score is closest to the PACS paper normalized PCA mixing metric.
- LSI-space score is a scATAC-aware analogue.
- UMAP-space score is a visualization-level approximation.
- Do not claim exact numerical comparability to the PACS paper unless dataset, batch definition, feature set, and preprocessing match.

## Output files

- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space/p56_10000_pca_space_scores.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space/p56_10000_before_pca_logNorm_embedding.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space/p56_10000_after_pca_logNorm_embedding.csv`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png/pdf`
