# P56 10000/10000 LSI-space paper-style batch mixing score

## Goal

This report computes the PACS-paper-style normalized batch mixing score in LSI space for the current main P56 PACS filtering setting.

## Coordinate-space interpretation

The PACS paper calculates the normalized batch mixing score in PCA space. For this scATAC pipeline, LSI space, especially LSI dimensions 2:30, is the closest analogue. UMAP-space scores are visualization-space approximations and are included only as reference.

## LSI-space scores

```text
                         space_name coordinate_space_type n_cells  k n_batches
1 P56 10000/10000 before, LSI space             LSI-space   13526 30         2
2  P56 10000/10000 after, LSI space             LSI-space   13526 30         2
  n_cell_types observed_batch_mixing_score expected_batch_mixing_score
1           15                  0.02549411                    0.490207
2           15                  0.16389423                    0.490207
  normalized_batch_mixing_score
1                    0.05200683
2                    0.33433678
                                                                                                                                                                             coordinate_columns_used
1 LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30
2 LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30
```

## LSI vs existing UMAP-space scores

```text
                           space_name coordinate_space_type
2  P56 10000/10000 before, UMAP space            UMAP-space
3   P56 10000/10000 after, UMAP space            UMAP-space
1   P56 10000/10000 before, LSI space             LSI-space
21   P56 10000/10000 after, LSI space             LSI-space
   observed_batch_mixing_score expected_batch_mixing_score
2                   0.01298733                    0.490207
3                   0.31899059                    0.490207
1                   0.02549411                    0.490207
21                  0.16389423                    0.490207
   normalized_batch_mixing_score
2                     0.02649357
3                     0.65072631
1                     0.05200683
21                    0.33433678
                                                                                                                                                                              coordinate_columns_used
2                                                                                                                                                                         before_umap_1;before_umap_2
3                                                                                                                                                                           after_umap_1;after_umap_2
1  LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30
21 LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30
```

## Interpretation

- Higher normalized batch mixing score indicates better batch mixing.
- If the LSI-space score increases from before to after, this supports PACS filtering in a space closer to the PACS paper PCA-space calculation.
- Do not directly compare the GEO UMAP-space reference score with these LSI-space PACS scores.

## Output files

- `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_10000_lsi_space_paper_style_scores.csv`
- `figures/mouse_kidney/gse157079_p56_10000_paper_style_lsi_vs_umap_score.png/pdf`
- Per-cell CSV files were written for before and after LSI-space scores.
