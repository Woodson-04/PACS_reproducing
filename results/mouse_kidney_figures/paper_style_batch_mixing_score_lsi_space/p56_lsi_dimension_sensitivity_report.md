# P56 LSI dimension sensitivity for paper-style batch mixing score

## Scores

```text
  setting dim_label coordinate_space_type n_cells  k n_batches n_cell_types
1  before  LSI_1:30             LSI-space   13526 30         2           15
2   after  LSI_1:30             LSI-space   13526 30         2           15
3  before  LSI_2:30             LSI-space   13526 30         2           15
4   after  LSI_2:30             LSI-space   13526 30         2           15
5  before  LSI_2:50             LSI-space   13526 30         2           15
6   after  LSI_2:50             LSI-space   13526 30         2           15
7  before  LSI_1:50             LSI-space   13526 30         2           15
8   after  LSI_1:50             LSI-space   13526 30         2           15
  observed_batch_mixing_score expected_batch_mixing_score
1                  0.02822958                    0.490207
2                  0.13738233                    0.490207
3                  0.02549411                    0.490207
4                  0.16389423                    0.490207
5                  0.03885110                    0.490207
6                  0.18283553                    0.490207
7                  0.03754005                    0.490207
8                  0.15135788                    0.490207
  normalized_batch_mixing_score
1                    0.05758707
2                    0.28025370
3                    0.05200683
4                    0.33433678
5                    0.07925448
6                    0.37297617
7                    0.07657999
8                    0.30876320
                                                                                                                                                                                                                                                                                                                               coordinate_columns_used
1                                                                                                                                             LSI_1;LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30
2                                                                                                                                             LSI_1;LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30
3                                                                                                                                                   LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30
4                                                                                                                                                   LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30
5       LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30;LSI_31;LSI_32;LSI_33;LSI_34;LSI_35;LSI_36;LSI_37;LSI_38;LSI_39;LSI_40;LSI_41;LSI_42;LSI_43;LSI_44;LSI_45;LSI_46;LSI_47;LSI_48;LSI_49;LSI_50
6       LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30;LSI_31;LSI_32;LSI_33;LSI_34;LSI_35;LSI_36;LSI_37;LSI_38;LSI_39;LSI_40;LSI_41;LSI_42;LSI_43;LSI_44;LSI_45;LSI_46;LSI_47;LSI_48;LSI_49;LSI_50
7 LSI_1;LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30;LSI_31;LSI_32;LSI_33;LSI_34;LSI_35;LSI_36;LSI_37;LSI_38;LSI_39;LSI_40;LSI_41;LSI_42;LSI_43;LSI_44;LSI_45;LSI_46;LSI_47;LSI_48;LSI_49;LSI_50
8 LSI_1;LSI_2;LSI_3;LSI_4;LSI_5;LSI_6;LSI_7;LSI_8;LSI_9;LSI_10;LSI_11;LSI_12;LSI_13;LSI_14;LSI_15;LSI_16;LSI_17;LSI_18;LSI_19;LSI_20;LSI_21;LSI_22;LSI_23;LSI_24;LSI_25;LSI_26;LSI_27;LSI_28;LSI_29;LSI_30;LSI_31;LSI_32;LSI_33;LSI_34;LSI_35;LSI_36;LSI_37;LSI_38;LSI_39;LSI_40;LSI_41;LSI_42;LSI_43;LSI_44;LSI_45;LSI_46;LSI_47;LSI_48;LSI_49;LSI_50
```

## Before-to-after improvement

```text
  dim_label     before     after improvement
1  LSI_1:30 0.05758707 0.2802537   0.2226666
3  LSI_2:30 0.05200683 0.3343368   0.2823300
5  LSI_2:50 0.07925448 0.3729762   0.2937217
7  LSI_1:50 0.07657999 0.3087632   0.2321832
```

## Interpretation

- Main LSI metric remains LSI_2:30.
- The purpose of this sensitivity check is to confirm that the before-to-after improvement is robust to including/excluding LSI_1 and using 30 vs 50 dimensions.
- If all dimension sets improve from before to after, excluding LSI_1 does not change the main conclusion.

## LSI-depth correlation summary

- before strongest Spearman: LSI_1 = -0.9983
- after strongest Spearman: LSI_1 = 0.9976
- LSI_1 correlations: before Spearman=-0.9983, Pearson=-0.9739; after Spearman=0.9976, Pearson=0.983

## Output files

- `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_dimension_sensitivity_scores.csv`
- `figures/mouse_kidney/gse157079_p56_lsi_dimension_sensitivity_batch_mixing.png/pdf`
