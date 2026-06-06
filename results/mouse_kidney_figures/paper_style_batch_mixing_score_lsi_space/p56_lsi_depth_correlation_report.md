# P56 LSI-depth correlation report

## Goal

This report checks whether LSI_1 is depth-associated in the current P56 PACS main result.

## Input files

- before: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_before_lsi_embedding.csv`
- after: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_after_lsi_embedding.csv`

## Strongest absolute correlations

- before strongest Spearman: LSI_1 = -0.9983
- after strongest Spearman: LSI_1 = 0.9976
- before strongest Pearson: LSI_1 = -0.9739
- after strongest Pearson: LSI_1 = 0.983

## LSI_1 correlations

- before LSI_1 Spearman: -0.9983; Pearson: -0.9739
- after LSI_1 Spearman: 0.9976; Pearson: 0.983

## Interpretation

LSI_1 is among the strongest depth-associated components. This supports the use of LSI_2:LSI_30 for batch mixing score calculation.

## Output files

- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_depth_correlation.csv`
- `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before.png/pdf`
- `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_after.png/pdf`
- `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before_after.png/pdf`
