# P56 PACS Batch-Filtering UMAP Replot Note

- date_time: 2026-06-04 21:07:34
- before_embedding_file: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/p56_before_lsi_umap_embedding.csv`
- after_embedding_file: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/p56_after_lsi_umap_embedding.csv`
- before_coordinate_columns: before_umap_1, before_umap_2
- after_coordinate_columns: after_umap_1, after_umap_2
- point_size: 0.6
- point_alpha: 0.95
- legend_point_size: 3

## Output figures overwritten

- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_before_by_batch.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_before_by_batch.pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_before_by_celltype.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_before_by_celltype.pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_after_by_batch.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_after_by_batch.pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_after_by_celltype.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_after_by_celltype.pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_four_panel.png`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_four_panel.pdf`

## Analysis note

This replot did not rerun PACS, LSI, UMAP, or matrix streaming.

After filtering, partial batch-associated structure is reduced, but residual batch separation remains. This is an initial P56-only top-peak PACS filtering result, not a final full correction.

Combined figure status: combined figure overwritten with patchwork
