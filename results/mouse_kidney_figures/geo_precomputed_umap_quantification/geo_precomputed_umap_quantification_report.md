# GEO precomputed UMAP quantitative reference

## Goal

This report evaluates the GEO-provided precomputed GSE157079 UMAP coordinates as a public reference embedding. It should not be interpreted as the PACS paper filtered UMAP unless directly documented.

## Global GSE157079 UMAP metrics

The global analysis uses sample as the sample/batch label and cell_type as the biological label. Age/developmental structure is also summarized separately because P0, P21, and P56 differences are biological rather than purely technical batch effects.

```text
                        stage batch_silhouette cell_type_silhouette
1 GEO precomputed UMAP global      -0.08034808            0.5640647
  normalized_batch_entropy same_batch_fraction same_celltype_fraction
1                0.8584297           0.3380468              0.9822315
  batch_prediction_accuracy age_group_silhouette same_age_fraction
1                  0.319482           -0.0497906         0.6205078
  normalized_age_entropy age_group_prediction_accuracy
1               0.707601                     0.6044732
```

## P56-only GEO UMAP metrics

This is the most directly comparable GEO reference for the P56_batch1 vs P56_batch2 PACS-filtered UMAP analyses.

```text
                            stage batch_silhouette cell_type_silhouette
1 GEO precomputed UMAP P56 subset      0.001408244              0.58692
  normalized_batch_entropy same_batch_fraction same_celltype_fraction
1                0.9256912           0.5396791              0.9753832
  batch_prediction_accuracy
1                 0.5468211
```

## P56 cell-type-stratified GEO metrics

```text
   cell_type n_cells                     batch_table batch_silhouette
1        CNT     560   P56_batch1=294;P56_batch2=266     0.0232929088
2        DCT     741   P56_batch1=307;P56_batch2=434    -0.0217298696
3       Endo    1143   P56_batch1=712;P56_batch2=431     0.0163954262
4         IC     409   P56_batch1=221;P56_batch2=188     0.0484640701
5     immune     185    P56_batch1=84;P56_batch2=101    -0.0033500453
6        LOH    1947  P56_batch1=1042;P56_batch2=905     0.0070420380
7         NP     143     P56_batch1=90;P56_batch2=53     0.0084047106
8     NP_LOH     161     P56_batch1=92;P56_batch2=69     0.0332522802
9         PC     478   P56_batch1=266;P56_batch2=212     0.0171049223
10      Podo     196    P56_batch1=157;P56_batch2=39     0.1722125932
11        PT    4757 P56_batch1=2459;P56_batch2=2298     0.0016895318
12    PT_out     237   P56_batch1=106;P56_batch2=131     0.0746758384
13       PT2    1482   P56_batch1=657;P56_batch2=825    -0.0002373697
14   stroma1     208   P56_batch1=107;P56_batch2=101     0.0222673698
15   stroma2     879   P56_batch1=535;P56_batch2=344    -0.0033382349
   normalized_batch_entropy same_batch_fraction
1                 0.9122842           0.5411310
2                 0.9441470           0.5254161
3                 0.9133006           0.5478857
4                 0.8898003           0.5638142
5                 0.8425051           0.6093694
6                 0.9536874           0.5172060
7                 0.8926840           0.5701632
8                 0.8897523           0.5608696
9                 0.8888814           0.5607392
10                0.7472928           0.7585034
11                0.9431363           0.5263541
12                0.9260663           0.5369902
13                0.9023097           0.5525641
14                0.9346779           0.5240385
15                0.9043295           0.5491088
```

## Comparison with P56 PACS-filtered UMAP

```text
                           setting batch_silhouette normalized_batch_entropy
1 GEO precomputed UMAP, P56 subset      0.001408244                0.9256912
2         P56 PACS 5000/5000 after      0.026612473                0.6956380
3       P56 PACS 10000/10000 after      0.013283103                0.7111057
4       P56 PACS 20000/10000 after      0.013948545                0.6878443
  same_batch_fraction batch_prediction_accuracy cell_type_silhouette
1           0.5396791                 0.5468211            0.5869200
2           0.6969269                 0.5899458            0.1593282
3           0.6810094                 0.5754066            0.3019770
4           0.7015674                 0.5320355            0.3025670
  same_celltype_fraction
1              0.9753832
2              0.6650426
3              0.7759180
4              0.7802923
```

Lower batch silhouette, lower same-batch fraction, and lower batch prediction accuracy indicate better batch mixing. Higher normalized batch entropy indicates better batch mixing. Higher cell type silhouette and same-celltype fraction indicate stronger biological cell type structure.

## Interpretation caveat

The GEO UMAP may reflect original atlas processing choices, possibly including integration or correction steps, but we should not claim it is the PACS paper batch-filtered UMAP unless directly documented.

## Use in presentation

我们将 GEO 预计算 UMAP 作为公开参考 embedding，用同一套 batch mixing 指标评估其混合程度，并与 P56 PACS-filtered UMAP 对照。

## Output files

- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/geo_precomputed_umap_quantification/geo_precomputed_umap_global_metrics.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/geo_precomputed_umap_quantification/geo_precomputed_umap_p56_metrics.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/geo_precomputed_umap_quantification/geo_precomputed_umap_p56_by_celltype_metrics.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/geo_precomputed_umap_quantification/geo_vs_p56_pacs_umap_metrics_comparison.csv`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_geo_precomputed_umap_mixing_metrics.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_geo_precomputed_umap_p56_mixing_metrics.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_geo_precomputed_umap_p56_by_celltype_entropy.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_geo_precomputed_umap_p56_by_celltype_samebatch.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png/pdf`
