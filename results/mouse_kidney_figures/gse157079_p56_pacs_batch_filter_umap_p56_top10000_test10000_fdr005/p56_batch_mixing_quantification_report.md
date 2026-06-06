# P56 Batch Mixing Quantification Report

This report quantifies existing before/after P56 PACS batch-filtering UMAP embeddings. It did not rerun PACS, MatrixMarket streaming, TF-IDF, LSI, or UMAP.

## Inputs

- before embedding: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/p56_before_lsi_umap_embedding.csv`
- after embedding: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/p56_after_lsi_umap_embedding.csv`
- before coordinate columns: before_umap_1, before_umap_2
- after coordinate columns: after_umap_1, after_umap_2
- k for kNN metrics: 30
- max_silhouette_n: 5000

## Global Metrics

```text
   stage batch_silhouette cell_type_silhouette normalized_batch_entropy
1 before        0.2024289            0.1857221                0.3041490
2  after        0.0132831            0.3019770                0.7111057
  same_batch_fraction same_celltype_fraction batch_prediction_accuracy
1           0.9870127              0.8446276                 0.7767373
2           0.6810094              0.7759180                 0.5754066
```

Metric interpretation:

- Lower batch silhouette indicates better batch mixing.
- Higher normalized batch entropy indicates better batch mixing.
- Lower same-batch fraction indicates better batch mixing.
- Higher same-celltype fraction and cell type silhouette suggest preserved biological structure.
- Lower batch prediction accuracy suggests weaker batch separability in UMAP space.

## Cell-Type-Stratified Metrics

```text
   cell_type n_cells                     batch_table before_batch_entropy
1        CNT     560   P56_batch1=294;P56_batch2=266            0.4290126
2        DCT     741   P56_batch1=307;P56_batch2=434            0.4735950
3       Endo    1143   P56_batch1=712;P56_batch2=431            0.5125635
4         IC     409   P56_batch1=221;P56_batch2=188            0.5824589
5     immune     185    P56_batch1=84;P56_batch2=101            0.5764773
6        LOH    1947  P56_batch1=1042;P56_batch2=905            0.3482285
7         NP     143     P56_batch1=90;P56_batch2=53            0.4674787
8     NP_LOH     161     P56_batch1=92;P56_batch2=69            0.4268670
9         PC     478   P56_batch1=266;P56_batch2=212            0.4810197
10      Podo     196    P56_batch1=157;P56_batch2=39            0.7286398
11        PT    4757 P56_batch1=2459;P56_batch2=2298            0.2787261
12    PT_out     237   P56_batch1=106;P56_batch2=131            0.8012531
13       PT2    1482   P56_batch1=657;P56_batch2=825            0.2415009
14   stroma1     208   P56_batch1=107;P56_batch2=101            0.7351047
15   stroma2     879   P56_batch1=535;P56_batch2=344            0.3927934
   after_batch_entropy before_same_batch_fraction after_same_batch_fraction
1            0.7360872                  0.9847619                 0.6605952
2            0.6748091                  0.9735942                 0.6898785
3            0.6126541                  0.9704287                 0.7586177
4            0.8300575                  0.9462918                 0.5958435
5            0.9469118                  0.8349550                 0.5371171
6            0.6903673                  0.9792672                 0.6923472
7            0.8060068                  0.9030303                 0.5939394
8            0.8429089                  0.9848861                 0.6016563
9            0.8031320                  0.9605300                 0.6231520
10           0.6246688                  0.9076531                 0.7171769
11           0.7281513                  0.9912690                 0.6692944
12           0.7704684                  0.8789030                 0.6390999
13           0.6759204                  0.9889564                 0.7149123
14           0.8600370                  0.8472756                 0.5788462
15           0.7015127                  0.9624194                 0.6831248
   before_batch_silhouette after_batch_silhouette
1               0.54918855             0.10719544
2               0.45815818             0.05546438
3               0.56797606             0.14713199
4               0.47168195             0.05290331
5               0.08152516             0.02670932
6               0.66160256             0.13940081
7               0.20444064             0.09872760
8               0.64976686             0.03160949
9               0.28926678             0.08967286
10              0.49993298             0.34868639
11              0.80502459             0.08280998
12              0.32896454             0.05413021
13              0.89254104             0.15229375
14              0.14719270             0.10125634
15              0.44793999             0.15037951
```

## PACS Peak Summary

  fdr_cutoff tested_peaks significant_batch_peaks retained_tested_peaks
1       0.05        10000                    6305                  3695

## PACS p-value/FDR Summary

```text
    p_value               fdr           
 Min.   :0.000e+00   Min.   :0.000e+00  
 1st Qu.:1.100e-07   1st Qu.:4.400e-07  
 Median :2.446e-03   Median :4.891e-03  
 Mean   :1.501e-01   Mean   :1.699e-01  
 3rd Qu.:1.661e-01   3rd Qu.:2.214e-01  
 Max.   :1.000e+00   Max.   :1.000e+00  
```

## Top 20 Most Significant Batch Peaks

```text
     peak_index seqnames     start       end width strand
429      294625     chrM     16101     16392   292      *
8261     147123    chr19  36918398  36919556  1159      *
2166     295348     chrX  13280132  13282295  2164      *
3693     115842    chr16  55933902  55935198  1297      *
9075      77294    chr13  65278230  65279549  1320      *
6692      65242    chr12  91384082  91384894   813      *
4042      89096    chr14  55671725  55672623   899      *
1184     110183    chr16  14159040  14160088  1049      *
6772     248696     chr7  44373723  44374830  1108      *
1884     139642    chr18  67448705  67449798  1094      *
1774     248317     chr7  38084780  38085566   787      *
6632     298309     chrX 103621991 103622879   889      *
7200     248593     chr7  43671576  43672353   778      *
1863     147210    chr19  37375754  37376692   939      *
6525     293972     chr9 120933302 120934064   763      *
5090      72389    chr13  35905710  35906739  1030      *
9567     246990     chr7  27731210  27732049   840      *
5928     119829    chr16  98081346  98083016  1671      *
1910     283170     chr9  48984813  48985852  1040      *
7157     180539     chr3  75956215  75957491  1277      *
                         name peak_detection top_peak_col       p_value
429          chrM:16101-16392           5225         9865 2.562207e-229
8261  chr19:36918398-36919556           1888         4853 1.761281e-117
2166   chrX:13280132-13282295           3537         9890 1.505370e-114
3693  chr16:55933902-55935198           2910         3781 2.438053e-101
9075  chr13:65278230-65279549           1772         2573  3.054652e-87
6692  chr12:91384082-91384894           2154         2237  3.401400e-86
4042  chr14:55671725-55672623           2794         2930  1.386842e-84
1184  chr16:14159040-14160088           4187         3620  7.806162e-82
6772   chr7:44373723-44374830           2138         8268  1.724548e-80
1884  chr18:67448705-67449798           3694         4542  8.459872e-78
1774   chr7:38084780-38085566           3759         8243  5.215331e-74
6632 chrX:103621991-103622879           2168         9957  3.938614e-70
7200   chr7:43671576-43672353           2060         8263  4.626810e-67
1863  chr19:37375754-37376692           3708         4859  1.726026e-66
6525 chr9:120933302-120934064           2189         9834  2.516312e-66
5090  chr13:35905710-35906739           2505         2422  9.558725e-65
9567   chr7:27731210-27732049           1709         8169  4.995257e-64
5928  chr16:98081346-98083016           2318         3871  5.760419e-64
1910   chr9:48984813-48985852           3676         9457  1.155847e-63
7157   chr3:75956215-75957491           2069         5867  6.783176e-63
               fdr is_batch_peak
429  2.562207e-225          TRUE
8261 8.806407e-114          TRUE
2166 5.017900e-111          TRUE
3693  6.095132e-98          TRUE
9075  6.109303e-84          TRUE
6692  5.669000e-83          TRUE
4042  1.981202e-81          TRUE
1184  9.757702e-79          TRUE
6772  1.916165e-77          TRUE
1884  8.459872e-75          TRUE
1774  4.741210e-71          TRUE
6632  3.282178e-67          TRUE
7200  3.559085e-64          TRUE
1863  1.232876e-63          TRUE
6525  1.677542e-63          TRUE
5090  5.974203e-62          TRUE
9567  2.938386e-61          TRUE
5928  3.200233e-61          TRUE
1910  6.083405e-61          TRUE
7157  3.391588e-60          TRUE
```

## Output Files

- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/p56_batch_mixing_metrics_summary.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/p56_batch_mixing_by_celltype.csv`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_batch_mixing_p56_top10000_test10000_fdr005_metrics_barplot.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_batch_mixing_p56_top10000_test10000_fdr005_by_celltype_entropy.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_batch_mixing_p56_top10000_test10000_fdr005_by_celltype_samebatch.png/pdf`

## Interpretation

Batch silhouette decreased, normalized batch entropy increased, and same-batch fraction decreased. These metrics support improved batch mixing after PACS filtering.
Cell type structure appears reasonably preserved by the same-celltype and cell type silhouette metrics.
Residual batch separation may remain; this quantification should be interpreted as an initial P56-only top-peak PACS filtering assessment, not a final full correction.
