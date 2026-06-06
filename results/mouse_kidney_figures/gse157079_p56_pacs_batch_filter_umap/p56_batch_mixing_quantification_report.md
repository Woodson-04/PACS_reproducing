# P56 Batch Mixing Quantification Report

This report quantifies existing before/after P56 PACS batch-filtering UMAP embeddings. It did not rerun PACS, MatrixMarket streaming, TF-IDF, LSI, or UMAP.

## Inputs

- before embedding: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap/p56_before_lsi_umap_embedding.csv`
- after embedding: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap/p56_after_lsi_umap_embedding.csv`
- before coordinate columns: before_umap_1, before_umap_2
- after coordinate columns: after_umap_1, after_umap_2
- k for kNN metrics: 30
- max_silhouette_n: 5000

## Global Metrics

```text
   stage batch_silhouette cell_type_silhouette normalized_batch_entropy
1 before       0.18305278            0.1755883                0.3036267
2  after       0.02661247            0.1593282                0.6956380
  same_batch_fraction same_celltype_fraction batch_prediction_accuracy
1           0.9865814              0.8488516                 0.7698374
2           0.6969269              0.6650426                 0.5899458
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
1        CNT     560   P56_batch1=294;P56_batch2=266            0.4530652
2        DCT     741   P56_batch1=307;P56_batch2=434            0.4880302
3       Endo    1143   P56_batch1=712;P56_batch2=431            0.4895370
4         IC     409   P56_batch1=221;P56_batch2=188            0.6298922
5     immune     185    P56_batch1=84;P56_batch2=101            0.5815106
6        LOH    1947  P56_batch1=1042;P56_batch2=905            0.3316718
7         NP     143     P56_batch1=90;P56_batch2=53            0.5138885
8     NP_LOH     161     P56_batch1=92;P56_batch2=69            0.4192464
9         PC     478   P56_batch1=266;P56_batch2=212            0.4725820
10      Podo     196    P56_batch1=157;P56_batch2=39            0.7447062
11        PT    4757 P56_batch1=2459;P56_batch2=2298            0.2695614
12    PT_out     237   P56_batch1=106;P56_batch2=131            0.7636267
13       PT2    1482   P56_batch1=657;P56_batch2=825            0.2574985
14   stroma1     208   P56_batch1=107;P56_batch2=101            0.7085141
15   stroma2     879   P56_batch1=535;P56_batch2=344            0.3977956
   after_batch_entropy before_same_batch_fraction after_same_batch_fraction
1            0.6782561                  0.9852381                 0.6847619
2            0.6502709                  0.9757085                 0.7082771
3            0.7255323                  0.9701954                 0.6626422
4            0.8034210                  0.9419723                 0.6192339
5            0.8937424                  0.8320721                 0.5657658
6            0.6872658                  0.9788735                 0.7053929
7            0.8440322                  0.9114219                 0.5986014
8            0.8004100                  0.9805383                 0.6360248
9            0.7252382                  0.9543236                 0.6759414
10           0.6054527                  0.9040816                 0.7641156
11           0.6593874                  0.9908906                 0.7197463
12           0.8050384                  0.8853727                 0.6992968
13           0.6460881                  0.9899010                 0.7393612
14           0.8524782                  0.8500000                 0.5860577
15           0.8179107                  0.9591581                 0.6099735
   before_batch_silhouette after_batch_silhouette
1               0.56926228             0.17890772
2               0.46357585             0.08293761
3               0.55983649             0.10766190
4               0.32611226             0.06632081
5               0.08179821             0.02613881
6               0.66419088             0.18711808
7               0.24135418             0.11615647
8               0.66746238             0.07778937
9               0.34344859             0.10730115
10              0.47060424             0.40031907
11              0.82063957             0.17172361
12              0.32676130             0.08094852
13              0.89812318             0.18324958
14              0.14884464             0.07057770
15              0.44989390             0.12202106
```

## PACS Peak Summary

  fdr_cutoff tested_peaks significant_batch_peaks retained_tested_peaks
1       0.05         5000                    3140                  1860

## PACS p-value/FDR Summary

```text
    p_value               fdr           
 Min.   :0.000e+00   Min.   :0.000e+00  
 1st Qu.:4.000e-08   1st Qu.:1.600e-07  
 Median :2.657e-03   Median :5.313e-03  
 Mean   :1.583e-01   Mean   :1.786e-01  
 3rd Qu.:1.747e-01   3rd Qu.:2.328e-01  
 Max.   :1.000e+00   Max.   :1.000e+00  
```

## Top 20 Most Significant Batch Peaks

```text
     peak_index seqnames     start       end width strand
429      294625     chrM     16101     16392   292      *
2166     295348     chrX  13280132  13282295  2164      *
3693     115842    chr16  55933902  55935198  1297      *
4042      89096    chr14  55671725  55672623   899      *
1184     110183    chr16  14159040  14160088  1049      *
1884     139642    chr18  67448705  67449798  1094      *
1774     248317     chr7  38084780  38085566   787      *
1863     147210    chr19  37375754  37376692   939      *
1910     283170     chr9  48984813  48985852  1040      *
2423      37568    chr11   5444293   5445534  1242      *
4828     197803     chr4  83417255  83418227   973      *
2690     236676     chr6  87671904  87672681   778      *
1410     182418     chr3  90051810  90053343  1534      *
3423      90603    chr14  65149466  65150464   999      *
2438     136057    chr18  42261916  42262659   744      *
1510     173506     chr2 173659163 173660161   999      *
3108     246943     chr7  27481557  27482096   540      *
1798     125904    chr17  46383314  46384081   768      *
3471     294795     chrX   7762493   7763042   550      *
1019     295510     chrX  18161864  18163542  1679      *
                         name peak_detection top_peak_col       p_value
429          chrM:16101-16392           5225         9865 2.562207e-229
2166   chrX:13280132-13282295           3537         9890 1.505370e-114
3693  chr16:55933902-55935198           2910         3781 2.438053e-101
4042  chr14:55671725-55672623           2794         2930  1.386842e-84
1184  chr16:14159040-14160088           4187         3620  7.806162e-82
1884  chr18:67448705-67449798           3694         4542  8.459872e-78
1774   chr7:38084780-38085566           3759         8243  5.215331e-74
1863  chr19:37375754-37376692           3708         4859  1.726026e-66
1910   chr9:48984813-48985852           3676         9457  1.155847e-63
2423    chr11:5444293-5445534           3407         1167  2.463066e-62
4828   chr4:83417255-83418227           2573         6388  6.974704e-61
2690   chr6:87671904-87672681           3286         7783  4.718266e-60
1410   chr3:90051810-90053343           3996         5930  6.689200e-60
3423  chr14:65149466-65150464           3003         2977  7.149969e-57
2438  chr18:42261916-42262659           3402         4474  1.164592e-56
1510 chr2:173659163-173660161           3928         5708  1.830730e-56
3108   chr7:27481557-27482096           3114         8161  5.235876e-56
1798  chr17:46383314-46384081           3748         4168  3.903916e-55
3471     chrX:7762493-7763042           2987         9875  7.353489e-55
1019   chrX:18161864-18163542           4350         9894  9.275636e-55
               fdr is_batch_peak
429  1.281104e-225          TRUE
2166 3.763425e-111          TRUE
3693  4.063421e-98          TRUE
4042  1.733552e-81          TRUE
1184  7.806162e-79          TRUE
1884  7.049893e-75          TRUE
1774  3.725236e-71          TRUE
1863  1.078766e-63          TRUE
1910  6.421372e-61          TRUE
2423  1.231533e-59          TRUE
4828  3.170320e-58          TRUE
2690  1.965944e-57          TRUE
1410  2.572769e-57          TRUE
3423  2.553560e-54          TRUE
2438  3.881974e-54          TRUE
1510  5.721031e-54          TRUE
3108  1.539964e-53          TRUE
1798  1.084421e-52          TRUE
3471  1.935129e-52          TRUE
1019  2.318909e-52          TRUE
```

## Output Files

- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap/p56_batch_mixing_metrics_summary.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap/p56_batch_mixing_by_celltype.csv`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_batch_mixing_metrics_barplot.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_batch_mixing_by_celltype_entropy.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_batch_mixing_by_celltype_samebatch.png/pdf`

## Interpretation

Batch silhouette decreased, normalized batch entropy increased, and same-batch fraction decreased. These metrics support improved batch mixing after PACS filtering.
Cell type structure metrics dropped noticeably. This may indicate possible overcorrection or loss of biological structure and should be checked visually.
Residual batch separation may remain; this quantification should be interpreted as an initial P56-only top-peak PACS filtering assessment, not a final full correction.
