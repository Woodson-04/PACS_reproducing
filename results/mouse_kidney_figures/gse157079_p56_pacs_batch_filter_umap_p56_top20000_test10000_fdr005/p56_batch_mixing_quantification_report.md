# P56 Batch Mixing Quantification Report

This report quantifies existing before/after P56 PACS batch-filtering UMAP embeddings. It did not rerun PACS, MatrixMarket streaming, TF-IDF, LSI, or UMAP.

## Inputs

- before embedding: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005/p56_before_lsi_umap_embedding.csv`
- after embedding: `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005/p56_after_lsi_umap_embedding.csv`
- before coordinate columns: before_umap_1, before_umap_2
- after coordinate columns: after_umap_1, after_umap_2
- k for kNN metrics: 30
- max_silhouette_n: 5000

## Global Metrics

```text
   stage batch_silhouette cell_type_silhouette normalized_batch_entropy
1 before       0.18373186            0.1889521                0.3101476
2  after       0.01394855            0.3025670                0.6878443
  same_batch_fraction same_celltype_fraction batch_prediction_accuracy
1           0.9903026              0.8871457                 0.7789552
2           0.7015674              0.7802923                 0.5320355
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
1        CNT     560   P56_batch1=294;P56_batch2=266            0.6353905
2        DCT     741   P56_batch1=307;P56_batch2=434            0.4159729
3       Endo    1143   P56_batch1=712;P56_batch2=431            0.3797377
4         IC     409   P56_batch1=221;P56_batch2=188            0.4477489
5     immune     185    P56_batch1=84;P56_batch2=101            0.5200823
6        LOH    1947  P56_batch1=1042;P56_batch2=905            0.2799344
7         NP     143     P56_batch1=90;P56_batch2=53            0.7852707
8     NP_LOH     161     P56_batch1=92;P56_batch2=69            0.4002844
9         PC     478   P56_batch1=266;P56_batch2=212            0.5349611
10      Podo     196    P56_batch1=157;P56_batch2=39            0.7656854
11        PT    4757 P56_batch1=2459;P56_batch2=2298            0.2486410
12    PT_out     237   P56_batch1=106;P56_batch2=131            0.6089727
13       PT2    1482   P56_batch1=657;P56_batch2=825            0.4095333
14   stroma1     208   P56_batch1=107;P56_batch2=101            0.6076812
15   stroma2     879   P56_batch1=535;P56_batch2=344            0.2899440
   after_batch_entropy before_same_batch_fraction after_same_batch_fraction
1            0.7182088                  0.9747024                 0.6780952
2            0.6537748                  0.9758884                 0.7200180
3            0.5903406                  0.9738991                 0.7994459
4            0.8318855                  0.9679707                 0.6079055
5            0.9222554                  0.8691892                 0.5592793
6            0.6898384                  0.9871597                 0.7071563
7            0.8029198                  0.7934732                 0.5976690
8            0.8321242                  0.9571429                 0.6198758
9            0.8251825                  0.9521618                 0.6041144
10           0.6199181                  0.9044218                 0.7102041
11           0.7053940                  0.9960970                 0.6814099
12           0.7394704                  0.9213783                 0.6454290
13           0.6274168                  0.9966037                 0.7419928
14           0.8472389                  0.8633013                 0.5850962
15           0.6303373                  0.9657565                 0.7706485
   before_batch_silhouette after_batch_silhouette
1               0.50864224             0.11625856
2               0.51484747             0.01925881
3               0.53691415             0.19296854
4               0.52410421             0.05361726
5               0.09518671             0.03300957
6               0.66873977             0.15486653
7               0.16294033             0.08496871
8               0.60780961             0.01463821
9               0.29869849             0.08693238
10              0.47397814             0.39674865
11              0.81682749             0.06646684
12              0.35535773             0.09027369
13              0.90017752             0.12820607
14              0.14948480             0.09729576
15              0.50193398             0.17177314
```

## PACS Peak Summary

  fdr_cutoff tested_peaks significant_batch_peaks retained_tested_peaks
1       0.05        10000                    6208                  3792

## PACS p-value/FDR Summary

```text
    p_value               fdr           
 Min.   :0.000e+00   Min.   :0.000e+00  
 1st Qu.:1.500e-07   1st Qu.:6.100e-07  
 Median :2.549e-03   Median :5.097e-03  
 Mean   :1.497e-01   Mean   :1.701e-01  
 3rd Qu.:1.766e-01   3rd Qu.:2.354e-01  
 Max.   :1.000e+00   Max.   :1.000e+00  
```

## Top 20 Most Significant Batch Peaks

```text
     peak_index seqnames     start       end width strand
429      294625     chrM     16101     16392   292      *
8261     147123    chr19  36918398  36919556  1159      *
1910     283170     chr9  48984813  48985852  1040      *
3693     115842    chr16  55933902  55935198  1297      *
2285      51839    chr11 102406999 102408205  1207      *
2166     295348     chrX  13280132  13282295  2164      *
9075      77294    chr13  65278230  65279549  1320      *
6632     298309     chrX 103621991 103622879   889      *
6692      65242    chr12  91384082  91384894   813      *
1184     110183    chr16  14159040  14160088  1049      *
6772     248696     chr7  44373723  44374830  1108      *
1884     139642    chr18  67448705  67449798  1094      *
4042      89096    chr14  55671725  55672623   899      *
1774     248317     chr7  38084780  38085566   787      *
6525     293972     chr9 120933302 120934064   763      *
3108     246943     chr7  27481557  27482096   540      *
1405      55216    chr11 119040728 119041562   835      *
1410     182418     chr3  90051810  90053343  1534      *
2423      37568    chr11   5444293   5445534  1242      *
3423      90603    chr14  65149466  65150464   999      *
                          name peak_detection top_peak_col       p_value
429           chrM:16101-16392           5225        19669 1.977001e-212
8261   chr19:36918398-36919556           1888         9718 1.192208e-105
1910    chr9:48984813-48985852           3676        18876  4.754754e-90
3693   chr16:55933902-55935198           2910         7624  6.737720e-82
2285 chr11:102406999-102408205           3474         3513  6.277043e-78
2166    chrX:13280132-13282295           3537        19720  8.695702e-78
9075   chr13:65278230-65279549           1772         5140  1.109581e-77
6632  chrX:103621991-103622879           2168        19895  2.572173e-77
6692   chr12:91384082-91384894           2154         4394  3.558559e-75
1184   chr16:14159040-14160088           4187         7274  1.202811e-72
6772    chr7:44373723-44374830           2138        16592  4.559411e-72
1884   chr18:67448705-67449798           3694         9106  7.571166e-72
4042   chr14:55671725-55672623           2794         5855  1.001692e-70
1774    chr7:38084780-38085566           3759        16559  1.033600e-68
6525  chr9:120933302-120934064           2189        19612  1.281699e-68
3108    chr7:27481557-27482096           3114        16411  1.537719e-65
1405 chr11:119040728-119041562           3998         3788  5.851741e-64
1410    chr3:90051810-90053343           3996        11954  1.239321e-61
2423     chr11:5444293-5445534           3407         2282  1.343809e-61
3423   chr14:65149466-65150464           3003         5950  4.740117e-58
               fdr is_batch_peak
429  1.977001e-208          TRUE
8261 5.961040e-102          TRUE
1910  1.584918e-86          TRUE
3693  1.684430e-78          TRUE
2285  1.255409e-74          TRUE
2166  1.449284e-74          TRUE
9075  1.585116e-74          TRUE
6632  3.215216e-74          TRUE
6692  3.953955e-72          TRUE
1184  1.202811e-69          TRUE
6772  4.144919e-69          TRUE
1884  6.309305e-69          TRUE
4042  7.705321e-68          TRUE
1774  7.382861e-66          TRUE
6525  8.544660e-66          TRUE
3108  9.610743e-63          TRUE
1405  3.442200e-61          TRUE
1410  6.885116e-59          TRUE
2423  7.072677e-59          TRUE
3423  2.370059e-55          TRUE
```

## Output Files

- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005/p56_batch_mixing_metrics_summary.csv`
- `/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005/p56_batch_mixing_by_celltype.csv`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_batch_mixing_p56_top20000_test10000_fdr005_metrics_barplot.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_batch_mixing_p56_top20000_test10000_fdr005_by_celltype_entropy.png/pdf`
- `/home/woodson/PACS_reproducing/figures/mouse_kidney/gse157079_p56_batch_mixing_p56_top20000_test10000_fdr005_by_celltype_samebatch.png/pdf`

## Interpretation

Batch silhouette decreased, normalized batch entropy increased, and same-batch fraction decreased. These metrics support improved batch mixing after PACS filtering.
Cell type structure metrics dropped noticeably. This may indicate possible overcorrection or loss of biological structure and should be checked visually.
Residual batch separation may remain; this quantification should be interpreted as an initial P56-only top-peak PACS filtering assessment, not a final full correction.
