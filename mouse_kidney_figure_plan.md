# Mouse Kidney Figure Reproduction Plan

## Stage 1: Figures Supported by Current Data

### 1. Cell Type Composition Bar Plot

- Data input: `x.sp_cluster2` from `data_for_test_for_t1e_power.rdata`.
- R packages: `ggplot2`.
- Output: `figures/mouse_kidney/cell_type_counts_barplot.png`.
- Can do now: yes.
- Article relation: overview/sample composition panel; not a full UMAP
  replacement.
- Difficulty: low.
- Expected visual: horizontal ordered bar plot of cell counts by annotated cell
  type.

### 2. PT vs LOH Accessibility Depth/QC Distribution

- Data input: `pmats`, `x.sp_cluster2`, `kidney_features_to_keep`,
  `r_by_ct_est$q_vec_new`.
- R packages: `Matrix`, `ggplot2`.
- Output: `figures/mouse_kidney/pt_loh_depth_distribution.png`.
- Can do now: yes.
- Article relation: QC/depth context for the PT vs LOH differential analysis.
- Difficulty: low.
- Expected visual: overlaid density curves for PT and LOH total accessibility
  counts on a log scale.

### 3. PACS P-Value QQ Plot

- Data input: `p_value_permuted_label[["our"]]` from
  `results/20260526_2318_large_baseline/pacs_kidney_notebook1_result.rds`.
- R packages: `ggplot2`.
- Output: `figures/mouse_kidney/pacs_permuted_qq_plot.png`.
- Can do now: yes.
- Article relation: diagnostic view of Type I error calibration.
- Difficulty: low.
- Expected visual: observed vs expected `-log10(p)` with diagonal reference.

### 4. PACS Type I Error / Power Benchmark Bar Plot

- Data input: `results/20260526_2318_large_baseline/summary.csv`.
- R packages: `ggplot2`.
- Output: `figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png`.
- Can do now: yes.
- Article relation: Notebook 1 benchmark summary.
- Difficulty: low.
- Expected visual: two-panel bar plot for Type I error and power across
  methods.

### 5. PT vs LOH DAR Volcano Plot

- Data input: actual-label PACS p-values from the large result RDS, plus PT/LOH
  count matrices from PACS_data.
- R packages: `Matrix`, `ggplot2`.
- Output: `figures/mouse_kidney/pt_loh_pacs_dar_volcano.png`.
- Can do now: yes, with effect size computed as a local accessibility
  difference or log fold-change proxy.
- Article relation: differential accessibility result overview, analogous to a
  DAR summary panel.
- Difficulty: medium.
- Expected visual: `log2FC` or accessibility difference against `-log10(p)`,
  with significant DARs highlighted.

### 6. Top DAR Peak Heatmap

- Data input: top actual-label PACS peaks and sampled PT/LOH matrices from the
  large result RDS/PACS_data.
- R packages: `Matrix`, `ggplot2`; optionally `pheatmap` or `ComplexHeatmap`.
- Output: `figures/mouse_kidney/top_dar_peak_heatmap.png`.
- Can do now: yes, using sampled cells/features and normalized counts.
- Article relation: differential accessibility pattern display.
- Difficulty: medium.
- Expected visual: heatmap of top DAR peaks across PT and LOH cells, grouped by
  cell type.

## Stage 2: Figures Requiring Extra Annotation or External Objects

### 1. DAR Peaks Linked Genes

- Data input: DAR peak coordinates, gene annotation, peak-to-gene or nearest
  gene mapping.
- R packages: `GenomicRanges`, `ChIPseeker` or custom annotation tools.
- Output: `figures/mouse_kidney/dar_linked_genes.png`.
- Can do now: not unless peak coordinates/gene annotation are confirmed.
- Article relation: biological interpretation of DARs.
- Difficulty: medium-high.
- Expected visual: ranked linked genes or pathway-associated DAR annotations.

### 2. Motif Enrichment

- Data input: DAR genomic ranges, genome build, motif database.
- R packages: `chromVAR`, `motifmatchr`, `JASPAR2020/2022`,
  `BSgenome.Mmusculus.*`.
- Output: `figures/mouse_kidney/dar_motif_enrichment.png`.
- Can do now: not without genome build and motif resources.
- Article relation: transcription factor/regulatory program analysis.
- Difficulty: high.
- Expected visual: enriched motif dot plot or ranked motif bar plot.

### 3. Genome Browser Style Peak Tracks

- Data input: fragment files or pseudobulk coverage tracks, peak coordinates,
  gene annotation.
- R packages: `Gviz`, `rtracklayer`, `GenomicRanges`; or external IGV tracks.
- Output: `figures/mouse_kidney/browser_style_peak_tracks.png`.
- Can do now: not with the current matrix-only data.
- Article relation: locus-level visual validation.
- Difficulty: high.
- Expected visual: PT/LOH accessibility tracks around marker loci with DAR peaks
  annotated.

### 4. Gene Activity / Regulatory Program Figures

- Data input: gene activity matrix, linked gene annotations, RNA expression data
  or the original integrated analysis object.
- R packages: `Seurat`, `Signac`, `ArchR`, `ggplot2`.
- Output: `figures/mouse_kidney/gene_activity_programs.png`.
- Can do now: not with current PACS_data alone.
- Article relation: downstream regulatory biology panels.
- Difficulty: high.
- Expected visual: heatmap or dot plot of gene activity/regulatory programs
  across kidney cell types.

## Recommended Order

1. Run `01_overview_plots.R` to establish the figure pipeline and visual style.
2. Add `02_pacs_dar_volcano_heatmap.R` for PT vs LOH DAR volcano and heatmap.
3. Inspect peak names to determine whether genomic coordinates are available.
4. Only then decide whether motif, linked-gene, or genome browser panels are
   feasible without additional downloaded data.
