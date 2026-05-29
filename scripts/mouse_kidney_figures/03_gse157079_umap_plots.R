#!/usr/bin/env Rscript

# Plot precomputed GSE157079 UMAP coordinates from the lightweight merged
# metadata table. This script does not read the large cell-by-peak matrix.
# Batch correction or recomputed UMAP would require a separate pipeline using
# the count matrix and is not performed here.

suppressPackageStartupMessages({
  library(ggplot2)
})

parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) return(defaults)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- sub("^--", "", key)
    if (!name %in% names(defaults)) stop("Unknown argument: --", name)
    if (i == length(args)) stop("Missing value for argument: --", name)
    defaults[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  defaults
}

params <- parse_args(list(
  metadata_csv = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv",
  out_dir = "/home/woodson/PACS_reproducing/figures/mouse_kidney"
))

if (!file.exists(params$metadata_csv)) {
  stop("Merged metadata CSV does not exist: ", params$metadata_csv)
}
dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, showProgress = FALSE, check.names = FALSE)))
  }
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

metadata <- read_table(params$metadata_csv)

colnames_text <- function(df) paste(names(df), collapse = ", ")

pick_col <- function(df, candidates, label, required = TRUE) {
  nms <- names(df)
  for (candidate in candidates) {
    hit <- nms[tolower(nms) == tolower(candidate)]
    if (length(hit) > 0) return(hit[[1]])
  }
  if (required) {
    stop(
      "Could not identify ", label, " column. Tried: ",
      paste(candidates, collapse = ", "),
      ". Available columns: ", colnames_text(df)
    )
  }
  NULL
}

pick_umap_cols <- function(df) {
  pairs <- list(
    c("umap_1", "umap_2"),
    c("UMAP_1", "UMAP_2"),
    c("umap.1", "umap.2"),
    c("UMAP.1", "UMAP.2"),
    c("umap-1", "umap-2"),
    c("UMAP-1", "UMAP-2"),
    c("umap1", "umap2"),
    c("UMAP1", "UMAP2")
  )
  nms_lower <- tolower(names(df))
  for (pair in pairs) {
    idx <- match(tolower(pair), nms_lower)
    if (all(!is.na(idx))) return(names(df)[idx])
  }
  candidates <- names(df)[grepl("umap", names(df), ignore.case = TRUE)]
  numeric_candidates <- candidates[vapply(df[candidates], function(x) {
    all(is.finite(suppressWarnings(as.numeric(x))[!is.na(x)]))
  }, logical(1))]
  if (length(numeric_candidates) >= 2) return(numeric_candidates[1:2])
  stop(
    "Could not identify two UMAP columns. Tried standard/fallback UMAP names. ",
    "Available columns: ", colnames_text(df)
  )
}

umap_cols <- pick_umap_cols(metadata)
celltype_col <- pick_col(
  metadata,
  c("cell_type", "celltype", "CellType", "clusters", "cluster", "annotation", "Annotation", "x.sp_cluster2", "V4"),
  "cell type/annotation"
)
sample_col <- pick_col(
  metadata,
  c("sample", "samples", "Sample", "Samples", "batch", "Batch", "condition", "Condition", "orig.ident", "donor", "Donor", "V3"),
  "sample/batch"
)
barcode_col <- pick_col(
  metadata,
  c("cell_barcode", "barcodes", "barcode", "Barcode", "cell", "cell_id", "V2"),
  "cell barcode",
  required = FALSE
)

plot_df <- data.frame(
  UMAP_1 = suppressWarnings(as.numeric(metadata[[umap_cols[[1]]]])),
  UMAP_2 = suppressWarnings(as.numeric(metadata[[umap_cols[[2]]]])),
  cell_type = as.factor(metadata[[celltype_col]]),
  stringsAsFactors = FALSE
)
plot_df$sample <- as.character(metadata[[sample_col]])
if (!is.null(barcode_col)) {
  plot_df$cell_barcode <- as.character(metadata[[barcode_col]])
}
plot_df <- plot_df[is.finite(plot_df$UMAP_1) & is.finite(plot_df$UMAP_2), ]
if (nrow(plot_df) == 0) {
  stop("No finite UMAP coordinates after parsing columns: ", paste(umap_cols, collapse = ", "))
}

sample_palette_preferred <- c(
  "P0_batch1" = "#FF0000",
  "P0_batch2" = "#0000FF",
  "P21_batch1" = "#00CC00",
  "P56_batch1" = "#FFD400",
  "P56_batch2" = "#A000FF"
)
sample_colors <- c("#FF0000", "#0000FF", "#00CC00", "#FFD400", "#A000FF", "#00FFFF", "#FF7F00")

celltype_palette_preferred <- c(
  "CNT" = "#66C7FF",
  "DCT" = "#D8896A",
  "Endo" = "#C69214",
  "IC" = "#B6A000",
  "immune" = "#7FB000",
  "LOH" = "#67B83F",
  "NP" = "#00A651",
  "NP_LOH" = "#00A98F",
  "PC" = "#00A6B8",
  "Podo" = "#00A7E1",
  "PT" = "#139DDF",
  "PT_out" = "#8F9BEF",
  "PT2" = "#B879E8",
  "stroma1" = "#DD70D6",
  "stroma2" = "#F06AA8"
)
celltype_colors <- c(
  "#66C7FF", "#D8896A", "#C69214", "#B6A000", "#7FB000",
  "#67B83F", "#00A651", "#00A98F", "#00A6B8", "#00A7E1",
  "#139DDF", "#8F9BEF", "#B879E8", "#DD70D6", "#F06AA8",
  "#F4A6C8", "#9ADBC5", "#E6C36A", "#BCA7FF", "#7DD3FC"
)

build_palette <- function(levels, preferred, fallback) {
  palette <- preferred[names(preferred) %in% levels]
  missing_levels <- setdiff(levels, names(palette))
  if (length(missing_levels) > 0) {
    extra_colors <- rep(fallback, length.out = length(missing_levels))
    names(extra_colors) <- missing_levels
    palette <- c(palette, extra_colors)
  }
  palette[levels]
}

theme_umap <- function() {
  theme_classic(base_size = 16) +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(face = "bold", size = 20),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.45, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

celltype_levels <- sort(unique(as.character(plot_df$cell_type)))
celltype_palette <- build_palette(celltype_levels, celltype_palette_preferred, celltype_colors)
plot_df$cell_type <- factor(plot_df$cell_type, levels = celltype_levels)

sample_levels <- sort(unique(as.character(plot_df$sample)))
sample_palette <- build_palette(sample_levels, sample_palette_preferred, sample_colors)
plot_df$sample <- factor(plot_df$sample, levels = sample_levels)

set.seed(1)
plot_df_celltype <- plot_df[sample(seq_len(nrow(plot_df))), ]
p_celltype <- ggplot(plot_df_celltype, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.64, alpha = 0.95, stroke = 0) +
  coord_equal() +
  scale_x_continuous(expand = expansion(mult = 0.015)) +
  scale_y_continuous(expand = expansion(mult = 0.015)) +
  scale_color_manual(values = celltype_palette, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cell type", title = "GSE157079 snATAC UMAP by cell type") +
  theme_umap()

ggsave(file.path(params$out_dir, "gse157079_umap_by_celltype.png"), p_celltype, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(params$out_dir, "gse157079_umap_by_celltype.pdf"), p_celltype, width = 8.5, height = 6.5)

set.seed(1)
plot_df_sample <- plot_df[sample(seq_len(nrow(plot_df))), ]
p_sample <- ggplot(plot_df_sample, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
  geom_point(size = 0.64, alpha = 0.95, stroke = 0) +
  coord_equal() +
  scale_x_continuous(expand = expansion(mult = 0.015)) +
  scale_y_continuous(expand = expansion(mult = 0.015)) +
  scale_color_manual(values = sample_palette, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Sample", title = "GSE157079 snATAC UMAP by sample") +
  theme_umap()

ggsave(file.path(params$out_dir, "gse157079_umap_by_sample.png"), p_sample, width = 8.5, height = 6.5, dpi = 300)
ggsave(file.path(params$out_dir, "gse157079_umap_by_sample.pdf"), p_sample, width = 8.5, height = 6.5)

cat("UMAP columns used: ", paste(umap_cols, collapse = ", "), "\n", sep = "")
cat("Cell type field used: ", celltype_col, "\n", sep = "")
cat("Sample field used: ", sample_col, "\n", sep = "")
if (!is.null(barcode_col)) cat("Barcode field detected: ", barcode_col, "\n", sep = "")
cat("Cell type levels: ", paste(celltype_levels, collapse = ", "), "\n", sep = "")
cat("Sample levels: ", paste(sample_levels, collapse = ", "), "\n", sep = "")
cat("Cells plotted: ", nrow(plot_df), "\n", sep = "")
cat("Saved UMAP figures to: ", normalizePath(params$out_dir, mustWork = FALSE), "\n", sep = "")
