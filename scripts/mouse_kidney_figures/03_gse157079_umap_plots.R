#!/usr/bin/env Rscript

# Plot GSE157079 UMAP coordinates from the lightweight merged metadata table.
# This script does not read the large cell-by-peak matrix.

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
  "sample/batch",
  required = FALSE
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
if (!is.null(sample_col)) {
  plot_df$sample <- as.factor(metadata[[sample_col]])
}
if (!is.null(barcode_col)) {
  plot_df$cell_barcode <- as.character(metadata[[barcode_col]])
}
plot_df <- plot_df[is.finite(plot_df$UMAP_1) & is.finite(plot_df$UMAP_2), ]
if (nrow(plot_df) == 0) {
  stop("No finite UMAP coordinates after parsing columns: ", paste(umap_cols, collapse = ", "))
}

discrete_palette <- function(n, option = "Dark 3") {
  if ("hcl.colors" %in% ls("package:grDevices")) {
    grDevices::hcl.colors(n, palette = option)
  } else {
    grDevices::rainbow(n)
  }
}

theme_umap <- function() {
  theme_classic(base_size = 12) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.title = element_text(face = "bold")
    )
}

cell_levels <- levels(plot_df$cell_type)
cell_pal <- setNames(discrete_palette(length(cell_levels), "Dark 3"), cell_levels)
p_celltype <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.18, alpha = 0.7, stroke = 0) +
  coord_equal() +
  scale_color_manual(values = cell_pal) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cell type", title = "GSE157079 snATAC UMAP by cell type") +
  theme_umap()

ggsave(file.path(params$out_dir, "gse157079_umap_by_celltype.png"), p_celltype, width = 8, height = 6, dpi = 300)
ggsave(file.path(params$out_dir, "gse157079_umap_by_celltype.pdf"), p_celltype, width = 8, height = 6)

if (!is.null(sample_col)) {
  sample_levels <- levels(plot_df$sample)
  sample_pal <- setNames(discrete_palette(length(sample_levels), "Set 2"), sample_levels)
  p_sample <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
    geom_point(size = 0.18, alpha = 0.7, stroke = 0) +
    coord_equal() +
    scale_color_manual(values = sample_pal) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Sample", title = "GSE157079 snATAC UMAP by sample") +
    theme_umap()
  ggsave(file.path(params$out_dir, "gse157079_umap_by_sample.png"), p_sample, width = 8, height = 6, dpi = 300)
  ggsave(file.path(params$out_dir, "gse157079_umap_by_sample.pdf"), p_sample, width = 8, height = 6)
  cat("Sample field used: ", sample_col, "\n", sep = "")
} else {
  cat("No sample/batch/condition field detected; sample UMAP skipped.\n")
}

cat("UMAP columns used: ", paste(umap_cols, collapse = ", "), "\n", sep = "")
cat("Cell type field used: ", celltype_col, "\n", sep = "")
if (!is.null(barcode_col)) cat("Barcode field detected: ", barcode_col, "\n", sep = "")
cat("Cells plotted: ", nrow(plot_df), "\n", sep = "")
cat("Saved UMAP figures to: ", normalizePath(params$out_dir, mustWork = FALSE), "\n", sep = "")
