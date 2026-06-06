#!/usr/bin/env Rscript

# Fast replot for existing P56 PACS batch-filtering UMAP embeddings.
# This script only reads existing CSV outputs and redraws figures with larger
# point and legend glyph sizes. It does not rerun matrix streaming, PACS,
# TF-IDF, LSI, or UMAP.

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
    value <- args[[i + 1]]
    old <- defaults[[name]]
    if (is.integer(old)) {
      defaults[[name]] <- as.integer(value)
    } else if (is.numeric(old)) {
      defaults[[name]] <- as.numeric(value)
    } else {
      defaults[[name]] <- value
    }
    i <- i + 2
  }
  defaults
}

params <- parse_args(list(
  result_dir = "",
  fig_dir = "figures/mouse_kidney",
  run_name = "",
  point_size = 0.55,
  point_alpha = 0.85,
  legend_point_size = 4,
  seed = 1L
))

use_run_name <- nzchar(params$run_name) && params$run_name != "default"
if (!nzchar(params$result_dir)) {
  params$result_dir <- if (use_run_name) {
    file.path("results/mouse_kidney_figures", paste0("gse157079_p56_pacs_batch_filter_umap_", params$run_name))
  } else {
    "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap"
  }
}
figure_prefix_name <- if (use_run_name) {
  paste0("gse157079_p56_pacs_batch_filter_", params$run_name)
} else {
  "gse157079_p56_pacs_batch_filter"
}

read_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, header = TRUE, showProgress = FALSE, check.names = FALSE)))
  }
  read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""), "\n", sep = "")
}

find_existing_file <- function(result_dir, candidates, label) {
  paths <- file.path(result_dir, candidates)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop("Could not find ", label, " embedding file. Tried: ", paste(paths, collapse = ", "))
  }
  hit[[1]]
}

detect_umap_cols <- function(df, preferred_prefix = NULL, label = "") {
  candidates <- list()
  if (!is.null(preferred_prefix)) {
    candidates[[length(candidates) + 1L]] <- paste0(preferred_prefix, c("_umap_1", "_umap_2"))
  }
  candidates[[length(candidates) + 1L]] <- c("lsi_umap_1", "lsi_umap_2")
  candidates[[length(candidates) + 1L]] <- c("umap_1", "umap_2")
  candidates[[length(candidates) + 1L]] <- c("UMAP_1", "UMAP_2")
  for (cols in candidates) {
    if (all(cols %in% names(df))) {
      log_msg(label, " coordinate columns detected: ", paste(cols, collapse = ", "))
      return(cols)
    }
  }
  stop(
    "Could not detect UMAP coordinate columns for ", label,
    ". Available columns: ", paste(names(df), collapse = ", ")
  )
}

prepare_df <- function(df, cols, label) {
  if (!"cell_type" %in% names(df)) {
    stop(label, " embedding is missing required column: cell_type")
  }
  if ("batch" %in% names(df)) {
    df$batch_plot <- as.character(df$batch)
  } else if ("sample" %in% names(df)) {
    df$batch_plot <- as.character(df$sample)
  } else {
    stop(label, " embedding is missing both batch and sample columns")
  }
  df$cell_type <- as.character(df$cell_type)
  df$plot_umap_1 <- as.numeric(df[[cols[[1]]]])
  df$plot_umap_2 <- as.numeric(df[[cols[[2]]]])
  if (anyNA(df$plot_umap_1) || anyNA(df$plot_umap_2)) {
    stop(label, " UMAP coordinates contain NA after numeric conversion")
  }
  expected_batches <- c("P56_batch1", "P56_batch2")
  df$batch_plot <- factor(df$batch_plot, levels = expected_batches)
  if (anyNA(df$batch_plot)) {
    missing <- sort(unique(as.character(df$batch_plot[is.na(df$batch_plot)])))
    stop(label, " has unexpected batch levels. Expected P56_batch1/P56_batch2.")
  }
  df$cell_type <- factor(df$cell_type, levels = sort(unique(df$cell_type)))
  df
}

make_palette <- function(levels, preferred, fallback) {
  pal <- preferred[names(preferred) %in% levels]
  missing <- setdiff(levels, names(pal))
  if (length(missing) > 0) {
    extra <- rep(fallback, length.out = length(missing))
    names(extra) <- missing
    pal <- c(pal, extra)
  }
  pal[levels]
}

batch_palette <- c(
  "P56_batch1" = "#FFD400",
  "P56_batch2" = "#A000FF"
)

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
celltype_fallback <- c(
  "#66C7FF", "#D8896A", "#C69214", "#B6A000", "#7FB000",
  "#67B83F", "#00A651", "#00A98F", "#00A6B8", "#00A7E1",
  "#139DDF", "#8F9BEF", "#B879E8", "#DD70D6", "#F06AA8",
  "#F4A6C8", "#9ADBC5", "#E6C36A", "#BCA7FF", "#7DD3FC"
)

plot_one <- function(df, color_col, palette, title, legend_title, seed) {
  set.seed(seed)
  plot_df <- df[sample(seq_len(nrow(df))), , drop = FALSE]
  ggplot(plot_df, aes(x = plot_umap_1, y = plot_umap_2, color = .data[[color_col]])) +
    geom_point(size = params$point_size, alpha = params$point_alpha, stroke = 0) +
    scale_color_manual(values = palette, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = params$legend_point_size, alpha = 1))) +
    coord_equal(expand = TRUE) +
    labs(x = "UMAP 1", y = "UMAP 2", color = legend_title, title = title) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 17),
      axis.title = element_text(size = 14),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.key.width = grid::unit(0.45, "cm")
    )
}

save_plot_pair <- function(plot, path_base, width = 8.5, height = 6.5) {
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, dpi = 300)
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height)
}

dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)
before_file <- find_existing_file(
  params$result_dir,
  c("p56_before_lsi_umap_embedding.csv", "before_lsi_umap_embedding.csv", "before_umap_embedding.csv"),
  "before"
)
after_file <- find_existing_file(
  params$result_dir,
  c("p56_after_lsi_umap_embedding.csv", "after_lsi_umap_embedding.csv", "after_umap_embedding.csv"),
  "after"
)

log_msg("Reading before embedding: ", before_file)
before_raw <- read_csv(before_file)
log_msg("Reading after embedding: ", after_file)
after_raw <- read_csv(after_file)

before_cols <- detect_umap_cols(before_raw, "before", "before")
after_cols <- detect_umap_cols(after_raw, "after", "after")
before_df <- prepare_df(before_raw, before_cols, "before")
after_df <- prepare_df(after_raw, after_cols, "after")

celltype_levels <- sort(unique(c(as.character(before_df$cell_type), as.character(after_df$cell_type))))
celltype_palette <- make_palette(celltype_levels, celltype_palette_preferred, celltype_fallback)

prefix <- file.path(params$fig_dir, figure_prefix_name)
p_before_batch <- plot_one(
  before_df,
  "batch_plot",
  batch_palette,
  "P56 before PACS filtering by batch",
  "batch",
  params$seed
)
p_before_celltype <- plot_one(
  before_df,
  "cell_type",
  celltype_palette,
  "P56 before PACS filtering by cell type",
  "cell type",
  params$seed
)
p_after_batch <- plot_one(
  after_df,
  "batch_plot",
  batch_palette,
  "P56 after PACS batch-peak filtering by batch",
  "batch",
  params$seed
)
p_after_celltype <- plot_one(
  after_df,
  "cell_type",
  celltype_palette,
  "P56 after PACS batch-peak filtering by cell type",
  "cell type",
  params$seed
)

outputs <- c(
  paste0(prefix, "_before_by_batch.png"),
  paste0(prefix, "_before_by_batch.pdf"),
  paste0(prefix, "_before_by_celltype.png"),
  paste0(prefix, "_before_by_celltype.pdf"),
  paste0(prefix, "_after_by_batch.png"),
  paste0(prefix, "_after_by_batch.pdf"),
  paste0(prefix, "_after_by_celltype.png"),
  paste0(prefix, "_after_by_celltype.pdf")
)

save_plot_pair(p_before_batch, paste0(prefix, "_before_by_batch"))
save_plot_pair(p_before_celltype, paste0(prefix, "_before_by_celltype"))
save_plot_pair(p_after_batch, paste0(prefix, "_after_by_batch"))
save_plot_pair(p_after_celltype, paste0(prefix, "_after_by_celltype"))

combined_status <- "combined figure skipped: neither patchwork nor cowplot is available"
if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- (p_before_batch + p_before_celltype) / (p_after_batch + p_after_celltype)
  ggsave(paste0(prefix, "_four_panel.png"), combined, width = 14, height = 10, dpi = 300)
  ggsave(paste0(prefix, "_four_panel.pdf"), combined, width = 14, height = 10)
  outputs <- c(outputs, paste0(prefix, "_four_panel.png"), paste0(prefix, "_four_panel.pdf"))
  combined_status <- "combined figure overwritten with patchwork"
} else if (requireNamespace("cowplot", quietly = TRUE)) {
  combined <- cowplot::plot_grid(p_before_batch, p_before_celltype, p_after_batch, p_after_celltype, ncol = 2)
  ggsave(paste0(prefix, "_four_panel.png"), combined, width = 14, height = 10, dpi = 300)
  ggsave(paste0(prefix, "_four_panel.pdf"), combined, width = 14, height = 10)
  outputs <- c(outputs, paste0(prefix, "_four_panel.png"), paste0(prefix, "_four_panel.pdf"))
  combined_status <- "combined figure overwritten with cowplot"
}

note <- c(
  "# P56 PACS Batch-Filtering UMAP Replot Note",
  "",
  paste0("- date_time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("- before_embedding_file: `", before_file, "`"),
  paste0("- after_embedding_file: `", after_file, "`"),
  paste0("- before_coordinate_columns: ", paste(before_cols, collapse = ", ")),
  paste0("- after_coordinate_columns: ", paste(after_cols, collapse = ", ")),
  paste0("- point_size: ", params$point_size),
  paste0("- point_alpha: ", params$point_alpha),
  paste0("- legend_point_size: ", params$legend_point_size),
  "",
  "## Output figures overwritten",
  "",
  paste0("- `", outputs, "`"),
  "",
  "## Analysis note",
  "",
  "This replot did not rerun PACS, LSI, UMAP, or matrix streaming.",
  "",
  "After filtering, partial batch-associated structure is reduced, but residual batch separation remains. This is an initial P56-only top-peak PACS filtering result, not a final full correction.",
  "",
  paste0("Combined figure status: ", combined_status)
)
writeLines(note, file.path(params$result_dir, "replot_note.md"))

log_msg("Four separate figures overwritten")
log_msg(combined_status)
log_msg("Replot note saved: ", file.path(params$result_dir, "replot_note.md"))
cat("before_embedding_file=", before_file, "\n", sep = "")
cat("after_embedding_file=", after_file, "\n", sep = "")
cat("before_coordinate_columns=", paste(before_cols, collapse = ","), "\n", sep = "")
cat("after_coordinate_columns=", paste(after_cols, collapse = ","), "\n", sep = "")
cat("four_separate_plots_overwritten=TRUE\n")
cat("combined_status=", combined_status, "\n", sep = "")
cat("heavy_analysis_rerun=FALSE\n")
