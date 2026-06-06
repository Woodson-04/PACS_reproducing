#!/usr/bin/env Rscript

# Compare PACS-paper-style normalized batch mixing scores across existing P56
# embeddings and the GEO-provided precomputed UMAP reference. This script only
# reads existing CSV files.

suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

cmd <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(script_arg) > 0) {
  sub("^--file=", "", script_arg[[1]])
} else {
  "scripts/mouse_kidney_figures/09b_compare_paper_style_batch_mixing_scores.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "09_compute_paper_style_batch_mixing_score.R"))

parse_args_09b <- function(defaults) {
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

params <- parse_args_09b(list(
  metadata_csv = "results/mouse_kidney_figures/gse157079_metadata_merged.csv",
  out_dir = "results/mouse_kidney_figures/paper_style_batch_mixing_score",
  fig_dir = "figures/mouse_kidney",
  k = 30L,
  seed = 1L
))

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

safe_compute <- function(df, space_name, out_stub, batch_col = "batch", coord_prefix = "") {
  message("Computing: ", space_name)
  result <- compute_score_df(
    df = df,
    space_name = space_name,
    batch_col = batch_col,
    celltype_col = "cell_type",
    k = params$k,
    coord_prefix = coord_prefix,
    seed = params$seed
  )
  out_csv <- file.path(params$out_dir, paste0(out_stub, ".csv"))
  write_score_outputs(result, out_csv)
  result$summary
}

rows <- list()
missing_inputs <- character()

metadata <- read_csv(params$metadata_csv)
p56_geo <- metadata[metadata$sample %in% c("P56_batch1", "P56_batch2"), , drop = FALSE]
rows[[length(rows) + 1L]] <- safe_compute(
  p56_geo,
  "GEO precomputed UMAP P56, UMAP space",
  "geo_precomputed_umap_p56",
  batch_col = "sample"
)

settings <- list(
  list(
    label = "P56 10000/10000 before",
    path = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/p56_before_lsi_umap_embedding.csv",
    stub = "p56_top10000_test10000_before"
  ),
  list(
    label = "P56 10000/10000 after",
    path = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/p56_after_lsi_umap_embedding.csv",
    stub = "p56_top10000_test10000_after"
  ),
  list(
    label = "P56 20000/10000 before",
    path = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005/p56_before_lsi_umap_embedding.csv",
    stub = "p56_top20000_test10000_before"
  ),
  list(
    label = "P56 20000/10000 after",
    path = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005/p56_after_lsi_umap_embedding.csv",
    stub = "p56_top20000_test10000_after"
  )
)

for (setting in settings) {
  if (!file.exists(setting$path)) {
    missing_inputs <- c(missing_inputs, setting$path)
    next
  }
  df <- read_csv(setting$path)
  rows[[length(rows) + 1L]] <- safe_compute(df, setting$label, setting$stub, batch_col = "batch")
}

summary_df <- do.call(rbind, rows)
summary_out <- file.path(params$out_dir, "paper_style_batch_mixing_scores.csv")
write_csv(summary_df, summary_out)

plot_df <- summary_df
plot_df$space_name <- factor(plot_df$space_name, levels = plot_df$space_name)
p <- ggplot(plot_df, aes(x = space_name, y = normalized_batch_mixing_score, fill = coordinate_space_type)) +
  geom_col(width = 0.68) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#666666") +
  scale_fill_manual(values = c("LSI-space" = "#2E86AB", "UMAP-space" = "#A000FF"), drop = FALSE) +
  labs(
    x = NULL,
    y = "Observed / expected batch mixing score",
    fill = "Coordinate space",
    title = "PACS-paper-style normalized batch mixing score"
  ) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))
plot_base <- file.path(params$fig_dir, "gse157079_paper_style_batch_mixing_score_comparison")
ggsave(paste0(plot_base, ".png"), p, width = 10, height = 5.8, dpi = 300)
ggsave(paste0(plot_base, ".pdf"), p, width = 10, height = 5.8)

all_umap <- all(summary_df$coordinate_space_type == "UMAP-space")
report_path <- file.path(params$out_dir, "paper_style_batch_mixing_score_report.md")
report <- c(
  "# PACS-paper-style normalized batch mixing score",
  "",
  "## Goal",
  "",
  "This report computes the PACS-paper-style normalized batch mixing score for existing P56 embeddings and the GEO-provided precomputed UMAP reference.",
  "",
  "Higher normalized batch mixing score indicates better batch mixing. A value near 1 means the observed different-batch neighbor fraction approaches the expectation from the cell_type-by-batch composition matrix.",
  "",
  "## Coordinate-space caveat",
  "",
  if (all_umap) {
    "No LSI coordinate columns were found in the available P56 embedding CSV files. Therefore, all scores here are UMAP-space approximations, not PCA/LSI-space scores. They should not be directly compared numerically to the PACS paper's PCA-space 0.122 -> 0.358 statement."
  } else {
    "LSI-space scores were computed where LSI columns were available. UMAP-space rows are visualization-space approximations and should not be directly compared numerically to PACS paper PCA-space scores."
  },
  "",
  "## Summary scores",
  "",
  "```text",
  paste(capture.output(print(summary_df)), collapse = "\n"),
  "```",
  "",
  "## Interpretation",
  "",
  "- The GEO row is a public GSE157079 atlas UMAP reference embedding, not PACS-filtered ground truth.",
  "- Compare before vs after rows within the same P56 setting to assess whether PACS filtering improves normalized batch mixing.",
  "- Because the current available P56 CSVs contain UMAP coordinates but no LSI_2:LSI_30 columns, the current comparison is UMAP-space unless future scripts save LSI coordinates.",
  "",
  "## Missing inputs",
  "",
  if (length(missing_inputs) == 0) "No expected input files were missing." else paste(paste0("- `", missing_inputs, "`"), collapse = "\n"),
  "",
  "## Output files",
  "",
  paste0("- `", summary_out, "`"),
  paste0("- `", plot_base, ".png/pdf`"),
  "- Per-cell CSV files were also written for each computed row."
)
writeLines(report, report_path)

message("Saved summary: ", summary_out)
message("Saved report: ", report_path)
message("Saved plot: ", plot_base, ".png/pdf")
print(summary_df)
