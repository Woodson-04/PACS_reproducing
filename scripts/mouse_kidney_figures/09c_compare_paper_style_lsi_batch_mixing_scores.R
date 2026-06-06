#!/usr/bin/env Rscript

# Compute PACS-paper-style normalized batch mixing scores in LSI space for the
# current main P56 setting with saved LSI coordinates, and compare them with the
# existing UMAP-space approximation.

suppressPackageStartupMessages({
  library(ggplot2)
})

cmd <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(script_arg) > 0) {
  sub("^--file=", "", script_arg[[1]])
} else {
  "scripts/mouse_kidney_figures/09c_compare_paper_style_lsi_batch_mixing_scores.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "09_compute_paper_style_batch_mixing_score.R"))

parse_args_09c <- function(defaults) {
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

params <- parse_args_09c(list(
  result_dir = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved",
  out_dir = "results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space",
  fig_dir = "figures/mouse_kidney",
  previous_umap_scores = "results/mouse_kidney_figures/paper_style_batch_mixing_score/paper_style_batch_mixing_scores.csv",
  k = 30L,
  seed = 1L
))

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

compute_one <- function(path, space_name, stub) {
  if (!file.exists(path)) stop("Missing LSI embedding file: ", path)
  df <- read_csv(path)
  result <- compute_score_df(
    df = df,
    space_name = space_name,
    batch_col = "batch",
    celltype_col = "cell_type",
    k = params$k,
    coord_prefix = "",
    seed = params$seed
  )
  out_csv <- file.path(params$out_dir, paste0(stub, ".csv"))
  write_score_outputs(result, out_csv)
  result$summary
}

before_path <- file.path(params$result_dir, "p56_before_lsi_embedding.csv")
after_path <- file.path(params$result_dir, "p56_after_lsi_embedding.csv")

before_lsi <- compute_one(
  before_path,
  "P56 10000/10000 before, LSI space",
  "p56_10000_before_lsi_space"
)
after_lsi <- compute_one(
  after_path,
  "P56 10000/10000 after, LSI space",
  "p56_10000_after_lsi_space"
)

lsi_scores <- rbind(before_lsi, after_lsi)
lsi_out <- file.path(params$out_dir, "p56_10000_lsi_space_paper_style_scores.csv")
write_csv(lsi_scores, lsi_out)

plot_rows <- list()
if (file.exists(params$previous_umap_scores)) {
  previous <- read_csv(params$previous_umap_scores)
  keep_names <- c("P56 10000/10000 before", "P56 10000/10000 after")
  previous <- previous[previous$space_name %in% keep_names, , drop = FALSE]
  previous$space_name <- paste0(previous$space_name, ", UMAP space")
  plot_rows[[length(plot_rows) + 1L]] <- previous
}
plot_rows[[length(plot_rows) + 1L]] <- lsi_scores
plot_df <- do.call(rbind, plot_rows)
plot_df$space_name <- factor(plot_df$space_name, levels = plot_df$space_name)

p <- ggplot(plot_df, aes(x = space_name, y = normalized_batch_mixing_score, fill = coordinate_space_type)) +
  geom_col(width = 0.68) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#666666") +
  scale_fill_manual(values = c("LSI-space" = "#2E86AB", "UMAP-space" = "#A000FF"), drop = FALSE) +
  labs(
    x = NULL,
    y = "Observed / expected batch mixing score",
    fill = "Coordinate space",
    title = "P56 10000/10000 paper-style batch mixing: LSI vs UMAP"
  ) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))
plot_base <- file.path(params$fig_dir, "gse157079_p56_10000_paper_style_lsi_vs_umap_score")
ggsave(paste0(plot_base, ".png"), p, width = 10, height = 5.8, dpi = 300)
ggsave(paste0(plot_base, ".pdf"), p, width = 10, height = 5.8)

report_path <- file.path(params$out_dir, "p56_10000_lsi_space_paper_style_score_report.md")
report <- c(
  "# P56 10000/10000 LSI-space paper-style batch mixing score",
  "",
  "## Goal",
  "",
  "This report computes the PACS-paper-style normalized batch mixing score in LSI space for the current main P56 PACS filtering setting.",
  "",
  "## Coordinate-space interpretation",
  "",
  "The PACS paper calculates the normalized batch mixing score in PCA space. For this scATAC pipeline, LSI space, especially LSI dimensions 2:30, is the closest analogue. UMAP-space scores are visualization-space approximations and are included only as reference.",
  "",
  "## LSI-space scores",
  "",
  "```text",
  paste(capture.output(print(lsi_scores)), collapse = "\n"),
  "```",
  "",
  "## LSI vs existing UMAP-space scores",
  "",
  "```text",
  paste(capture.output(print(plot_df[, c("space_name", "coordinate_space_type", "observed_batch_mixing_score", "expected_batch_mixing_score", "normalized_batch_mixing_score", "coordinate_columns_used")])), collapse = "\n"),
  "```",
  "",
  "## Interpretation",
  "",
  "- Higher normalized batch mixing score indicates better batch mixing.",
  "- If the LSI-space score increases from before to after, this supports PACS filtering in a space closer to the PACS paper PCA-space calculation.",
  "- Do not directly compare the GEO UMAP-space reference score with these LSI-space PACS scores.",
  "",
  "## Output files",
  "",
  paste0("- `", lsi_out, "`"),
  paste0("- `", plot_base, ".png/pdf`"),
  "- Per-cell CSV files were written for before and after LSI-space scores."
)
writeLines(report, report_path)

message("Saved LSI-space scores: ", lsi_out)
message("Saved report: ", report_path)
message("Saved plot: ", plot_base, ".png/pdf")
print(lsi_scores)
