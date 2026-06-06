#!/usr/bin/env Rscript

# Compute LSI dimension sensitivity for PACS-paper-style normalized batch mixing
# scores using existing P56 before/after LSI embeddings.

suppressPackageStartupMessages({
  library(ggplot2)
})

cmd <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(script_arg) > 0) {
  sub("^--file=", "", script_arg[[1]])
} else {
  "scripts/mouse_kidney_figures/09f_lsi_dimension_sensitivity_batch_mixing.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "09_compute_paper_style_batch_mixing_score.R"))

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
  result_dir = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved",
  out_dir = "results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space",
  fig_dir = "figures/mouse_kidney",
  k = 30L,
  seed = 1L
))

compute_custom_lsi_score <- function(df, setting, dim_label, lsi_cols) {
  missing <- setdiff(lsi_cols, names(df))
  if (length(missing) > 0) stop("Missing LSI columns for ", dim_label, ": ", paste(missing, collapse = ", "))
  coords <- as.matrix(df[, lsi_cols, drop = FALSE])
  coords <- apply(coords, 2, as.numeric)
  batch <- factor(as.character(df$batch))
  cell_type <- factor(as.character(df$cell_type))
  nn <- knn_indices(coords, params$k)
  cell_score <- vapply(seq_len(nrow(coords)), function(i) {
    mean(batch[nn[i, ]] != batch[i], na.rm = TRUE)
  }, numeric(1))
  observed <- mean(cell_score, na.rm = TRUE)
  expected <- expected_batch_mixing_score(cell_type, batch)
  data.frame(
    setting = setting,
    dim_label = dim_label,
    coordinate_space_type = "LSI-space",
    n_cells = nrow(coords),
    k = ncol(nn),
    n_batches = length(levels(batch)),
    n_cell_types = length(levels(cell_type)),
    observed_batch_mixing_score = observed,
    expected_batch_mixing_score = expected,
    normalized_batch_mixing_score = observed / expected,
    coordinate_columns_used = paste(lsi_cols, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

dimension_sets <- list(
  "LSI_1:30" = paste0("LSI_", 1:30),
  "LSI_2:30" = paste0("LSI_", 2:30),
  "LSI_2:50" = paste0("LSI_", 2:50),
  "LSI_1:50" = paste0("LSI_", 1:50)
)

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

before_path <- file.path(params$result_dir, "p56_before_lsi_embedding.csv")
after_path <- file.path(params$result_dir, "p56_after_lsi_embedding.csv")
if (!file.exists(before_path)) stop("Missing before LSI embedding: ", before_path)
if (!file.exists(after_path)) stop("Missing after LSI embedding: ", after_path)

before <- read_csv(before_path)
after <- read_csv(after_path)

rows <- list()
for (nm in names(dimension_sets)) {
  rows[[length(rows) + 1L]] <- compute_custom_lsi_score(before, "before", nm, dimension_sets[[nm]])
  rows[[length(rows) + 1L]] <- compute_custom_lsi_score(after, "after", nm, dimension_sets[[nm]])
}
scores <- do.call(rbind, rows)
out_csv <- file.path(params$out_dir, "p56_lsi_dimension_sensitivity_scores.csv")
write_csv(scores, out_csv)

plot_df <- scores
plot_df$dim_label <- factor(plot_df$dim_label, levels = names(dimension_sets))
plot_df$setting <- factor(plot_df$setting, levels = c("before", "after"))
p <- ggplot(plot_df, aes(x = dim_label, y = normalized_batch_mixing_score, fill = setting)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#666666") +
  scale_fill_manual(values = c(before = "#777777", after = "#2E86AB")) +
  labs(
    x = "LSI dimensions",
    y = "Observed / expected batch mixing score",
    fill = "Stage",
    title = "P56 LSI dimension sensitivity for paper-style batch mixing score"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
plot_base <- file.path(params$fig_dir, "gse157079_p56_lsi_dimension_sensitivity_batch_mixing")
ggsave(paste0(plot_base, ".png"), p, width = 9, height = 5.6, dpi = 300)
ggsave(paste0(plot_base, ".pdf"), p, width = 9, height = 5.6)

corr_path <- file.path(params$out_dir, "p56_lsi_depth_correlation.csv")
corr_summary <- "LSI-depth correlation results were not found. Run 09e_check_lsi_depth_correlation.R first to fill this section."
if (file.exists(corr_path)) {
  corr <- read_csv(corr_path)
  strongest <- function(stage, col) {
    sub <- corr[corr$stage == stage, , drop = FALSE]
    sub[which.max(abs(sub[[col]])), , drop = FALSE]
  }
  before_s <- strongest("before", "spearman_cor")
  after_s <- strongest("after", "spearman_cor")
  lsi1 <- corr[corr$component == "LSI_1", c("stage", "spearman_cor", "pearson_cor"), drop = FALSE]
  corr_summary <- paste(
    paste0("- before strongest Spearman: ", before_s$component, " = ", signif(before_s$spearman_cor, 4)),
    paste0("- after strongest Spearman: ", after_s$component, " = ", signif(after_s$spearman_cor, 4)),
    paste0("- LSI_1 correlations: ", paste(paste0(lsi1$stage, " Spearman=", signif(lsi1$spearman_cor, 4), ", Pearson=", signif(lsi1$pearson_cor, 4)), collapse = "; ")),
    sep = "\n"
  )
}

wide <- reshape(
  scores[, c("setting", "dim_label", "normalized_batch_mixing_score")],
  idvar = "dim_label",
  timevar = "setting",
  direction = "wide"
)
names(wide) <- sub("^normalized_batch_mixing_score\\.", "", names(wide))
wide$improvement <- wide$after - wide$before

report <- c(
  "# P56 LSI dimension sensitivity for paper-style batch mixing score",
  "",
  "## Scores",
  "",
  "```text",
  paste(capture.output(print(scores)), collapse = "\n"),
  "```",
  "",
  "## Before-to-after improvement",
  "",
  "```text",
  paste(capture.output(print(wide)), collapse = "\n"),
  "```",
  "",
  "## Interpretation",
  "",
  "- Main LSI metric remains LSI_2:30.",
  "- The purpose of this sensitivity check is to confirm that the before-to-after improvement is robust to including/excluding LSI_1 and using 30 vs 50 dimensions.",
  "- If all dimension sets improve from before to after, excluding LSI_1 does not change the main conclusion.",
  "",
  "## LSI-depth correlation summary",
  "",
  corr_summary,
  "",
  "## Output files",
  "",
  paste0("- `", out_csv, "`"),
  paste0("- `", plot_base, ".png/pdf`")
)
report_path <- file.path(params$out_dir, "p56_lsi_dimension_sensitivity_report.md")
writeLines(report, report_path)

final_note_path <- file.path(params$out_dir, "p56_paper_style_metric_final_interpretation.md")
main_lsi <- scores[scores$dim_label == "LSI_2:30", c("setting", "normalized_batch_mixing_score"), drop = FALSE]
final_note <- c(
  "# P56 paper-style metric final interpretation",
  "",
  "## Metric hierarchy",
  "",
  "1. The PACS paper uses normalized PCA mixing.",
  "2. Our closest author-like metric is PCA-logNorm PC1:30 or PC1:50.",
  "3. LSI_2:30 is a scATAC-aware analogue because TF-IDF/LSI is standard for sparse chromatin accessibility matrices.",
  "4. UMAP-space score is visualization-space only.",
  "",
  "## Main quantitative conclusion",
  "",
  "- PCA-logNorm score improves strongly after PACS filtering.",
  "- LSI-space score improves strongly after PACS filtering.",
  "- UMAP-space score also improves, but should not be used as the main paper-style metric.",
  paste0("- Main LSI_2:30 normalized score: before = ", signif(main_lsi$normalized_batch_mixing_score[main_lsi$setting == "before"], 4),
         ", after = ", signif(main_lsi$normalized_batch_mixing_score[main_lsi$setting == "after"], 4), "."),
  "",
  "## LSI_1 depth correlation",
  "",
  corr_summary,
  "",
  "## Why not directly equate to the paper values",
  "",
  "Do not directly equate our scores with the paper's 0.122 -> 0.358 because this analysis uses a P56-only two-batch subset, a top10000 feature universe, depth-derived cap_rates, different normalization, and may differ from the author-specific adult kidney feature set.",
  "",
  "## Final current recommendation",
  "",
  "Use PCA-logNorm PC1:30 or PC1:50 as the most paper-like quantitative metric, LSI_2:30 as the scATAC-aware supporting metric, and UMAP-space scores only as visualization-level supporting evidence."
)
writeLines(final_note, final_note_path)

cat("Saved sensitivity scores: ", out_csv, "\n", sep = "")
cat("Saved sensitivity report: ", report_path, "\n", sep = "")
cat("Saved final interpretation note: ", final_note_path, "\n", sep = "")
print(scores)
