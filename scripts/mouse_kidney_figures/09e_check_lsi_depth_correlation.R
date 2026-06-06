#!/usr/bin/env Rscript

# Check correlation between P56 LSI components and cell depth for the current
# main PACS-filtered result. This script only reads existing LSI embedding CSVs.

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
  result_dir = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved",
  out_dir = "results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space",
  fig_dir = "figures/mouse_kidney"
))

read_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, header = TRUE, showProgress = FALSE, check.names = FALSE)))
  }
  read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

detect_lsi_cols <- function(df) {
  cols <- grep("^LSI_[0-9]+$", names(df), value = TRUE)
  if (length(cols) == 0) stop("No LSI columns detected")
  idx <- as.integer(sub("^LSI_", "", cols))
  cols[order(idx)]
}

detect_depth_col <- function(df, preferred) {
  candidates <- unique(c(preferred, "cell_depth", "pilot_depth", "depth", "full_matrix_depth"))
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    stop("Could not detect depth column. Tried: ", paste(candidates, collapse = ", "))
  }
  hit[[1]]
}

compute_correlations <- function(path, stage, preferred_depth_col) {
  df <- read_csv(path)
  lsi_cols <- detect_lsi_cols(df)
  depth_col <- detect_depth_col(df, preferred_depth_col)
  depth <- as.numeric(df[[depth_col]])
  rows <- lapply(lsi_cols, function(col) {
    x <- as.numeric(df[[col]])
    data.frame(
      stage = stage,
      component = col,
      component_index = as.integer(sub("^LSI_", "", col)),
      depth_col = depth_col,
      pearson_cor = suppressWarnings(cor(x, depth, method = "pearson", use = "complete.obs")),
      spearman_cor = suppressWarnings(cor(x, depth, method = "spearman", use = "complete.obs")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

plot_stage <- function(df, stage, out_base) {
  sub <- df[df$stage == stage, , drop = FALSE]
  long <- rbind(
    data.frame(component_index = sub$component_index, method = "Pearson", value = sub$pearson_cor),
    data.frame(component_index = sub$component_index, method = "Spearman", value = sub$spearman_cor)
  )
  p <- ggplot(long, aes(x = component_index, y = value, color = method)) +
    geom_hline(yintercept = 0, color = "#999999", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    labs(
      x = "LSI component",
      y = "Correlation with depth",
      color = "Method",
      title = paste0("P56 ", stage, " LSI-depth correlation")
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, ".png"), p, width = 8.5, height = 5.2, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 8.5, height = 5.2)
}

plot_combined <- function(df, out_base) {
  long <- rbind(
    data.frame(stage = df$stage, component_index = df$component_index, method = "Pearson", value = df$pearson_cor),
    data.frame(stage = df$stage, component_index = df$component_index, method = "Spearman", value = df$spearman_cor)
  )
  p <- ggplot(long, aes(x = component_index, y = value, color = stage)) +
    geom_hline(yintercept = 0, color = "#999999", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ method, nrow = 1) +
    labs(
      x = "LSI component",
      y = "Correlation with depth",
      color = "Stage",
      title = "P56 before/after LSI-depth correlation"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, ".png"), p, width = 10, height = 5.2, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 10, height = 5.2)
}

strongest_row <- function(df, stage, method_col) {
  sub <- df[df$stage == stage, , drop = FALSE]
  sub[which.max(abs(sub[[method_col]])), , drop = FALSE]
}

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

before_path <- file.path(params$result_dir, "p56_before_lsi_embedding.csv")
after_path <- file.path(params$result_dir, "p56_after_lsi_embedding.csv")
if (!file.exists(before_path)) stop("Missing before LSI embedding: ", before_path)
if (!file.exists(after_path)) stop("Missing after LSI embedding: ", after_path)

before <- compute_correlations(before_path, "before", "before_depth")
after <- compute_correlations(after_path, "after", "after_depth")
cor_df <- rbind(before, after)

out_csv <- file.path(params$out_dir, "p56_lsi_depth_correlation.csv")
write.csv(cor_df, out_csv, row.names = FALSE)

plot_stage(cor_df, "before", file.path(params$fig_dir, "gse157079_p56_lsi_depth_correlation_before"))
plot_stage(cor_df, "after", file.path(params$fig_dir, "gse157079_p56_lsi_depth_correlation_after"))
plot_combined(cor_df, file.path(params$fig_dir, "gse157079_p56_lsi_depth_correlation_before_after"))

before_spear <- strongest_row(cor_df, "before", "spearman_cor")
after_spear <- strongest_row(cor_df, "after", "spearman_cor")
before_pear <- strongest_row(cor_df, "before", "pearson_cor")
after_pear <- strongest_row(cor_df, "after", "pearson_cor")

lsi1_before <- cor_df[cor_df$stage == "before" & cor_df$component == "LSI_1", , drop = FALSE]
lsi1_after <- cor_df[cor_df$stage == "after" & cor_df$component == "LSI_1", , drop = FALSE]
lsi1_dominant <- identical(before_spear$component[[1]], "LSI_1") || identical(after_spear$component[[1]], "LSI_1") ||
  identical(before_pear$component[[1]], "LSI_1") || identical(after_pear$component[[1]], "LSI_1")

interpretation <- if (lsi1_dominant) {
  "LSI_1 is among the strongest depth-associated components. This supports the use of LSI_2:LSI_30 for batch mixing score calculation."
} else {
  "LSI_1 was excluded following common scATAC LSI practice; depth correlation check did not show LSI_1 as uniquely dominant."
}

report <- c(
  "# P56 LSI-depth correlation report",
  "",
  "## Goal",
  "",
  "This report checks whether LSI_1 is depth-associated in the current P56 PACS main result.",
  "",
  "## Input files",
  "",
  paste0("- before: `", before_path, "`"),
  paste0("- after: `", after_path, "`"),
  "",
  "## Strongest absolute correlations",
  "",
  paste0("- before strongest Spearman: ", before_spear$component, " = ", signif(before_spear$spearman_cor, 4)),
  paste0("- after strongest Spearman: ", after_spear$component, " = ", signif(after_spear$spearman_cor, 4)),
  paste0("- before strongest Pearson: ", before_pear$component, " = ", signif(before_pear$pearson_cor, 4)),
  paste0("- after strongest Pearson: ", after_pear$component, " = ", signif(after_pear$pearson_cor, 4)),
  "",
  "## LSI_1 correlations",
  "",
  paste0("- before LSI_1 Spearman: ", signif(lsi1_before$spearman_cor, 4), "; Pearson: ", signif(lsi1_before$pearson_cor, 4)),
  paste0("- after LSI_1 Spearman: ", signif(lsi1_after$spearman_cor, 4), "; Pearson: ", signif(lsi1_after$pearson_cor, 4)),
  "",
  "## Interpretation",
  "",
  interpretation,
  "",
  "## Output files",
  "",
  paste0("- `", out_csv, "`"),
  "- `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before.png/pdf`",
  "- `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_after.png/pdf`",
  "- `figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before_after.png/pdf`"
)
report_path <- file.path(params$out_dir, "p56_lsi_depth_correlation_report.md")
writeLines(report, report_path)

cat("Saved correlation CSV: ", out_csv, "\n", sep = "")
cat("Saved report: ", report_path, "\n", sep = "")
print(cor_df[cor_df$component_index <= 5, ])
