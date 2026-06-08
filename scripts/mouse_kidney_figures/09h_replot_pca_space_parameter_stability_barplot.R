#!/usr/bin/env Rscript

# Replot P56 PACS PCA-logNorm parameter stability from an existing score CSV.
# This script does not rerun PACS, UMAP, MatrixMarket streaming, TF-IDF/LSI, or
# PCA scoring. It only reads p56_pca_parameter_stability_scores.csv and creates
# a report-ready PCA-logNorm bar plot.

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
  score_csv = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space_parameter_stability/p56_pca_parameter_stability_scores.csv",
  out_dir = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space_parameter_stability",
  fig_dir = "/home/woodson/PACS_reproducing/figures/mouse_kidney",
  report_fig_dir = "/home/woodson/PACS_reproducing/report/pdf_report/materials/figures"
))

if (!file.exists(params$score_csv)) {
  stop("Missing score CSV: ", params$score_csv)
}

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(params$score_csv, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c(
  "setting_label",
  "stage",
  "pc_dims",
  "coordinate_space_type",
  "normalized_batch_mixing_score"
)
missing_cols <- setdiff(required_cols, names(scores))
if (length(missing_cols) > 0) {
  stop("Score CSV is missing required columns: ", paste(missing_cols, collapse = ", "))
}

plot_df <- scores[scores$coordinate_space_type == "PCA-logNorm", , drop = FALSE]
if (nrow(plot_df) == 0) {
  stop("No PCA-logNorm rows found in score CSV.")
}

setting_levels <- c("5000/5000", "10000/10000", "20000/10000")
stage_levels <- c("before", "after")
plot_df <- plot_df[plot_df$setting_label %in% setting_levels & plot_df$stage %in% stage_levels, , drop = FALSE]
plot_df$setting_label <- factor(plot_df$setting_label, levels = setting_levels)
plot_df$stage <- factor(plot_df$stage, levels = stage_levels)
plot_df$pc_dims <- factor(
  paste0("PC1:", plot_df$pc_dims),
  levels = paste0("PC1:", sort(unique(as.integer(scores$pc_dims))))
)
plot_df$label <- sprintf("%.3f", plot_df$normalized_batch_mixing_score)

plot_df <- plot_df[order(plot_df$pc_dims, plot_df$setting_label, plot_df$stage), , drop = FALSE]

stage_colors <- c(
  before = "#8E9AAF",
  after = "#2F80ED"
)

p <- ggplot(plot_df, aes(
  x = setting_label,
  y = normalized_batch_mixing_score,
  fill = stage
)) +
  geom_col(
    position = position_dodge(width = 0.76),
    width = 0.64,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.76),
    vjust = -0.35,
    size = 3.2,
    color = "#263238"
  ) +
  facet_wrap(~ pc_dims, nrow = 1) +
  scale_fill_manual(values = stage_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, max(plot_df$normalized_batch_mixing_score, na.rm = TRUE) * 1.16),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "P56 PACS PCA-logNorm parameter stability",
    x = "PACS setting",
    y = "Normalized batch mixing score",
    fill = "Stage"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    strip.background = element_rect(fill = "#F3F5F7", color = "#D0D7DE"),
    strip.text = element_text(face = "bold", size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top",
    panel.spacing.x = unit(0.6, "cm")
  )

plot_png <- file.path(params$fig_dir, "gse157079_p56_pca_space_parameter_stability_scores_replot.png")
plot_pdf <- file.path(params$fig_dir, "gse157079_p56_pca_space_parameter_stability_scores_replot.pdf")
ggsave(plot_png, p, width = 8.6, height = 4.4, dpi = 300)
ggsave(plot_pdf, p, width = 8.6, height = 4.4)

copy_script <- "/home/woodson/PACS_reproducing/report/pdf_report/COPY_PCA_STABILITY_REPLOT_FIGURE.sh"
copy_ok <- FALSE
dir.create(params$report_fig_dir, recursive = TRUE, showWarnings = FALSE)
if (dir.exists(params$report_fig_dir)) {
  copy_ok <- file.copy(
    plot_png,
    file.path(params$report_fig_dir, basename(plot_png)),
    overwrite = TRUE
  )
}
if (!isTRUE(copy_ok)) {
  writeLines(c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "cd /home/woodson/PACS_reproducing",
    "mkdir -p report/pdf_report/materials/figures",
    "cp -v figures/mouse_kidney/gse157079_p56_pca_space_parameter_stability_scores_replot.png report/pdf_report/materials/figures/"
  ), copy_script)
}

note_path <- file.path(params$out_dir, "p56_pca_parameter_stability_replot_note.md")
pc30 <- plot_df[plot_df$pc_dims == "PC1:30", c("setting_label", "stage", "normalized_batch_mixing_score"), drop = FALSE]
note <- c(
  "# P56 PCA-logNorm Parameter Stability Replot Note",
  "",
  "## Input",
  "",
  paste0("- CSV: `", params$score_csv, "`"),
  "",
  "## Outputs",
  "",
  paste0("- PNG: `", plot_png, "`"),
  paste0("- PDF: `", plot_pdf, "`"),
  paste0("- Report material PNG: `", file.path(params$report_fig_dir, basename(plot_png)), "`"),
  "",
  "## Scope",
  "",
  "- This replot contains only PCA-logNorm parameter stability scores.",
  "- It does not include GEO reference data.",
  "- It does not include UMAP-space auxiliary metrics.",
  "- It reads existing PCA score CSV output and does not rerun PACS, UMAP, MatrixMarket streaming, TF-IDF/LSI, or PCA scoring.",
  "",
  "## Axis Label Handling",
  "",
  "- The x-axis labels use `angle = 35, hjust = 1, vjust = 1`, so labels tilt toward the lower right.",
  "",
  "## PC1:30 Summary",
  "",
  "```text",
  paste(capture.output(print(pc30)), collapse = "\n"),
  "```",
  "",
  "## Interpretation",
  "",
  "All three settings show after > before, supporting that PACS filtering improves batch mixing in PCA-logNorm space across parameter choices.",
  "",
  "The 5000/5000 after score is highest, but it retains the fewest peaks and may be more aggressive. The 10000/10000 setting remains the main display setting because it provides a more balanced tradeoff among batch mixing, retained features, and cell-type preservation."
)
writeLines(note, note_path)

section_path <- file.path(params$out_dir, "section10_pca_parameter_stability_replacement.html")
wide <- reshape(
  plot_df[, c("setting_label", "stage", "pc_dims", "normalized_batch_mixing_score")],
  idvar = c("setting_label", "pc_dims"),
  timevar = "stage",
  direction = "wide"
)
wide <- wide[order(wide$pc_dims, wide$setting_label), , drop = FALSE]
rows_html <- apply(wide, 1, function(r) {
  paste0(
    "<tr><td>", r["pc_dims"], "</td>",
    "<td>", r["setting_label"], "</td>",
    "<td>", sprintf("%.4f", as.numeric(r["normalized_batch_mixing_score.before"])), "</td>",
    "<td>", sprintf("%.4f", as.numeric(r["normalized_batch_mixing_score.after"])), "</td></tr>"
  )
})
section <- c(
  "<h2>10. GEO Reference 与参数稳定性</h2>",
  "<h3>10.1 GEO UMAP-space reference</h3>",
  "<p>GEO precomputed UMAP 可作为 public atlas/reference embedding，但它可能包含 atlas-level preprocessing 或 integration/correction steps，不应解释为 PACS-filtered ground truth。</p>",
  "<table><thead><tr><th>setting</th><th>UMAP-space normalized score</th></tr></thead><tbody>",
  "<tr><td>P56 10000/10000 before</td><td>0.02649</td></tr>",
  "<tr><td>P56 10000/10000 after</td><td>0.65073</td></tr>",
  "<tr><td>GEO P56 UMAP-space</td><td>0.93903</td></tr>",
  "</tbody></table>",
  "<h3>10.2 PCA-logNorm parameter stability</h3>",
  "<p>参数稳定性部分只使用 PCA-logNorm paper-style normalized batch mixing score，不包含 GEO reference，也不包含 UMAP-space auxiliary metrics。</p>",
  "<table><thead><tr><th>PC dimensions</th><th>setting</th><th>before</th><th>after</th></tr></thead><tbody>",
  rows_html,
  "</tbody></table>",
  "<img class=\"fig-small\" alt=\"P56 PCA-logNorm parameter stability\" src=\"materials/figures/gse157079_p56_pca_space_parameter_stability_scores_replot.png\" />",
  "<p>三组 setting 均显示 after &gt; before，说明 PACS filtering 对 PCA-logNorm space 中 batch mixing 的改善在参数选择上具有稳定性。5000/5000 after score 最高，但 retained peaks 最少，可能更激进；10000/10000 是当前主展示设置，因为它在 batch mixing 与 retained features / cell-type preservation 之间更平衡。</p>"
)
writeLines(section, section_path)

cat("Saved PNG: ", plot_png, "\n", sep = "")
cat("Saved PDF: ", plot_pdf, "\n", sep = "")
cat("Copied report material PNG: ", copy_ok, "\n", sep = "")
if (!isTRUE(copy_ok)) cat("Saved copy helper: ", copy_script, "\n", sep = "")
cat("Saved note: ", note_path, "\n", sep = "")
cat("Saved section HTML: ", section_path, "\n", sep = "")
