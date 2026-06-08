#!/usr/bin/env Rscript

# PCA-logNorm paper-style normalized batch mixing scores for existing P56 PACS
# parameter settings. This script reuses saved sparse matrices and retained
# peak indices. It does not rerun PACS, UMAP, MatrixMarket streaming, TF-IDF, or
# LSI.

suppressPackageStartupMessages({
  library(Matrix)
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
  base_results_dir = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures",
  out_dir = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space_parameter_stability",
  fig_dir = "/home/woodson/PACS_reproducing/figures/mouse_kidney",
  k = 30L,
  pcs_list = "20,30,50",
  seed = 1L
))

for (pkg in c("data.table", "irlba", "FNN")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package ", pkg, " is required for PCA parameter stability scoring.")
  }
}

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""), "\n", sep = "")
}

pcs_values <- as.integer(strsplit(params$pcs_list, ",")[[1]])
pcs_values <- sort(unique(pcs_values[!is.na(pcs_values)]))
if (length(pcs_values) == 0) stop("pcs_list did not contain valid integers")

settings <- data.frame(
  setting_label = c("5000/5000", "10000/10000", "20000/10000"),
  result_dir_name = c(
    "gse157079_p56_pacs_batch_filter_umap",
    "gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved",
    "gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005"
  ),
  n_top_peaks = c(10000L, 10000L, 20000L),
  tested_peaks = c(5000L, 10000L, 10000L),
  fdr_cutoff = c(0.05, 0.05, 0.05),
  significant_batch_peaks = c(3140L, 6305L, 6208L),
  retained_peaks = c(1860L, 3695L, 3792L),
  stringsAsFactors = FALSE
)
settings$result_dir <- file.path(params$base_results_dir, settings$result_dir_name)

required_files <- c(
  counts = "p56_counts_top_peaks_sparse.rds",
  metadata = "p56_metadata.csv",
  retained = "p56_retained_peak_indices.csv"
)

read_csv <- function(path) {
  as.data.frame(data.table::fread(path, header = TRUE, showProgress = FALSE, check.names = FALSE))
}

audit_setting <- function(row) {
  present <- sapply(required_files, function(f) file.exists(file.path(row$result_dir, f)))
  retained_has_col <- FALSE
  if (present[["retained"]]) {
    retained_head <- read_csv(file.path(row$result_dir, required_files[["retained"]]))
    retained_has_col <- "top_peak_col" %in% names(retained_head)
  }
  data.frame(
    setting_label = row$setting_label,
    result_dir = row$result_dir,
    n_top_peaks = row$n_top_peaks,
    tested_peaks = row$tested_peaks,
    fdr_cutoff = row$fdr_cutoff,
    significant_batch_peaks = row$significant_batch_peaks,
    retained_peaks = row$retained_peaks,
    counts_rds_exists = present[["counts"]],
    retained_indices_exists = present[["retained"]],
    metadata_exists = present[["metadata"]],
    retained_has_top_peak_col = retained_has_col,
    can_compute_directly = all(present) && retained_has_col,
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, lapply(seq_len(nrow(settings)), function(i) audit_setting(settings[i, ])))
write.csv(audit, file.path(params$out_dir, "pca_parameter_stability_input_audit.csv"), row.names = FALSE)

normalize_log_sparse <- function(counts) {
  counts <- as(counts, "dgCMatrix")
  counts@x <- rep(1, length(counts@x))
  depth <- as.numeric(Matrix::rowSums(counts))
  keep <- depth > 0
  if (!all(keep)) {
    counts <- counts[keep, , drop = FALSE]
    depth <- depth[keep]
  }
  median_depth <- stats::median(depth)
  norm <- Diagonal(x = median_depth / depth) %*% counts
  norm@x <- log1p(norm@x)
  list(matrix = as(norm, "dgCMatrix"), keep_cells = keep, depth = depth)
}

run_pca <- function(mat, n_pcs, seed) {
  set.seed(seed)
  n_pcs <- min(n_pcs, nrow(mat) - 1L, ncol(mat) - 1L)
  if (n_pcs < 2L) stop("Matrix too small for PCA")
  pca <- irlba::prcomp_irlba(mat, n = n_pcs, center = TRUE, scale. = FALSE)
  scores <- pca$x
  colnames(scores) <- paste0("PC_", seq_len(ncol(scores)))
  scores
}

expected_batch_mixing <- function(meta) {
  tab <- table(meta$cell_type, meta$batch)
  n <- nrow(meta)
  total <- 0
  for (a in rownames(tab)) {
    row_total <- sum(tab[a, ])
    for (b in colnames(tab)) {
      m_ab <- tab[a, b]
      if (row_total > 0) {
        total <- total + m_ab * ((row_total - m_ab) / row_total)
      }
    }
  }
  as.numeric(total / n)
}

score_embedding <- function(embedding, meta, k) {
  if (nrow(embedding) != nrow(meta)) {
    stop("Embedding rows do not match metadata rows")
  }
  nn <- FNN::get.knn(embedding, k = k)$nn.index
  batch <- as.character(meta$batch)
  neighbor_batch <- matrix(batch[nn], nrow = nrow(nn), ncol = ncol(nn))
  cell_scores <- rowMeans(neighbor_batch != batch)
  obs <- mean(cell_scores)
  exp <- expected_batch_mixing(meta)
  list(
    observed = obs,
    expected = exp,
    normalized = obs / exp,
    per_cell = cell_scores
  )
}

compute_setting <- function(row) {
  if (!audit$can_compute_directly[audit$setting_label == row$setting_label]) {
    log_msg("Skipping unavailable setting: ", row$setting_label)
    return(NULL)
  }

  log_msg("Loading setting ", row$setting_label)
  counts <- readRDS(file.path(row$result_dir, required_files[["counts"]]))
  metadata <- read_csv(file.path(row$result_dir, required_files[["metadata"]]))
  retained <- read_csv(file.path(row$result_dir, required_files[["retained"]]))
  if (!"top_peak_col" %in% names(retained)) {
    stop(row$setting_label, ": p56_retained_peak_indices.csv lacks top_peak_col")
  }
  if (!"batch" %in% names(metadata)) {
    if ("sample" %in% names(metadata)) {
      metadata$batch <- metadata$sample
    } else {
      stop(row$setting_label, ": metadata lacks batch/sample column")
    }
  }
  if (!"cell_type" %in% names(metadata)) {
    stop(row$setting_label, ": metadata lacks cell_type column")
  }

  before_counts <- counts
  after_cols <- as.integer(retained$top_peak_col)
  after_cols <- after_cols[!is.na(after_cols)]
  after_cols <- after_cols[after_cols >= 1 & after_cols <= ncol(counts)]
  after_cols <- sort(unique(after_cols))
  after_counts <- counts[, after_cols, drop = FALSE]
  if (ncol(after_counts) != row$retained_peaks) {
    warning(row$setting_label, ": after matrix columns ", ncol(after_counts), " != expected retained peaks ", row$retained_peaks)
  }

  stages <- list(before = before_counts, after = after_counts)
  rows <- list()
  row_id <- 0L
  max_pc <- max(pcs_values)

  for (stage in names(stages)) {
    log_msg(row$setting_label, " ", stage, ": normalization")
    norm <- normalize_log_sparse(stages[[stage]])
    meta_stage <- metadata[norm$keep_cells, , drop = FALSE]
    log_msg(row$setting_label, " ", stage, ": PCA")
    pca_scores <- run_pca(norm$matrix, max_pc, params$seed)

    for (pc_dim in pcs_values) {
      use_dim <- min(pc_dim, ncol(pca_scores))
      emb <- pca_scores[, seq_len(use_dim), drop = FALSE]
      sc <- score_embedding(emb, meta_stage, params$k)
      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        setting_label = row$setting_label,
        result_dir = row$result_dir,
        n_top_peaks = row$n_top_peaks,
        tested_peaks = row$tested_peaks,
        fdr_cutoff = row$fdr_cutoff,
        significant_batch_peaks = row$significant_batch_peaks,
        retained_peaks = row$retained_peaks,
        stage = stage,
        coordinate_space_type = "PCA-logNorm",
        pc_dims = use_dim,
        n_cells = nrow(meta_stage),
        n_features = ncol(stages[[stage]]),
        observed_batch_mixing_score = sc$observed,
        expected_batch_mixing_score = sc$expected,
        normalized_batch_mixing_score = sc$normalized,
        k = params$k,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

available_settings <- settings[audit$can_compute_directly, , drop = FALSE]
if (nrow(available_settings) == 0) {
  stop("No settings have all required files for direct PCA-logNorm scoring")
}

score_list <- lapply(seq_len(nrow(available_settings)), function(i) compute_setting(available_settings[i, ]))
scores <- do.call(rbind, score_list)
scores_out <- file.path(params$out_dir, "p56_pca_parameter_stability_scores.csv")
write.csv(scores, scores_out, row.names = FALSE)

plot_df <- scores
plot_df$setting_label <- factor(plot_df$setting_label, levels = c("5000/5000", "10000/10000", "20000/10000"))
plot_df$stage <- factor(plot_df$stage, levels = c("before", "after"))
plot_df$pc_label <- paste0("PC1:", plot_df$pc_dims)

p <- ggplot(plot_df, aes(x = setting_label, y = normalized_batch_mixing_score, color = stage, group = stage)) +
  geom_point(size = 2.8) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ pc_label, nrow = 1) +
  scale_color_manual(values = c(before = "#7A8793", after = "#1F77B4")) +
  labs(
    title = "P56 PCA-logNorm paper-style batch mixing parameter stability",
    x = "PACS setting",
    y = "Normalized batch mixing score",
    color = "Stage"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    strip.background = element_rect(fill = "#eef2f6", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

plot_png <- file.path(params$fig_dir, "gse157079_p56_pca_space_parameter_stability_scores.png")
plot_pdf <- file.path(params$fig_dir, "gse157079_p56_pca_space_parameter_stability_scores.pdf")
ggsave(plot_png, p, width = 8.5, height = 4.4, dpi = 300)
ggsave(plot_pdf, p, width = 8.5, height = 4.4)

format_table <- function(x) paste(capture.output(print(x, row.names = FALSE)), collapse = "\n")
pc30 <- scores[scores$pc_dims == 30, ]
pc50 <- scores[scores$pc_dims == 50, ]

report_path <- file.path(params$out_dir, "p56_pca_parameter_stability_report.md")
report <- c(
  "# P56 PCA-logNorm Parameter Stability Report",
  "",
  "## Input Audit Summary",
  "",
  "```text",
  format_table(audit),
  "```",
  "",
  "## Computed Settings",
  "",
  paste("- ", available_settings$setting_label, collapse = "\n"),
  "",
  "## Score Table",
  "",
  "```text",
  format_table(scores[, c("setting_label", "stage", "pc_dims", "n_cells", "n_features", "normalized_batch_mixing_score")]),
  "```",
  "",
  "## PC1:30 Comparison",
  "",
  "```text",
  format_table(pc30[, c("setting_label", "stage", "normalized_batch_mixing_score")]),
  "```",
  "",
  "## PC1:50 Comparison",
  "",
  "```text",
  format_table(pc50[, c("setting_label", "stage", "normalized_batch_mixing_score")]),
  "```",
  "",
  "## Interpretation",
  "",
  "- PCA-logNorm parameter stability is more consistent with the PACS paper-style normalized PCA mixing metric than UMAP-space auxiliary metrics.",
  "- If before-to-after improvement holds across settings, this supports PACS filtering beyond visualization-space UMAP.",
  "- The preferred setting should balance batch mixing improvement and biological structure preservation.",
  "- Compare this PCA-space analysis with the existing UMAP/LSI summaries before updating Section 10.",
  "",
  "## Output Files",
  "",
  paste0("- `", scores_out, "`"),
  paste0("- `", plot_png, "`"),
  paste0("- `", plot_pdf, "`")
)
writeLines(report, report_path)

section_path <- file.path(params$out_dir, "section10_pca_parameter_stability_replacement.html")
pc30_wide <- reshape(
  pc30[, c("setting_label", "stage", "normalized_batch_mixing_score")],
  idvar = "setting_label",
  timevar = "stage",
  direction = "wide"
)
pc50_wide <- reshape(
  pc50[, c("setting_label", "stage", "normalized_batch_mixing_score")],
  idvar = "setting_label",
  timevar = "stage",
  direction = "wide"
)
merged <- merge(pc30_wide, pc50_wide, by = "setting_label", suffixes = c("_pc30", "_pc50"))

html_rows <- apply(merged, 1, function(r) {
  paste0(
    "<tr><td>", r["setting_label"], "</td>",
    "<td>", sprintf("%.4f", as.numeric(r["normalized_batch_mixing_score.before_pc30"])), "</td>",
    "<td>", sprintf("%.4f", as.numeric(r["normalized_batch_mixing_score.after_pc30"])), "</td>",
    "<td>", sprintf("%.4f", as.numeric(r["normalized_batch_mixing_score.before_pc50"])), "</td>",
    "<td>", sprintf("%.4f", as.numeric(r["normalized_batch_mixing_score.after_pc50"])), "</td></tr>"
  )
})
section_html <- c(
  "<h2>10. GEO Reference 与 PCA-space 参数稳定性</h2>",
  "<p>GEO precomputed UMAP 可作为 public atlas/reference embedding，但它可能包含 atlas-level preprocessing 或 integration/correction steps，不应解释为 PACS-filtered ground truth。</p>",
  "<table><thead><tr><th>setting</th><th>UMAP-space normalized score</th></tr></thead><tbody>",
  "<tr><td>P56 10000/10000 before</td><td>0.02649</td></tr>",
  "<tr><td>P56 10000/10000 after</td><td>0.65073</td></tr>",
  "<tr><td>GEO P56 UMAP-space</td><td>0.93903</td></tr>",
  "</tbody></table>",
  "<p>下面的 PCA-logNorm 参数稳定性更接近 PACS paper-style normalized PCA mixing metric，因此应优先用于参数稳定性表述。</p>",
  "<table><thead><tr><th>setting</th><th>PC1:30 before</th><th>PC1:30 after</th><th>PC1:50 before</th><th>PC1:50 after</th></tr></thead><tbody>",
  html_rows,
  "</tbody></table>",
  "<img class=\"fig-small\" alt=\"P56 PCA-space parameter stability\" src=\"materials/figures/gse157079_p56_pca_space_parameter_stability_scores.png\" />",
  "<p>若各 setting 均显示 after 高于 before，则说明 PACS filtering 的 batch mixing 改善不仅存在于 UMAP visualization space，也存在于更接近论文 normalized PCA mixing 的 PCA-logNorm space。</p>",
  "<p>GEO UMAP 仍作为 public reference embedding 单独解释，不作为 PACS-filtered ground truth。</p>"
)
writeLines(section_html, section_path)

log_msg("Saved scores: ", scores_out)
log_msg("Saved plot: ", plot_png)
log_msg("Saved report: ", report_path)
log_msg("Saved Section 10 replacement HTML: ", section_path)
