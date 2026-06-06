#!/usr/bin/env Rscript

# Compute a PCA-logNorm-space version of the PACS-paper-style normalized batch
# mixing score for the current P56 main setting. This script reads existing
# sparse RDS outputs and does not rerun PACS or MatrixMarket streaming.

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

cmd <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", cmd, value = TRUE)
script_path <- if (length(script_arg) > 0) {
  sub("^--file=", "", script_arg[[1]])
} else {
  "scripts/mouse_kidney_figures/09d_compute_pca_space_batch_mixing_score.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "09_compute_paper_style_batch_mixing_score.R"))

parse_args_09d <- function(defaults) {
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

params <- parse_args_09d(list(
  result_dir = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved",
  out_dir = "results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space",
  fig_dir = "figures/mouse_kidney",
  seed = 1L,
  pcs_list = "20,30,50",
  k = 30L
))

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""), "\n", sep = "")
}

parse_pcs_list <- function(x) {
  vals <- as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
  vals <- vals[!is.na(vals) & vals > 0]
  if (length(vals) == 0) stop("--pcs_list must contain at least one positive integer")
  sort(unique(vals))
}

required_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
  path
}

log_normalize_sparse <- function(counts) {
  counts <- as(counts, "dgCMatrix")
  counts@x <- rep(1, length(counts@x))
  depth <- Matrix::rowSums(counts)
  if (any(depth <= 0)) stop("Cannot log-normalize matrix with empty cells")
  scale_factor <- stats::median(as.numeric(depth))
  norm <- Diagonal(x = scale_factor / as.numeric(depth)) %*% counts
  norm <- as(norm, "dgCMatrix")
  norm@x <- log1p(norm@x)
  list(matrix = norm, depth = as.numeric(depth), scale_factor = scale_factor)
}

compute_pca_embedding <- function(counts, metadata, setting, max_pcs, out_csv) {
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package irlba is required for sparse PCA-logNorm embeddings.")
  }
  log_msg(setting, ": binarizing and log-normalizing sparse matrix")
  norm <- log_normalize_sparse(counts)
  nv <- min(max_pcs, nrow(norm$matrix) - 1L, ncol(norm$matrix) - 1L)
  if (nv < 2L) stop(setting, ": not enough rows/columns for PCA")
  log_msg(setting, ": PCA start with nv=", nv, " on dim=", paste(dim(norm$matrix), collapse = " x "))
  center_vec <- Matrix::colMeans(norm$matrix)
  set.seed(params$seed)
  pca <- irlba::prcomp_irlba(norm$matrix, n = nv, center = center_vec, scale. = FALSE)
  pcs <- as.data.frame(pca$x, check.names = FALSE)
  colnames(pcs) <- paste0("PC_", seq_len(ncol(pcs)))
  embedding <- metadata
  embedding[[paste0(setting, "_depth")]] <- norm$depth
  embedding <- cbind(embedding, pcs)
  write_csv(embedding, out_csv)
  log_msg(setting, ": saved PCA embedding: ", out_csv)
  list(
    embedding = embedding,
    n_features = ncol(counts),
    n_cells = nrow(counts),
    pca_dim = dim(pca$x),
    scale_factor = norm$scale_factor,
    out_csv = out_csv
  )
}

compute_paper_score_from_embedding <- function(embedding, setting, n_features, pcs) {
  pc_cols <- paste0("PC_", seq_len(pcs))
  if (!all(pc_cols %in% names(embedding))) {
    stop(setting, ": missing requested PC columns for PC1:", pcs)
  }
  coords <- as.matrix(embedding[, pc_cols, drop = FALSE])
  batch <- factor(as.character(embedding$batch))
  cell_type <- factor(as.character(embedding$cell_type))
  nn <- knn_indices(coords, params$k)
  cell_score <- vapply(seq_len(nrow(coords)), function(i) {
    mean(batch[nn[i, ]] != batch[i], na.rm = TRUE)
  }, numeric(1))
  observed <- mean(cell_score, na.rm = TRUE)
  expected <- expected_batch_mixing_score(cell_type, batch)
  data.frame(
    setting = setting,
    coordinate_space_type = "PCA-logNorm",
    pc_dims = pcs,
    n_cells = nrow(coords),
    n_features = n_features,
    observed_batch_mixing_score = observed,
    expected_batch_mixing_score = expected,
    normalized_batch_mixing_score = observed / expected,
    k = ncol(nn),
    coordinate_columns_used = paste(pc_cols, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

read_optional_previous_scores <- function() {
  rows <- list()
  lsi_path <- "results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_10000_lsi_space_paper_style_scores.csv"
  umap_path <- "results/mouse_kidney_figures/paper_style_batch_mixing_score/paper_style_batch_mixing_scores.csv"
  if (file.exists(lsi_path)) {
    lsi <- read_csv(lsi_path)
    lsi$setting <- gsub(", LSI space", "", lsi$space_name)
    lsi$pc_dims <- NA_integer_
    rows[[length(rows) + 1L]] <- data.frame(
      setting = lsi$setting,
      coordinate_space_type = lsi$coordinate_space_type,
      pc_dims = lsi$pc_dims,
      observed_batch_mixing_score = lsi$observed_batch_mixing_score,
      expected_batch_mixing_score = lsi$expected_batch_mixing_score,
      normalized_batch_mixing_score = lsi$normalized_batch_mixing_score,
      stringsAsFactors = FALSE
    )
  }
  if (file.exists(umap_path)) {
    umap <- read_csv(umap_path)
    umap <- umap[umap$space_name %in% c("P56 10000/10000 before", "P56 10000/10000 after"), , drop = FALSE]
    umap$setting <- paste0(umap$space_name, ", UMAP")
    umap$pc_dims <- NA_integer_
    rows[[length(rows) + 1L]] <- data.frame(
      setting = umap$setting,
      coordinate_space_type = umap$coordinate_space_type,
      pc_dims = umap$pc_dims,
      observed_batch_mixing_score = umap$observed_batch_mixing_score,
      expected_batch_mixing_score = umap$expected_batch_mixing_score,
      normalized_batch_mixing_score = umap$normalized_batch_mixing_score,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

plot_score_comparison <- function(pca_scores, out_base) {
  pca_plot <- data.frame(
    setting = paste0(pca_scores$setting, " PC1:", pca_scores$pc_dims),
    coordinate_space_type = pca_scores$coordinate_space_type,
    normalized_batch_mixing_score = pca_scores$normalized_batch_mixing_score,
    stringsAsFactors = FALSE
  )
  previous <- read_optional_previous_scores()
  if (nrow(previous) > 0) {
    previous_plot <- data.frame(
      setting = previous$setting,
      coordinate_space_type = previous$coordinate_space_type,
      normalized_batch_mixing_score = previous$normalized_batch_mixing_score,
      stringsAsFactors = FALSE
    )
    plot_df <- rbind(pca_plot, previous_plot)
  } else {
    plot_df <- pca_plot
  }
  plot_df$setting <- factor(plot_df$setting, levels = plot_df$setting)
  p <- ggplot(plot_df, aes(x = setting, y = normalized_batch_mixing_score, fill = coordinate_space_type)) +
    geom_col(width = 0.68) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#666666") +
    scale_fill_manual(values = c("PCA-logNorm" = "#2E86AB", "LSI-space" = "#00A651", "UMAP-space" = "#A000FF"), drop = FALSE) +
    labs(
      x = NULL,
      y = "Observed / expected batch mixing score",
      fill = "Coordinate space",
      title = "P56 10000/10000 paper-style batch mixing: PCA-logNorm vs LSI vs UMAP"
    ) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, ".png"), p, width = 12, height = 6.5, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 12, height = 6.5)
  plot_df
}

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- required_file(file.path(params$result_dir, "p56_counts_top_peaks_sparse.rds"), "before sparse counts RDS")
metadata_path <- required_file(file.path(params$result_dir, "p56_metadata.csv"), "P56 metadata CSV")
retained_path <- required_file(file.path(params$result_dir, "p56_retained_peak_indices.csv"), "retained peak CSV")
top_peak_path <- required_file(file.path(params$result_dir, "p56_top_peak_indices.csv"), "top peak CSV")

pcs_values <- parse_pcs_list(params$pcs_list)
max_pcs <- max(pcs_values)

log_msg("Loading sparse counts: ", counts_path)
counts <- readRDS(counts_path)
counts <- as(counts, "dgCMatrix")
metadata <- read_csv(metadata_path)
retained <- read_csv(retained_path)
top_peaks <- read_csv(top_peak_path)

if (nrow(metadata) != nrow(counts)) {
  stop("Metadata rows do not match counts rows: metadata=", nrow(metadata), ", counts=", nrow(counts))
}
if (!"batch" %in% names(metadata) && "sample" %in% names(metadata)) {
  metadata$batch <- metadata$sample
}
if (!all(c("batch", "cell_type") %in% names(metadata))) {
  stop("Metadata must contain batch/sample and cell_type columns")
}

if ("top_peak_col" %in% names(retained)) {
  retained_cols <- as.integer(retained$top_peak_col)
  mapping_note <- "After matrix was subset using p56_retained_peak_indices.csv column `top_peak_col`."
} else if (all(c("peak_index", "top_peak_col") %in% names(top_peaks)) && "peak_index" %in% names(retained)) {
  retained_cols <- match(retained$peak_index, top_peaks$peak_index)
  mapping_note <- "After matrix was subset by matching retained peak_index to p56_top_peak_indices.csv peak_index."
} else {
  stop("Could not map retained peaks to before matrix columns. Need retained top_peak_col or peak_index plus top peak mapping.")
}
if (anyNA(retained_cols) || any(retained_cols < 1L) || any(retained_cols > ncol(counts))) {
  stop("Retained peak column mapping contains NA or out-of-range values")
}
retained_cols <- sort(unique(retained_cols))
after_counts <- counts[, retained_cols, drop = FALSE]
log_msg("Before counts dim=", paste(dim(counts), collapse = " x "), "; after counts dim=", paste(dim(after_counts), collapse = " x "))

before_pca <- compute_pca_embedding(
  counts,
  metadata,
  "before",
  max_pcs,
  file.path(params$out_dir, "p56_10000_before_pca_logNorm_embedding.csv")
)
after_pca <- compute_pca_embedding(
  after_counts,
  metadata,
  "after",
  max_pcs,
  file.path(params$out_dir, "p56_10000_after_pca_logNorm_embedding.csv")
)

score_rows <- list()
for (pcs in pcs_values) {
  score_rows[[length(score_rows) + 1L]] <- compute_paper_score_from_embedding(before_pca$embedding, "before", before_pca$n_features, pcs)
  score_rows[[length(score_rows) + 1L]] <- compute_paper_score_from_embedding(after_pca$embedding, "after", after_pca$n_features, pcs)
}
pca_scores <- do.call(rbind, score_rows)
score_out <- file.path(params$out_dir, "p56_10000_pca_space_scores.csv")
write_csv(pca_scores, score_out)

plot_base <- file.path(params$fig_dir, "gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison")
plot_df <- plot_score_comparison(pca_scores, plot_base)

report_path <- file.path(params$out_dir, "p56_10000_pca_space_score_report.md")
report <- c(
  "# P56 10000/10000 PCA-logNorm paper-style batch mixing score",
  "",
  "## Goal",
  "",
  "This report computes a PCA-logNorm-space version of the PACS-paper-style normalized batch mixing score for the current P56 main setting.",
  "",
  "## Inputs",
  "",
  paste0("- result_dir: `", params$result_dir, "`"),
  paste0("- before sparse counts: `", counts_path, "`"),
  paste0("- metadata: `", metadata_path, "`"),
  paste0("- retained peaks: `", retained_path, "`"),
  "",
  "## Matrix loading and retained peak mapping",
  "",
  paste0("- before matrix dim: ", paste(dim(counts), collapse = " x ")),
  paste0("- after matrix dim: ", paste(dim(after_counts), collapse = " x ")),
  paste0("- retained peak columns: ", length(retained_cols)),
  paste0("- mapping: ", mapping_note),
  "",
  "## PCA-logNorm method",
  "",
  "- Counts were binarized.",
  "- Cells were library-size normalized to the median cell depth.",
  "- Values were log1p transformed.",
  "- Sparse truncated PCA was computed with `irlba::prcomp_irlba`.",
  "",
  "## PCA-space scores",
  "",
  "```text",
  paste(capture.output(print(pca_scores)), collapse = "\n"),
  "```",
  "",
  "## PCA vs LSI vs UMAP comparison",
  "",
  "```text",
  paste(capture.output(print(plot_df)), collapse = "\n"),
  "```",
  "",
  "## Interpretation",
  "",
  "- PCA-space score is closest to the PACS paper normalized PCA mixing metric.",
  "- LSI-space score is a scATAC-aware analogue.",
  "- UMAP-space score is a visualization-level approximation.",
  "- Do not claim exact numerical comparability to the PACS paper unless dataset, batch definition, feature set, and preprocessing match.",
  "",
  "## Output files",
  "",
  paste0("- `", score_out, "`"),
  paste0("- `", before_pca$out_csv, "`"),
  paste0("- `", after_pca$out_csv, "`"),
  paste0("- `", plot_base, ".png/pdf`")
)
writeLines(report, report_path)

log_msg("Saved PCA-space scores: ", score_out)
log_msg("Saved report: ", report_path)
log_msg("Saved plot: ", plot_base, ".png/pdf")
print(pca_scores)
