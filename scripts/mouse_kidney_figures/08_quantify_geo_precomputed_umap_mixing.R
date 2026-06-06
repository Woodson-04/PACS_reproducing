#!/usr/bin/env Rscript

# Quantify the GEO-provided precomputed GSE157079 UMAP coordinates as a
# reference embedding. This script does not rerun PACS, MatrixMarket streaming,
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
  metadata_csv = "results/mouse_kidney_figures/gse157079_metadata_merged.csv",
  out_dir = "results/mouse_kidney_figures/geo_precomputed_umap_quantification",
  fig_dir = "figures/mouse_kidney",
  k = 30L,
  seed = 1L,
  max_silhouette_n = 5000L
))

read_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, header = TRUE, showProgress = FALSE, check.names = FALSE)))
  }
  read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""), "\n", sep = "")
}

require_columns <- function(df, cols, label) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(label, " missing required columns: ", paste(missing, collapse = ", "),
         ". Available columns: ", paste(names(df), collapse = ", "))
  }
}

make_age_group <- function(sample) {
  out <- ifelse(grepl("^P0", sample), "P0",
    ifelse(grepl("^P21", sample), "P21",
      ifelse(grepl("^P56", sample), "P56", NA_character_)))
  factor(out, levels = c("P0", "P21", "P56"))
}

mean_silhouette <- function(coords, labels, seed, max_n = 5000L) {
  labels <- factor(labels)
  keep <- !is.na(labels)
  coords <- coords[keep, , drop = FALSE]
  labels <- labels[keep]
  if (length(unique(labels)) < 2 || nrow(coords) < 3) return(NA_real_)
  if (nrow(coords) > max_n) {
    set.seed(seed)
    idx <- sort(sample(seq_len(nrow(coords)), max_n))
    coords <- coords[idx, , drop = FALSE]
    labels <- labels[idx]
  }
  if (length(unique(labels)) < 2) return(NA_real_)
  d <- dist(coords)
  sil <- cluster::silhouette(as.integer(labels), d)
  mean(sil[, "sil_width"], na.rm = TRUE)
}

knn_indices_bruteforce <- function(coords, k, chunk_size = 500L) {
  n <- nrow(coords)
  k <- min(k, n - 1L)
  out <- matrix(NA_integer_, nrow = n, ncol = k)
  x2 <- rowSums(coords^2)
  for (st in seq(1L, n, by = chunk_size)) {
    en <- min(n, st + chunk_size - 1L)
    block <- coords[st:en, , drop = FALSE]
    d2 <- outer(rowSums(block^2), x2, "+") - 2 * (block %*% t(coords))
    for (ii in seq_len(nrow(d2))) {
      global_i <- st + ii - 1L
      d2[ii, global_i] <- Inf
      out[global_i, ] <- order(d2[ii, ], decreasing = FALSE)[seq_len(k)]
    }
    log_msg("brute-force kNN processed rows ", st, "-", en, " / ", n)
  }
  out
}

knn_indices <- function(coords, k) {
  if (requireNamespace("FNN", quietly = TRUE)) {
    log_msg("Using FNN::get.knn for kNN metrics")
    return(FNN::get.knn(coords, k = k)$nn.index)
  }
  if (requireNamespace("RANN", quietly = TRUE)) {
    log_msg("Using RANN::nn2 for kNN metrics")
    idx <- RANN::nn2(coords, coords, k = k + 1L)$nn.idx
    return(idx[, -1L, drop = FALSE])
  }
  if (nrow(coords) > 10000) {
    stop("FNN or RANN is required for kNN metrics on more than 10000 cells.")
  }
  log_msg("FNN/RANN unavailable; using chunked brute-force kNN")
  knn_indices_bruteforce(coords, k = k)
}

shannon_entropy <- function(x) {
  tab <- table(x)
  p <- as.numeric(tab) / sum(tab)
  -sum(p * log(p))
}

knn_label_metrics <- function(coords, label, k, prefix) {
  nn <- knn_indices(coords, k)
  label <- factor(label)
  n_label <- length(levels(label))
  entropy <- numeric(nrow(coords))
  same_fraction <- numeric(nrow(coords))
  for (i in seq_len(nrow(coords))) {
    nb <- nn[i, ]
    nb_label <- label[nb]
    entropy[i] <- if (n_label > 1) shannon_entropy(nb_label) / log(n_label) else NA_real_
    same_fraction[i] <- mean(nb_label == label[i], na.rm = TRUE)
  }
  out <- list()
  out[[paste0("normalized_", prefix, "_entropy")]] <- mean(entropy, na.rm = TRUE)
  out[[paste0("same_", prefix, "_fraction")]] <- mean(same_fraction, na.rm = TRUE)
  out
}

prediction_accuracy <- function(df, label_col, seed) {
  label <- factor(df[[label_col]])
  if (length(levels(label)) < 2) return(NA_real_)
  set.seed(seed)
  idx <- sample(seq_len(nrow(df)))
  train_n <- floor(0.7 * nrow(df))
  train <- idx[seq_len(train_n)]
  test <- idx[(train_n + 1L):length(idx)]
  train_df <- data.frame(label = label[train], umap_1 = df$umap_1[train], umap_2 = df$umap_2[train])
  test_df <- data.frame(label = label[test], umap_1 = df$umap_1[test], umap_2 = df$umap_2[test])
  if (length(levels(label)) == 2) {
    fit <- glm(label ~ umap_1 + umap_2, data = train_df, family = binomial())
    prob <- predict(fit, newdata = test_df, type = "response")
    pred <- ifelse(prob >= 0.5, levels(label)[2], levels(label)[1])
    return(mean(pred == as.character(test_df$label), na.rm = TRUE))
  }
  if (requireNamespace("nnet", quietly = TRUE)) {
    fit <- nnet::multinom(label ~ umap_1 + umap_2, data = train_df, trace = FALSE, maxit = 200)
    pred <- predict(fit, newdata = test_df)
    return(mean(as.character(pred) == as.character(test_df$label), na.rm = TRUE))
  }
  centroids <- aggregate(cbind(umap_1, umap_2) ~ label, data = train_df, FUN = mean)
  pred <- vapply(seq_len(nrow(test_df)), function(i) {
    d2 <- (centroids$umap_1 - test_df$umap_1[i])^2 + (centroids$umap_2 - test_df$umap_2[i])^2
    as.character(centroids$label[which.min(d2)])
  }, character(1))
  mean(pred == as.character(test_df$label), na.rm = TRUE)
}

global_metrics <- function(df, stage, batch_col, bio_col, k, seed, max_silhouette_n) {
  coords <- as.matrix(df[, c("umap_1", "umap_2")])
  batch_knn <- knn_label_metrics(coords, df[[batch_col]], k, "batch")
  bio_knn <- knn_label_metrics(coords, df[[bio_col]], k, "celltype")
  data.frame(
    stage = stage,
    batch_silhouette = mean_silhouette(coords, df[[batch_col]], seed, max_silhouette_n),
    cell_type_silhouette = mean_silhouette(coords, df[[bio_col]], seed, max_silhouette_n),
    normalized_batch_entropy = batch_knn$normalized_batch_entropy,
    same_batch_fraction = batch_knn$same_batch_fraction,
    same_celltype_fraction = bio_knn$same_celltype_fraction,
    batch_prediction_accuracy = prediction_accuracy(df, batch_col, seed),
    stringsAsFactors = FALSE
  )
}

age_metrics <- function(df, k, seed, max_silhouette_n) {
  coords <- as.matrix(df[, c("umap_1", "umap_2")])
  age_knn <- knn_label_metrics(coords, df$age_group, k, "age")
  data.frame(
    stage = "GEO precomputed UMAP global age reference",
    age_group_silhouette = mean_silhouette(coords, df$age_group, seed, max_silhouette_n),
    same_age_fraction = age_knn$same_age_fraction,
    normalized_age_entropy = age_knn$normalized_age_entropy,
    age_group_prediction_accuracy = prediction_accuracy(df, "age_group", seed),
    stringsAsFactors = FALSE
  )
}

celltype_p56_metrics <- function(df, k, seed, max_silhouette_n) {
  cts <- sort(unique(as.character(df$cell_type)))
  rows <- list()
  for (ct in cts) {
    sub <- df[df$cell_type == ct, , drop = FALSE]
    if (nrow(sub) < 50 || length(unique(sub$sample)) < 2) next
    coords <- as.matrix(sub[, c("umap_1", "umap_2")])
    km <- knn_label_metrics(coords, sub$sample, min(k, nrow(sub) - 1L), "batch")
    rows[[length(rows) + 1L]] <- data.frame(
      cell_type = ct,
      n_cells = nrow(sub),
      batch_table = paste(names(table(sub$sample)), as.integer(table(sub$sample)), sep = "=", collapse = ";"),
      batch_silhouette = mean_silhouette(coords, sub$sample, seed, max_silhouette_n),
      normalized_batch_entropy = km$normalized_batch_entropy,
      same_batch_fraction = km$same_batch_fraction,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

plot_global_metrics <- function(df, metrics, title, out_base) {
  long <- do.call(rbind, lapply(names(metrics), function(col) {
    data.frame(metric = metrics[[col]], value = df[[col]], stringsAsFactors = FALSE)
  }))
  p <- ggplot(long, aes(x = metric, y = value)) +
    geom_col(width = 0.68, fill = "#2E86AB") +
    labs(x = NULL, y = "Metric value", title = title) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, ".png"), p, width = 9, height = 5.5, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 9, height = 5.5)
}

plot_celltype_metric <- function(df, value_col, ylab, title, out_base) {
  if (nrow(df) == 0) return(FALSE)
  df$cell_type <- factor(df$cell_type, levels = df$cell_type[order(df[[value_col]], decreasing = TRUE)])
  p <- ggplot(df, aes(x = cell_type, y = .data[[value_col]])) +
    geom_col(width = 0.68, fill = "#2E86AB") +
    labs(x = "Cell type", y = ylab, title = title) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, ".png"), p, width = 9, height = 5.5, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 9, height = 5.5)
  TRUE
}

read_pacs_after <- function(path, setting) {
  if (!file.exists(path)) return(NULL)
  df <- read_csv(path)
  after <- df[df$stage == "after", , drop = FALSE]
  if (nrow(after) == 0) return(NULL)
  data.frame(
    setting = setting,
    batch_silhouette = after$batch_silhouette,
    normalized_batch_entropy = after$normalized_batch_entropy,
    same_batch_fraction = after$same_batch_fraction,
    batch_prediction_accuracy = after$batch_prediction_accuracy,
    cell_type_silhouette = after$cell_type_silhouette,
    same_celltype_fraction = after$same_celltype_fraction,
    stringsAsFactors = FALSE
  )
}

plot_comparison <- function(df, out_base) {
  metric_names <- c(
    batch_silhouette = "Batch silhouette",
    normalized_batch_entropy = "Normalized batch entropy",
    same_batch_fraction = "Same-batch fraction",
    batch_prediction_accuracy = "Batch prediction accuracy",
    cell_type_silhouette = "Cell type silhouette",
    same_celltype_fraction = "Same-celltype fraction"
  )
  long <- do.call(rbind, lapply(names(metric_names), function(col) {
    data.frame(setting = df$setting, metric = metric_names[[col]], value = df[[col]], stringsAsFactors = FALSE)
  }))
  long$setting <- factor(long$setting, levels = df$setting)
  p <- ggplot(long, aes(x = metric, y = value, fill = setting)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.68) +
    scale_fill_manual(values = c("#6B7280", "#2E86AB", "#00A651", "#A000FF"), drop = FALSE) +
    labs(x = NULL, y = "Metric value", fill = "Embedding", title = "GEO precomputed UMAP vs P56 PACS-filtered UMAP metrics") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, ".png"), p, width = 11, height = 6, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 11, height = 6)
}

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading merged metadata")
df <- read_csv(params$metadata_csv)
require_columns(df, c("row_index", "sample", "cell_type", "umap_1", "umap_2"), "metadata")
df$sample <- factor(as.character(df$sample))
df$cell_type <- factor(as.character(df$cell_type))
df$age_group <- make_age_group(as.character(df$sample))
df$umap_1 <- as.numeric(df$umap_1)
df$umap_2 <- as.numeric(df$umap_2)
if (anyNA(df$umap_1) || anyNA(df$umap_2)) stop("GEO UMAP coordinates contain NA")
if (anyNA(df$age_group)) stop("Could not map all samples to age_group")

log_msg("Computing global GEO UMAP metrics")
global_geo <- global_metrics(
  df,
  "GEO precomputed UMAP global",
  "sample",
  "cell_type",
  params$k,
  params$seed,
  params$max_silhouette_n
)
age_geo <- age_metrics(df, params$k, params$seed, params$max_silhouette_n)

log_msg("Computing P56-only GEO UMAP metrics")
p56 <- df[df$sample %in% c("P56_batch1", "P56_batch2"), , drop = FALSE]
p56$sample <- droplevels(p56$sample)
p56$cell_type <- droplevels(p56$cell_type)
p56_geo <- global_metrics(
  p56,
  "GEO precomputed UMAP P56 subset",
  "sample",
  "cell_type",
  params$k,
  params$seed,
  params$max_silhouette_n
)
by_celltype <- celltype_p56_metrics(p56, params$k, params$seed, params$max_silhouette_n)

global_out <- file.path(params$out_dir, "geo_precomputed_umap_global_metrics.csv")
p56_out <- file.path(params$out_dir, "geo_precomputed_umap_p56_metrics.csv")
by_celltype_out <- file.path(params$out_dir, "geo_precomputed_umap_p56_by_celltype_metrics.csv")
write.csv(cbind(global_geo, age_geo[, -1, drop = FALSE]), global_out, row.names = FALSE)
write.csv(p56_geo, p56_out, row.names = FALSE)
write.csv(by_celltype, by_celltype_out, row.names = FALSE)

plot_global_metrics(
  cbind(global_geo, age_geo[, -1, drop = FALSE]),
  c(
    batch_silhouette = "Sample silhouette",
    normalized_batch_entropy = "Normalized sample entropy",
    same_batch_fraction = "Same-sample fraction",
    batch_prediction_accuracy = "Sample prediction accuracy",
    cell_type_silhouette = "Cell type silhouette",
    same_celltype_fraction = "Same-celltype fraction",
    age_group_silhouette = "Age silhouette",
    same_age_fraction = "Same-age fraction",
    age_group_prediction_accuracy = "Age prediction accuracy"
  ),
  "GEO precomputed UMAP global reference metrics",
  file.path(params$fig_dir, "gse157079_geo_precomputed_umap_mixing_metrics")
)

plot_global_metrics(
  p56_geo,
  c(
    batch_silhouette = "P56 batch silhouette",
    normalized_batch_entropy = "Normalized batch entropy",
    same_batch_fraction = "Same-batch fraction",
    batch_prediction_accuracy = "Batch prediction accuracy",
    cell_type_silhouette = "Cell type silhouette",
    same_celltype_fraction = "Same-celltype fraction"
  ),
  "GEO precomputed UMAP P56 reference metrics",
  file.path(params$fig_dir, "gse157079_geo_precomputed_umap_p56_mixing_metrics")
)

entropy_plot_ok <- plot_celltype_metric(
  by_celltype,
  "normalized_batch_entropy",
  "Normalized batch entropy",
  "GEO precomputed UMAP P56 by-celltype batch entropy",
  file.path(params$fig_dir, "gse157079_geo_precomputed_umap_p56_by_celltype_entropy")
)
samebatch_plot_ok <- plot_celltype_metric(
  by_celltype,
  "same_batch_fraction",
  "Same-batch fraction",
  "GEO precomputed UMAP P56 by-celltype same-batch fraction",
  file.path(params$fig_dir, "gse157079_geo_precomputed_umap_p56_by_celltype_samebatch")
)

pacs_paths <- list(
  "P56 PACS 5000/5000 after" = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap/p56_batch_mixing_metrics_summary.csv",
  "P56 PACS 10000/10000 after" = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005/p56_batch_mixing_metrics_summary.csv",
  "P56 PACS 20000/10000 after" = "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top20000_test10000_fdr005/p56_batch_mixing_metrics_summary.csv"
)
comparison_rows <- list(data.frame(
  setting = "GEO precomputed UMAP, P56 subset",
  batch_silhouette = p56_geo$batch_silhouette,
  normalized_batch_entropy = p56_geo$normalized_batch_entropy,
  same_batch_fraction = p56_geo$same_batch_fraction,
  batch_prediction_accuracy = p56_geo$batch_prediction_accuracy,
  cell_type_silhouette = p56_geo$cell_type_silhouette,
  same_celltype_fraction = p56_geo$same_celltype_fraction,
  stringsAsFactors = FALSE
))
for (nm in names(pacs_paths)) {
  row <- read_pacs_after(pacs_paths[[nm]], nm)
  if (!is.null(row)) comparison_rows[[length(comparison_rows) + 1L]] <- row
}
comparison <- do.call(rbind, comparison_rows)
comparison_out <- file.path(params$out_dir, "geo_vs_p56_pacs_umap_metrics_comparison.csv")
write.csv(comparison, comparison_out, row.names = FALSE)
plot_comparison(comparison, file.path(params$fig_dir, "gse157079_geo_vs_p56_pacs_umap_metrics_comparison"))

report_path <- file.path(params$out_dir, "geo_precomputed_umap_quantification_report.md")
report <- c(
  "# GEO precomputed UMAP quantitative reference",
  "",
  "## Goal",
  "",
  "This report evaluates the GEO-provided precomputed GSE157079 UMAP coordinates as a public reference embedding. It should not be interpreted as the PACS paper filtered UMAP unless directly documented.",
  "",
  "## Global GSE157079 UMAP metrics",
  "",
  "The global analysis uses sample as the sample/batch label and cell_type as the biological label. Age/developmental structure is also summarized separately because P0, P21, and P56 differences are biological rather than purely technical batch effects.",
  "",
  "```text",
  paste(capture.output(print(cbind(global_geo, age_geo[, -1, drop = FALSE]))), collapse = "\n"),
  "```",
  "",
  "## P56-only GEO UMAP metrics",
  "",
  "This is the most directly comparable GEO reference for the P56_batch1 vs P56_batch2 PACS-filtered UMAP analyses.",
  "",
  "```text",
  paste(capture.output(print(p56_geo)), collapse = "\n"),
  "```",
  "",
  "## P56 cell-type-stratified GEO metrics",
  "",
  "```text",
  paste(capture.output(print(by_celltype)), collapse = "\n"),
  "```",
  "",
  "## Comparison with P56 PACS-filtered UMAP",
  "",
  "```text",
  paste(capture.output(print(comparison)), collapse = "\n"),
  "```",
  "",
  "Lower batch silhouette, lower same-batch fraction, and lower batch prediction accuracy indicate better batch mixing. Higher normalized batch entropy indicates better batch mixing. Higher cell type silhouette and same-celltype fraction indicate stronger biological cell type structure.",
  "",
  "## Interpretation caveat",
  "",
  "The GEO UMAP may reflect original atlas processing choices, possibly including integration or correction steps, but we should not claim it is the PACS paper batch-filtered UMAP unless directly documented.",
  "",
  "## Use in presentation",
  "",
  "我们将 GEO 预计算 UMAP 作为公开参考 embedding，用同一套 batch mixing 指标评估其混合程度，并与 P56 PACS-filtered UMAP 对照。",
  "",
  "## Output files",
  "",
  paste0("- `", global_out, "`"),
  paste0("- `", p56_out, "`"),
  paste0("- `", by_celltype_out, "`"),
  paste0("- `", comparison_out, "`"),
  paste0("- `", file.path(params$fig_dir, "gse157079_geo_precomputed_umap_mixing_metrics.png/pdf"), "`"),
  paste0("- `", file.path(params$fig_dir, "gse157079_geo_precomputed_umap_p56_mixing_metrics.png/pdf"), "`"),
  paste0("- `", file.path(params$fig_dir, "gse157079_geo_precomputed_umap_p56_by_celltype_entropy.png/pdf"), "`"),
  paste0("- `", file.path(params$fig_dir, "gse157079_geo_precomputed_umap_p56_by_celltype_samebatch.png/pdf"), "`"),
  paste0("- `", file.path(params$fig_dir, "gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png/pdf"), "`")
)
writeLines(report, report_path)

log_msg("Saved global GEO metrics: ", global_out)
log_msg("Saved P56 GEO metrics: ", p56_out)
log_msg("Saved P56 by-celltype GEO metrics: ", by_celltype_out)
log_msg("Saved GEO vs PACS comparison: ", comparison_out)
log_msg("Saved report: ", report_path)
cat("geo_p56_batch_silhouette=", p56_geo$batch_silhouette, "\n", sep = "")
cat("geo_p56_normalized_batch_entropy=", p56_geo$normalized_batch_entropy, "\n", sep = "")
cat("geo_p56_same_batch_fraction=", p56_geo$same_batch_fraction, "\n", sep = "")
cat("geo_p56_batch_prediction_accuracy=", p56_geo$batch_prediction_accuracy, "\n", sep = "")
cat("geo_p56_cell_type_silhouette=", p56_geo$cell_type_silhouette, "\n", sep = "")
cat("geo_p56_same_celltype_fraction=", p56_geo$same_celltype_fraction, "\n", sep = "")
cat("entropy_plot_generated=", entropy_plot_ok, "\n", sep = "")
cat("samebatch_plot_generated=", samebatch_plot_ok, "\n", sep = "")
