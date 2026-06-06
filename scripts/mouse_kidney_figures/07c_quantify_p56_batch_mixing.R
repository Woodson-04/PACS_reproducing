#!/usr/bin/env Rscript

# Quantify P56 before/after PACS batch-filtering UMAP results.
# This script only reads existing CSV outputs. It does not rerun PACS, matrix
# streaming, TF-IDF, LSI, or UMAP.

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
  k = 30L,
  seed = 1L,
  max_silhouette_n = 5000L
))

use_run_name <- nzchar(params$run_name) && params$run_name != "default"
if (!nzchar(params$result_dir)) {
  params$result_dir <- if (use_run_name) {
    file.path("results/mouse_kidney_figures", paste0("gse157079_p56_pacs_batch_filter_umap_", params$run_name))
  } else {
    "results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap"
  }
}
metric_prefix_name <- if (use_run_name) {
  paste0("gse157079_p56_batch_mixing_", params$run_name)
} else {
  "gse157079_p56_batch_mixing"
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

find_existing_file <- function(candidates, label) {
  paths <- file.path(params$result_dir, candidates)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) stop("Could not find ", label, " file. Tried: ", paste(paths, collapse = ", "))
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
      log_msg(label, " coordinate columns: ", paste(cols, collapse = ", "))
      return(cols)
    }
  }
  stop("Could not detect UMAP coordinate columns for ", label, ". Available columns: ", paste(names(df), collapse = ", "))
}

prepare_embedding <- function(df, cols, label) {
  if (!"cell_type" %in% names(df)) stop(label, " embedding lacks cell_type")
  if ("batch" %in% names(df)) {
    batch <- df$batch
  } else if ("sample" %in% names(df)) {
    batch <- df$sample
  } else {
    stop(label, " embedding lacks both batch and sample")
  }
  out <- data.frame(
    row_index = if ("row_index" %in% names(df)) df$row_index else seq_len(nrow(df)),
    cell_barcode = if ("cell_barcode" %in% names(df)) as.character(df$cell_barcode) else paste0("cell_", seq_len(nrow(df))),
    batch = factor(as.character(batch)),
    cell_type = factor(as.character(df$cell_type)),
    umap_1 = as.numeric(df[[cols[[1]]]]),
    umap_2 = as.numeric(df[[cols[[2]]]]),
    stringsAsFactors = FALSE
  )
  if (anyNA(out$umap_1) || anyNA(out$umap_2)) stop(label, " UMAP coordinates contain NA")
  out
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
    log_msg("kNN processed rows ", st, "-", en, " / ", n)
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
  log_msg("FNN/RANN unavailable; using chunked brute-force kNN")
  knn_indices_bruteforce(coords, k = k)
}

shannon_entropy <- function(x) {
  tab <- table(x)
  p <- as.numeric(tab) / sum(tab)
  -sum(p * log(p))
}

knn_metrics <- function(coords, batch, cell_type, k) {
  nn <- knn_indices(coords, k)
  batch <- factor(batch)
  cell_type <- factor(cell_type)
  n_batch <- length(levels(batch))
  batch_entropy <- numeric(nrow(coords))
  same_batch <- numeric(nrow(coords))
  same_celltype <- numeric(nrow(coords))
  for (i in seq_len(nrow(coords))) {
    nb <- nn[i, ]
    nb_batch <- batch[nb]
    nb_ct <- cell_type[nb]
    batch_entropy[i] <- if (n_batch > 1) shannon_entropy(nb_batch) / log(n_batch) else NA_real_
    same_batch[i] <- mean(nb_batch == batch[i], na.rm = TRUE)
    same_celltype[i] <- mean(nb_ct == cell_type[i], na.rm = TRUE)
  }
  list(
    normalized_batch_entropy = mean(batch_entropy, na.rm = TRUE),
    same_batch_fraction = mean(same_batch, na.rm = TRUE),
    same_celltype_fraction = mean(same_celltype, na.rm = TRUE),
    per_cell = data.frame(
      row_index = seq_len(nrow(coords)),
      normalized_batch_entropy = batch_entropy,
      same_batch_fraction = same_batch,
      same_celltype_fraction = same_celltype
    )
  )
}

batch_prediction_accuracy <- function(df, seed) {
  batch <- factor(df$batch)
  if (length(levels(batch)) != 2) return(NA_real_)
  set.seed(seed)
  idx <- sample(seq_len(nrow(df)))
  train_n <- floor(0.7 * nrow(df))
  train <- idx[seq_len(train_n)]
  test <- idx[(train_n + 1L):length(idx)]
  train_df <- data.frame(batch = batch[train], umap_1 = df$umap_1[train], umap_2 = df$umap_2[train])
  test_df <- data.frame(batch = batch[test], umap_1 = df$umap_1[test], umap_2 = df$umap_2[test])
  fit <- glm(batch ~ umap_1 + umap_2, data = train_df, family = binomial())
  prob <- predict(fit, newdata = test_df, type = "response")
  pred <- ifelse(prob >= 0.5, levels(batch)[2], levels(batch)[1])
  mean(pred == as.character(test_df$batch), na.rm = TRUE)
}

global_metrics <- function(df, stage, k, seed, max_silhouette_n) {
  coords <- as.matrix(df[, c("umap_1", "umap_2")])
  km <- knn_metrics(coords, df$batch, df$cell_type, k)
  data.frame(
    stage = stage,
    batch_silhouette = mean_silhouette(coords, df$batch, seed, max_silhouette_n),
    cell_type_silhouette = mean_silhouette(coords, df$cell_type, seed, max_silhouette_n),
    normalized_batch_entropy = km$normalized_batch_entropy,
    same_batch_fraction = km$same_batch_fraction,
    same_celltype_fraction = km$same_celltype_fraction,
    batch_prediction_accuracy = batch_prediction_accuracy(df, seed),
    stringsAsFactors = FALSE
  )
}

celltype_metrics_one <- function(before_df, after_df, k, seed, max_silhouette_n) {
  cts <- sort(intersect(unique(as.character(before_df$cell_type)), unique(as.character(after_df$cell_type))))
  rows <- list()
  for (ct in cts) {
    b <- before_df[before_df$cell_type == ct, , drop = FALSE]
    a <- after_df[after_df$cell_type == ct, , drop = FALSE]
    if (nrow(b) < 50 || length(unique(b$batch)) < 2 || length(unique(a$batch)) < 2) next
    bkm <- knn_metrics(as.matrix(b[, c("umap_1", "umap_2")]), b$batch, b$cell_type, min(k, nrow(b) - 1L))
    akm <- knn_metrics(as.matrix(a[, c("umap_1", "umap_2")]), a$batch, a$cell_type, min(k, nrow(a) - 1L))
    rows[[length(rows) + 1L]] <- data.frame(
      cell_type = ct,
      n_cells = nrow(b),
      batch_table = paste(names(table(b$batch)), as.integer(table(b$batch)), sep = "=", collapse = ";"),
      before_batch_entropy = bkm$normalized_batch_entropy,
      after_batch_entropy = akm$normalized_batch_entropy,
      before_same_batch_fraction = bkm$same_batch_fraction,
      after_same_batch_fraction = akm$same_batch_fraction,
      before_batch_silhouette = mean_silhouette(as.matrix(b[, c("umap_1", "umap_2")]), b$batch, seed, max_silhouette_n),
      after_batch_silhouette = mean_silhouette(as.matrix(a[, c("umap_1", "umap_2")]), a$batch, seed, max_silhouette_n),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

long_global_for_plot <- function(summary_df) {
  metrics <- c(
    batch_silhouette = "Batch silhouette",
    normalized_batch_entropy = "Normalized batch entropy",
    same_batch_fraction = "Same-batch fraction",
    same_celltype_fraction = "Same-celltype fraction",
    batch_prediction_accuracy = "Batch prediction accuracy"
  )
  out <- do.call(rbind, lapply(names(metrics), function(m) {
    data.frame(stage = summary_df$stage, metric = metrics[[m]], value = summary_df[[m]], stringsAsFactors = FALSE)
  }))
  out$stage <- factor(out$stage, levels = c("before", "after"))
  out
}

plot_global_metrics <- function(summary_df, out_base) {
  df <- long_global_for_plot(summary_df)
  p <- ggplot(df, aes(x = metric, y = value, fill = stage)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    scale_fill_manual(values = c(before = "#777777", after = "#2E86AB")) +
    labs(x = NULL, y = "Metric value", fill = "Stage", title = "P56 before/after PACS filtering batch-mixing metrics") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, "_metrics_barplot.png"), p, width = 9, height = 5.5, dpi = 300)
  ggsave(paste0(out_base, "_metrics_barplot.pdf"), p, width = 9, height = 5.5)
}

plot_celltype_metric <- function(df, before_col, after_col, ylab, title, out_base) {
  if (nrow(df) == 0) return(FALSE)
  long <- rbind(
    data.frame(cell_type = df$cell_type, stage = "before", value = df[[before_col]], stringsAsFactors = FALSE),
    data.frame(cell_type = df$cell_type, stage = "after", value = df[[after_col]], stringsAsFactors = FALSE)
  )
  long$stage <- factor(long$stage, levels = c("before", "after"))
  long$cell_type <- factor(long$cell_type, levels = df$cell_type[order(df[[before_col]], decreasing = TRUE)])
  p <- ggplot(long, aes(x = cell_type, y = value, fill = stage)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    scale_fill_manual(values = c(before = "#777777", after = "#2E86AB")) +
    labs(x = "Cell type", y = ylab, fill = "Stage", title = title) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))
  ggsave(paste0(out_base, ".png"), p, width = 9, height = 5.5, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 9, height = 5.5)
  TRUE
}

dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)
metric_prefix <- file.path(params$fig_dir, metric_prefix_name)
before_file <- find_existing_file(c("p56_before_lsi_umap_embedding.csv", "before_lsi_umap_embedding.csv"), "before embedding")
after_file <- find_existing_file(c("p56_after_lsi_umap_embedding.csv", "after_lsi_umap_embedding.csv"), "after embedding")
pacs_file <- file.path(params$result_dir, "p56_pacs_batch_peak_results.csv")
batch_summary_file <- file.path(params$result_dir, "p56_batch_peak_summary.csv")

log_msg("Reading embeddings")
before_raw <- read_csv(before_file)
after_raw <- read_csv(after_file)
before_cols <- detect_umap_cols(before_raw, "before", "before")
after_cols <- detect_umap_cols(after_raw, "after", "after")
before_df <- prepare_embedding(before_raw, before_cols, "before")
after_df <- prepare_embedding(after_raw, after_cols, "after")

if (nrow(before_df) != nrow(after_df)) stop("Before and after embeddings have different row counts")
if (!all(before_df$row_index == after_df$row_index)) {
  log_msg("Warning: row_index order differs between before and after; metrics are still computed by each file independently")
}

log_msg("Computing global before metrics")
before_metrics <- global_metrics(before_df, "before", params$k, params$seed, params$max_silhouette_n)
log_msg("Computing global after metrics")
after_metrics <- global_metrics(after_df, "after", params$k, params$seed, params$max_silhouette_n)
summary_df <- rbind(before_metrics, after_metrics)

log_msg("Computing cell-type-stratified metrics")
by_celltype <- celltype_metrics_one(before_df, after_df, params$k, params$seed, params$max_silhouette_n)

metrics_out <- file.path(params$result_dir, "p56_batch_mixing_metrics_summary.csv")
celltype_out <- file.path(params$result_dir, "p56_batch_mixing_by_celltype.csv")
write.csv(summary_df, metrics_out, row.names = FALSE)
write.csv(by_celltype, celltype_out, row.names = FALSE)

pacs_df <- if (file.exists(pacs_file)) read_csv(pacs_file) else data.frame()
batch_summary <- if (file.exists(batch_summary_file)) read_csv(batch_summary_file) else data.frame()
top20 <- data.frame()
if (nrow(pacs_df) > 0 && "p_value" %in% names(pacs_df)) {
  top20 <- pacs_df[order(pacs_df$p_value), , drop = FALSE]
  top20 <- head(top20, 20)
}

plot_global_metrics(summary_df, metric_prefix)
entropy_plot_ok <- plot_celltype_metric(
  by_celltype,
  "before_batch_entropy",
  "after_batch_entropy",
  "Normalized batch entropy",
  "P56 cell-type-level batch entropy before/after PACS filtering",
  paste0(metric_prefix, "_by_celltype_entropy")
)
samebatch_plot_ok <- plot_celltype_metric(
  by_celltype,
  "before_same_batch_fraction",
  "after_same_batch_fraction",
  "Same-batch fraction",
  "P56 cell-type-level same-batch fraction before/after PACS filtering",
  paste0(metric_prefix, "_by_celltype_samebatch")
)

metric_delta <- setNames(summary_df[summary_df$stage == "after", -1] - summary_df[summary_df$stage == "before", -1], names(summary_df)[-1])
improved_batch <- isTRUE(metric_delta$batch_silhouette < 0) &&
  isTRUE(metric_delta$normalized_batch_entropy > 0) &&
  isTRUE(metric_delta$same_batch_fraction < 0)
celltype_drop <- isTRUE(metric_delta$same_celltype_fraction < -0.10) ||
  isTRUE(metric_delta$cell_type_silhouette < -0.10)

interpretation <- c()
if (improved_batch) {
  interpretation <- c(interpretation, "Batch silhouette decreased, normalized batch entropy increased, and same-batch fraction decreased. These metrics support improved batch mixing after PACS filtering.")
} else {
  interpretation <- c(interpretation, "The batch-mixing metrics do not all move in the expected direction. This suggests improvement is weak, mixed, or parameter-dependent.")
}
if (celltype_drop) {
  interpretation <- c(interpretation, "Cell type structure metrics dropped noticeably. This may indicate possible overcorrection or loss of biological structure and should be checked visually.")
} else {
  interpretation <- c(interpretation, "Cell type structure appears reasonably preserved by the same-celltype and cell type silhouette metrics.")
}
interpretation <- c(interpretation, "Residual batch separation may remain; this quantification should be interpreted as an initial P56-only top-peak PACS filtering assessment, not a final full correction.")

report_path <- file.path(params$result_dir, "p56_batch_mixing_quantification_report.md")
report <- c(
  "# P56 Batch Mixing Quantification Report",
  "",
  "This report quantifies existing before/after P56 PACS batch-filtering UMAP embeddings. It did not rerun PACS, MatrixMarket streaming, TF-IDF, LSI, or UMAP.",
  "",
  "## Inputs",
  "",
  paste0("- before embedding: `", before_file, "`"),
  paste0("- after embedding: `", after_file, "`"),
  paste0("- before coordinate columns: ", paste(before_cols, collapse = ", ")),
  paste0("- after coordinate columns: ", paste(after_cols, collapse = ", ")),
  paste0("- k for kNN metrics: ", params$k),
  paste0("- max_silhouette_n: ", params$max_silhouette_n),
  "",
  "## Global Metrics",
  "",
  "```text",
  paste(capture.output(print(summary_df)), collapse = "\n"),
  "```",
  "",
  "Metric interpretation:",
  "",
  "- Lower batch silhouette indicates better batch mixing.",
  "- Higher normalized batch entropy indicates better batch mixing.",
  "- Lower same-batch fraction indicates better batch mixing.",
  "- Higher same-celltype fraction and cell type silhouette suggest preserved biological structure.",
  "- Lower batch prediction accuracy suggests weaker batch separability in UMAP space.",
  "",
  "## Cell-Type-Stratified Metrics",
  "",
  "```text",
  paste(capture.output(print(by_celltype)), collapse = "\n"),
  "```",
  "",
  "## PACS Peak Summary",
  "",
  if (nrow(batch_summary) > 0) paste(capture.output(print(batch_summary)), collapse = "\n") else "No p56_batch_peak_summary.csv found.",
  "",
  "## PACS p-value/FDR Summary",
  "",
  "```text",
  if (nrow(pacs_df) > 0 && all(c("p_value", "fdr") %in% names(pacs_df))) {
    paste(capture.output(print(summary(pacs_df[, c("p_value", "fdr")]))), collapse = "\n")
  } else {
    "PACS p_value/fdr columns unavailable."
  },
  "```",
  "",
  "## Top 20 Most Significant Batch Peaks",
  "",
  "```text",
  if (nrow(top20) > 0) paste(capture.output(print(top20)), collapse = "\n") else "Top peak table unavailable.",
  "```",
  "",
  "## Output Files",
  "",
  paste0("- `", metrics_out, "`"),
  paste0("- `", celltype_out, "`"),
  paste0("- `", metric_prefix, "_metrics_barplot.png/pdf`"),
  paste0("- `", metric_prefix, "_by_celltype_entropy.png/pdf`"),
  paste0("- `", metric_prefix, "_by_celltype_samebatch.png/pdf`"),
  "",
  "## Interpretation",
  "",
  interpretation
)
writeLines(report, report_path)

log_msg("Saved metrics summary: ", metrics_out)
log_msg("Saved cell-type metrics: ", celltype_out)
log_msg("Saved quantification report: ", report_path)
cat("before_batch_silhouette=", before_metrics$batch_silhouette, "\n", sep = "")
cat("after_batch_silhouette=", after_metrics$batch_silhouette, "\n", sep = "")
cat("before_normalized_batch_entropy=", before_metrics$normalized_batch_entropy, "\n", sep = "")
cat("after_normalized_batch_entropy=", after_metrics$normalized_batch_entropy, "\n", sep = "")
cat("before_same_batch_fraction=", before_metrics$same_batch_fraction, "\n", sep = "")
cat("after_same_batch_fraction=", after_metrics$same_batch_fraction, "\n", sep = "")
cat("before_same_celltype_fraction=", before_metrics$same_celltype_fraction, "\n", sep = "")
cat("after_same_celltype_fraction=", after_metrics$same_celltype_fraction, "\n", sep = "")
cat("before_batch_prediction_accuracy=", before_metrics$batch_prediction_accuracy, "\n", sep = "")
cat("after_batch_prediction_accuracy=", after_metrics$batch_prediction_accuracy, "\n", sep = "")
cat("entropy_plot_generated=", entropy_plot_ok, "\n", sep = "")
cat("samebatch_plot_generated=", samebatch_plot_ok, "\n", sep = "")
