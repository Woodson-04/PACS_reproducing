#!/usr/bin/env Rscript

# Compute a PACS-paper-style normalized batch mixing score.
#
# Cell-level observed score: proportion of k nearest neighbors with different
# batch identities.
# Data-level observed score: mean cell-level score.
# Expected score: computed from the cell_type-by-batch count matrix.
# Normalized score: observed / expected.
#
# The PACS paper reports this score in PCA space. For the current scATAC
# pipeline, LSI space is the closest analogue when LSI coordinates are present.
# UMAP-space scores are reported as visualization-space approximations.

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

read_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, header = TRUE, showProgress = FALSE, check.names = FALSE)))
  }
  read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
}

detect_coord_cols <- function(df, coord_prefix = "") {
  names_df <- names(df)
  if (nzchar(coord_prefix)) {
    pref_lsi <- paste0(coord_prefix, "_lsi_", 2:30)
    pref_lsi_caps <- paste0(coord_prefix, "_LSI_", 2:30)
    pref_umap <- paste0(coord_prefix, "_umap_", 1:2)
    pref_umap_caps <- paste0(coord_prefix, "_UMAP_", 1:2)
    if (all(pref_lsi %in% names_df)) return(list(cols = pref_lsi, type = "LSI-space"))
    if (all(pref_lsi_caps %in% names_df)) return(list(cols = pref_lsi_caps, type = "LSI-space"))
    if (all(pref_umap %in% names_df)) return(list(cols = pref_umap, type = "UMAP-space"))
    if (all(pref_umap_caps %in% names_df)) return(list(cols = pref_umap_caps, type = "UMAP-space"))
  }

  lsi_caps <- paste0("LSI_", 2:30)
  lsi_lower <- paste0("lsi_", 2:30)
  before_lsi <- paste0("before_lsi_", 2:30)
  after_lsi <- paste0("after_lsi_", 2:30)
  if (all(lsi_caps %in% names_df)) return(list(cols = lsi_caps, type = "LSI-space"))
  if (all(lsi_lower %in% names_df)) return(list(cols = lsi_lower, type = "LSI-space"))
  if (all(before_lsi %in% names_df)) return(list(cols = before_lsi, type = "LSI-space"))
  if (all(after_lsi %in% names_df)) return(list(cols = after_lsi, type = "LSI-space"))

  candidates <- list(
    c("before_umap_1", "before_umap_2"),
    c("after_umap_1", "after_umap_2"),
    c("lsi_umap_1", "lsi_umap_2"),
    c("umap_1", "umap_2"),
    c("UMAP_1", "UMAP_2")
  )
  for (cols in candidates) {
    if (all(cols %in% names_df)) return(list(cols = cols, type = "UMAP-space"))
  }
  stop("Could not detect coordinate columns. Available columns: ", paste(names_df, collapse = ", "))
}

knn_indices <- function(coords, k) {
  k <- min(k, nrow(coords) - 1L)
  if (k < 1L) stop("Need at least two cells for kNN")
  if (requireNamespace("FNN", quietly = TRUE)) {
    return(FNN::get.knn(coords, k = k)$nn.index)
  }
  if (requireNamespace("RANN", quietly = TRUE)) {
    idx <- RANN::nn2(coords, coords, k = k + 1L)$nn.idx
    return(idx[, -1L, drop = FALSE])
  }
  stop("Package FNN or RANN is required for nearest-neighbor batch mixing scores.")
}

expected_batch_mixing_score <- function(cell_type, batch) {
  tab <- table(cell_type, batch)
  total <- sum(tab)
  expected_sum <- 0
  for (i in seq_len(nrow(tab))) {
    row_total <- sum(tab[i, ])
    if (row_total == 0) next
    for (j in seq_len(ncol(tab))) {
      m_ij <- tab[i, j]
      expected_ij <- (row_total - m_ij) / row_total
      expected_sum <- expected_sum + m_ij * expected_ij
    }
  }
  as.numeric(expected_sum / total)
}

compute_score_df <- function(df,
                             space_name,
                             batch_col = "batch",
                             celltype_col = "cell_type",
                             k = 30L,
                             coord_prefix = "",
                             seed = 1L) {
  if (!batch_col %in% names(df) && "sample" %in% names(df)) {
    batch_col <- "sample"
  }
  missing <- setdiff(c(batch_col, celltype_col), names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "),
         ". Available columns: ", paste(names(df), collapse = ", "))
  }
  coord <- detect_coord_cols(df, coord_prefix = coord_prefix)
  coords <- as.matrix(df[, coord$cols, drop = FALSE])
  coords <- apply(coords, 2, as.numeric)
  if (anyNA(coords)) stop("Coordinate columns contain NA after numeric conversion")

  batch <- factor(as.character(df[[batch_col]]))
  cell_type <- factor(as.character(df[[celltype_col]]))
  keep <- !is.na(batch) & !is.na(cell_type) & rowSums(is.na(coords)) == 0
  coords <- coords[keep, , drop = FALSE]
  batch <- droplevels(batch[keep])
  cell_type <- droplevels(cell_type[keep])
  if (length(unique(batch)) < 2) stop("Need at least two batch identities")
  if (length(unique(cell_type)) < 1) stop("Need at least one cell type")

  set.seed(seed)
  nn <- knn_indices(coords, k)
  k_eff <- ncol(nn)
  cell_score <- vapply(seq_len(nrow(coords)), function(i) {
    mean(batch[nn[i, ]] != batch[i], na.rm = TRUE)
  }, numeric(1))
  observed <- mean(cell_score, na.rm = TRUE)
  expected <- expected_batch_mixing_score(cell_type, batch)
  normalized <- observed / expected

  summary <- data.frame(
    space_name = space_name,
    coordinate_space_type = coord$type,
    n_cells = nrow(coords),
    k = k_eff,
    n_batches = length(levels(batch)),
    n_cell_types = length(levels(cell_type)),
    observed_batch_mixing_score = observed,
    expected_batch_mixing_score = expected,
    normalized_batch_mixing_score = normalized,
    coordinate_columns_used = paste(coord$cols, collapse = ";"),
    stringsAsFactors = FALSE
  )
  per_cell <- data.frame(
    row_index = if ("row_index" %in% names(df)) df$row_index[keep] else seq_len(nrow(coords)),
    batch = as.character(batch),
    cell_type = as.character(cell_type),
    cell_observed_mixing_score = cell_score,
    stringsAsFactors = FALSE
  )
  list(summary = summary, per_cell = per_cell)
}

write_score_outputs <- function(result, out_csv) {
  write_csv(result$summary, out_csv)
  per_cell_csv <- sub("\\.csv$", "_per_cell.csv", out_csv)
  if (identical(per_cell_csv, out_csv)) per_cell_csv <- paste0(out_csv, "_per_cell.csv")
  write_csv(result$per_cell, per_cell_csv)
  per_cell_csv
}

main <- function() {
  params <- parse_args(list(
    input_csv = "",
    out_csv = "",
    space_name = "",
    batch_col = "batch",
    celltype_col = "cell_type",
    k = 30L,
    coord_prefix = "",
    seed = 1L
  ))
  if (!nzchar(params$input_csv)) stop("--input_csv is required")
  if (!nzchar(params$out_csv)) stop("--out_csv is required")
  if (!nzchar(params$space_name)) params$space_name <- basename(params$input_csv)
  df <- read_csv(params$input_csv)
  result <- compute_score_df(
    df = df,
    space_name = params$space_name,
    batch_col = params$batch_col,
    celltype_col = params$celltype_col,
    k = params$k,
    coord_prefix = params$coord_prefix,
    seed = params$seed
  )
  per_cell_csv <- write_score_outputs(result, params$out_csv)
  cat("Saved summary: ", params$out_csv, "\n", sep = "")
  cat("Saved per-cell scores: ", per_cell_csv, "\n", sep = "")
  print(result$summary)
}

if (!interactive()) {
  cmd <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", cmd, value = TRUE)
  script_name <- if (length(script_arg) > 0) basename(sub("^--file=", "", script_arg[[1]])) else ""
  if (identical(script_name, "09_compute_paper_style_batch_mixing_score.R")) {
    main()
  }
}
