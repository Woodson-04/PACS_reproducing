#!/usr/bin/env Rscript

# P56-specific PACS batch-effect peak filtering demonstration.
#
# This script builds a P56-only cell x peak sparse matrix from the GSE157079
# MatrixMarket file, computes a before-filtering TF-IDF/LSI/UMAP, tests
# P56_batch1 vs P56_batch2 batch-associated peaks with PACS while adjusting for
# cell_type, removes significant batch peaks, and recomputes after-filtering
# TF-IDF/LSI/UMAP. This is a first PACS paper-style demonstration, not the
# final full-dataset analysis.

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(PACS)
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
  gse_dir = "/home/woodson/biostatistic/pacs/GSE157079",
  metadata_csv = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv",
  out_dir = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures",
  fig_dir = "/home/woodson/PACS_reproducing/figures/mouse_kidney",
  matrix_file = "",
  n_top_peaks = 20000L,
  max_pacs_peaks = 20000L,
  fdr_cutoff = 0.05,
  seed = 1L,
  chunk_lines = 100000L,
  progress_every = 5000000L,
  min_cells_per_celltype = 30L,
  min_cells_per_batch_per_celltype = 5L,
  pacs_peaks_per_round = 1000L,
  run_name = ""
))

for (pkg in c("data.table", "irlba", "uwot")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package ", pkg, " is required for this script.")
  }
}

if (!nzchar(params$matrix_file)) {
  params$matrix_file <- file.path(params$gse_dir, "GSE157079_snATAC_cell_by_peak_matrix.txt.gz")
}
peak_file <- file.path(params$gse_dir, "GSE157079_snATAC_peak_list.csv.gz")
if (!file.exists(params$matrix_file)) stop("Missing matrix file: ", params$matrix_file)
if (!file.exists(peak_file)) stop("Missing peak list file: ", peak_file)
if (!file.exists(params$metadata_csv)) stop("Missing metadata CSV: ", params$metadata_csv)

use_run_name <- nzchar(params$run_name) && params$run_name != "default"
result_dir_name <- if (use_run_name) {
  paste0("gse157079_p56_pacs_batch_filter_umap_", params$run_name)
} else {
  "gse157079_p56_pacs_batch_filter_umap"
}
figure_prefix_name <- if (use_run_name) {
  paste0("gse157079_p56_pacs_batch_filter_", params$run_name)
} else {
  "gse157079_p56_pacs_batch_filter"
}
run_dir <- file.path(params$out_dir, result_dir_name)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""), "\n", sep = "")
}

read_csv <- function(path) {
  as.data.frame(data.table::fread(path, header = TRUE, showProgress = FALSE, check.names = FALSE))
}

standardize_peak_list <- function(df) {
  required <- c("seqnames", "start", "end", "name")
  if (!all(required %in% names(df))) {
    expected <- c("peak_index", "seqnames", "start", "end", "width", "strand", "name")
    if (ncol(df) < length(expected)) stop("Peak list has unexpected columns: ", paste(names(df), collapse = ", "))
    df <- df[, seq_along(expected), drop = FALSE]
    names(df) <- expected
  } else if (!"peak_index" %in% names(df)) {
    first_name <- names(df)[[1]]
    if (is.na(first_name) || first_name == "" || grepl("^V1$|^\\.\\.\\.1$|^X$", first_name)) {
      names(df)[[1]] <- "peak_index"
    } else {
      df$peak_index <- seq_len(nrow(df))
    }
  }
  df$peak_index <- as.integer(df$peak_index)
  df
}

read_matrix_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  first_line <- readLines(con, n = 1, warn = FALSE)
  dims_line <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "%")) next
    dims_line <- line
    break
  }
  dims <- scan(text = dims_line, quiet = TRUE)
  if (length(first_line) != 1 || !grepl("^%%MatrixMarket", first_line)) stop("Matrix file does not start with MatrixMarket header")
  if (length(dims) < 3) stop("MatrixMarket dimensions line is malformed: ", dims_line)
  list(first_line = first_line, dims_line = dims_line, n_cells = as.integer(dims[1]), n_peaks = as.integer(dims[2]), n_nonzero = as.numeric(dims[3]))
}

open_coordinate_stream <- function(path) {
  con <- gzfile(path, open = "rt")
  invisible(readLines(con, n = 1, warn = FALSE))
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) {
      close(con)
      stop("Could not find MatrixMarket dimensions line")
    }
    if (!startsWith(line, "%")) break
  }
  con
}

read_coordinate_chunk <- function(con, n) {
  lines <- readLines(con, n = n, warn = FALSE)
  if (length(lines) == 0) return(NULL)
  lines <- lines[nzchar(lines)]
  lines <- lines[!startsWith(lines, "%")]
  if (length(lines) == 0) return(data.table::data.table(i = integer(), j = integer(), x = numeric()))
  data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE, col.names = c("i", "j", "x"), showProgress = FALSE)
}

filter_p56_metadata <- function(metadata) {
  p56 <- metadata[metadata$sample %in% c("P56_batch1", "P56_batch2"), , drop = FALSE]
  p56$batch <- factor(p56$sample, levels = c("P56_batch1", "P56_batch2"))
  total_tab <- table(p56$cell_type)
  cross_tab <- table(p56$cell_type, p56$batch)
  keep_types <- names(total_tab)[total_tab >= params$min_cells_per_celltype]
  keep_types <- keep_types[keep_types %in% rownames(cross_tab)[apply(cross_tab[keep_types, , drop = FALSE], 1, function(x) all(x >= params$min_cells_per_batch_per_celltype))]]
  removed <- setdiff(sort(unique(p56$cell_type)), keep_types)
  filtered <- p56[p56$cell_type %in% keep_types, , drop = FALSE]
  filtered <- filtered[order(filtered$row_index), , drop = FALSE]
  if (nrow(filtered) < 100) stop("Cell type filtering retained too few P56 cells: ", nrow(filtered))
  list(metadata = filtered, removed_cell_types = removed, pre_filter = p56, cross_tab = cross_tab)
}

first_pass_p56 <- function(path, selected_cells, dims, chunk_lines, progress_every) {
  cell_map <- integer(dims$n_cells)
  cell_map[selected_cells] <- seq_along(selected_cells)
  cell_depth <- numeric(length(selected_cells))
  peak_detection <- numeric(dims$n_peaks)
  processed <- 0
  retained <- 0
  next_progress <- progress_every
  con <- open_coordinate_stream(path)
  on.exit(close(con), add = TRUE)
  log_msg("First matrix pass start for selected P56 cells")
  repeat {
    dt <- read_coordinate_chunk(con, chunk_lines)
    if (is.null(dt)) break
    if (nrow(dt) == 0) next
    processed <- processed + nrow(dt)
    mapped_i <- cell_map[dt$i]
    keep <- mapped_i > 0L
    if (any(keep)) {
      cell_depth <- cell_depth + tabulate(mapped_i[keep], nbins = length(selected_cells))
      peak_detection <- peak_detection + tabulate(dt$j[keep], nbins = dims$n_peaks)
      retained <- retained + sum(keep)
    }
    if (processed >= next_progress) {
      log_msg("First pass processed ", processed, " lines; P56 retained entries ", retained)
      next_progress <- next_progress + progress_every
    }
  }
  log_msg("First pass completed: processed ", processed, "; selected-cell entries ", retained)
  list(cell_depth = cell_depth, peak_detection = peak_detection, processed = processed, retained = retained)
}

second_pass_p56_matrix <- function(path, selected_cells, selected_peaks, dims, chunk_lines, progress_every) {
  cell_map <- integer(dims$n_cells)
  cell_map[selected_cells] <- seq_along(selected_cells)
  peak_map <- integer(dims$n_peaks)
  peak_map[selected_peaks] <- seq_along(selected_peaks)
  i_chunks <- list()
  j_chunks <- list()
  x_chunks <- list()
  chunk_id <- 0L
  processed <- 0
  retained <- 0
  next_progress <- progress_every
  con <- open_coordinate_stream(path)
  on.exit(close(con), add = TRUE)
  log_msg("Second matrix pass start: building P56 cell x top-peak sparse matrix")
  repeat {
    dt <- read_coordinate_chunk(con, chunk_lines)
    if (is.null(dt)) break
    if (nrow(dt) == 0) next
    processed <- processed + nrow(dt)
    mapped_i <- cell_map[dt$i]
    mapped_j <- peak_map[dt$j]
    keep <- mapped_i > 0L & mapped_j > 0L
    if (any(keep)) {
      chunk_id <- chunk_id + 1L
      i_chunks[[chunk_id]] <- mapped_i[keep]
      j_chunks[[chunk_id]] <- mapped_j[keep]
      x_chunks[[chunk_id]] <- dt$x[keep]
      retained <- retained + sum(keep)
    }
    if (processed >= next_progress) {
      log_msg("Second pass processed ", processed, " lines; retained ", retained, " nonzeros")
      next_progress <- next_progress + progress_every
    }
  }
  log_msg("Second pass completed: processed ", processed, "; retained ", retained)
  if (retained == 0) stop("No nonzeros retained in P56 top-peak matrix")
  mat <- sparseMatrix(i = unlist(i_chunks, use.names = FALSE), j = unlist(j_chunks, use.names = FALSE), x = unlist(x_chunks, use.names = FALSE), dims = c(length(selected_cells), length(selected_peaks)))
  list(matrix = as(mat, "dgCMatrix"), processed = processed, retained = retained)
}

compute_lsi_umap <- function(counts, metadata, prefix, seed) {
  counts <- as(counts, "dgCMatrix")
  counts@x <- rep(1, length(counts@x))
  depth <- Matrix::rowSums(counts)
  peak_detection <- Matrix::colSums(counts)
  if (any(depth <= 0) || any(peak_detection <= 0)) stop(prefix, ": empty rows/columns passed into TF-IDF")
  log_msg(prefix, " TF-IDF start")
  tf <- Diagonal(x = 1 / as.numeric(depth)) %*% counts
  idf <- log(1 + nrow(counts) / as.numeric(peak_detection))
  tfidf <- tf %*% Diagonal(x = idf)
  tfidf@x <- log1p(tfidf@x * 1e4)
  nv <- min(50L, nrow(tfidf) - 1L, ncol(tfidf) - 1L)
  if (nv < 2L) stop(prefix, ": not enough rows/columns for LSI")
  log_msg(prefix, " LSI start with nv=", nv)
  svd_res <- irlba::irlba(tfidf, nv = nv)
  lsi <- svd_res$u %*% diag(svd_res$d, nrow = length(svd_res$d))
  colnames(lsi) <- paste0("LSI_", seq_len(ncol(lsi)))
  depth_cor <- apply(lsi, 2, function(x) suppressWarnings(cor(x, as.numeric(depth), method = "spearman")))
  umap_dims <- if (ncol(lsi) >= 3) 2:min(30L, ncol(lsi)) else seq_len(ncol(lsi))
  log_msg(prefix, " UMAP start using LSI dims ", paste(umap_dims, collapse = ","))
  set.seed(seed)
  umap <- uwot::umap(lsi[, umap_dims, drop = FALSE], n_neighbors = 30, min_dist = 0.3, metric = "cosine", verbose = TRUE, ret_model = FALSE)
  colnames(umap) <- c(paste0(prefix, "_umap_1"), paste0(prefix, "_umap_2"))
  out <- metadata
  out[[paste0(prefix, "_depth")]] <- as.numeric(depth)
  out[[paste0(prefix, "_umap_1")]] <- umap[, 1]
  out[[paste0(prefix, "_umap_2")]] <- umap[, 2]
  lsi_out <- metadata
  lsi_out[[paste0(prefix, "_depth")]] <- as.numeric(depth)
  lsi_out <- cbind(lsi_out, as.data.frame(lsi, check.names = FALSE))
  list(
    metadata = out,
    lsi_embedding = lsi_out,
    tfidf_dim = dim(tfidf),
    lsi_dim = dim(lsi),
    umap_dims = umap_dims,
    depth_cor = depth_cor
  )
}

get_pacs_fun <- function(name) {
  if (exists(name, where = asNamespace("PACS"), inherits = FALSE)) get(name, envir = asNamespace("PACS")) else NULL
}

pacs_test_sparse_local_fixed <- function(covariate_meta.data, formula_full, formula_null, pic_matrix, n_peaks_per_round = 1000L, T_proportion_cutoff = 0.2, cap_rates, label = "") {
  pacs_test_cumu <- get_pacs_fun("pacs_test_cumu")
  pacs_test_logit <- get_pacs_fun("pacs_test_logit")
  if (is.null(pacs_test_cumu) || is.null(pacs_test_logit)) stop("PACS pacs_test_cumu/logit functions are unavailable")
  if (!inherits(pic_matrix, "Matrix")) stop(label, ": PACS pic_matrix must be sparse Matrix peaks x cells")
  n_cell <- ncol(pic_matrix)
  n_peaks <- nrow(pic_matrix)
  p_names <- rownames(pic_matrix)
  if (nrow(covariate_meta.data) != n_cell) stop(label, ": metadata rows do not match PACS matrix cells")
  if (length(cap_rates) != n_cell) stop(label, ": cap_rates length does not match PACS matrix cells")
  if (is.null(p_names) || length(p_names) != n_peaks) {
    p_names <- paste0("peak_", seq_len(n_peaks))
    rownames(pic_matrix) <- p_names
  }
  if (is.null(colnames(pic_matrix))) colnames(pic_matrix) <- paste0("cell_", seq_len(n_cell))
  n_peaks_per_round <- min(as.integer(n_peaks_per_round), n_peaks)
  log_msg("PACS local wrapper ", label, ": dim=", paste(dim(pic_matrix), collapse = " x "), "; n_peaks_per_round=", n_peaks_per_round)
  pic_matrix_2 <- pic_matrix
  pic_matrix_2@x[pic_matrix_2@x == 1] <- 0
  pic_matrix_2 <- Matrix::drop0(pic_matrix_2)
  pic_matrix_2@x <- rep(1, length(pic_matrix_2@x))
  pic_matrixbin <- pic_matrix
  pic_matrixbin@x <- rep(1, length(pic_matrixbin@x))
  rs <- Matrix::rowSums(pic_matrixbin)
  rs2 <- Matrix::rowSums(pic_matrix_2)
  p_2 <- rs2 / rs
  p_2[is.na(p_2)] <- 0
  f_sel <- names(p_2)[p_2 >= T_proportion_cutoff]
  f_b_sel <- names(p_2)[p_2 < T_proportion_cutoff]
  p_cumu <- list()
  p_logit <- list()
  if (length(f_sel) > 0) {
    mat <- pic_matrix[f_sel, , drop = FALSE]
    n_iter <- ceiling(nrow(mat) / n_peaks_per_round)
    for (jj in seq_len(n_iter)) {
      st <- (jj - 1L) * n_peaks_per_round + 1L
      en <- min(nrow(mat), jj * n_peaks_per_round)
      dense <- as.matrix(mat[st:en, , drop = FALSE])
      log_msg(label, " PACS cumulative block ", jj, "/", n_iter, " dim=", paste(dim(dense), collapse = " x "))
      p_cumu[[jj]] <- pacs_test_cumu(covariate_meta.data = covariate_meta.data, max_T = 2, formula_full = formula_full, formula_null = formula_null, pic_matrix = dense, cap_rates = cap_rates, n_cores = 1)
    }
  }
  if (length(f_b_sel) > 0) {
    mat <- pic_matrixbin[f_b_sel, , drop = FALSE]
    n_iter <- ceiling(nrow(mat) / n_peaks_per_round)
    for (jj in seq_len(n_iter)) {
      st <- (jj - 1L) * n_peaks_per_round + 1L
      en <- min(nrow(mat), jj * n_peaks_per_round)
      dense <- as.matrix(mat[st:en, , drop = FALSE])
      log_msg(label, " PACS logit block ", jj, "/", n_iter, " dim=", paste(dim(dense), collapse = " x "))
      p_logit[[jj]] <- pacs_test_logit(covariate_meta.data = covariate_meta.data, formula_full = formula_full, formula_null = formula_null, pic_matrix = dense, cap_rates = cap_rates, n_cores = 1)
    }
  }
  p_val_cumu <- if (length(p_cumu)) unlist(lapply(p_cumu, function(x) x$pacs_p_val), use.names = TRUE) else numeric(0)
  p_val_logit <- if (length(p_logit)) unlist(lapply(p_logit, function(x) x$pacs_p_val), use.names = TRUE) else numeric(0)
  p_val <- c(p_val_cumu, p_val_logit)[p_names]
  if (length(p_val) != n_peaks || !identical(names(p_val), p_names) || anyNA(p_val)) stop(label, ": PACS p-value merge failed")
  list(pacs_converged = NULL, pacs_p_val = p_val)
}

sample_palette <- c("P56_batch1" = "#FFD400", "P56_batch2" = "#A000FF")
celltype_palette_preferred <- c("CNT" = "#66C7FF", "DCT" = "#D8896A", "Endo" = "#C69214", "IC" = "#B6A000", "immune" = "#7FB000", "LOH" = "#67B83F", "NP" = "#00A651", "NP_LOH" = "#00A98F", "PC" = "#00A6B8", "Podo" = "#00A7E1", "PT" = "#139DDF", "PT_out" = "#8F9BEF", "PT2" = "#B879E8", "stroma1" = "#DD70D6", "stroma2" = "#F06AA8")
celltype_fallback <- c("#66C7FF", "#D8896A", "#C69214", "#B6A000", "#7FB000", "#67B83F", "#00A651", "#00A98F", "#00A6B8", "#00A7E1", "#139DDF", "#8F9BEF", "#B879E8", "#DD70D6", "#F06AA8")
make_palette <- function(levels, preferred, fallback) {
  pal <- preferred[names(preferred) %in% levels]
  missing <- setdiff(levels, names(pal))
  if (length(missing)) {
    extra <- rep(fallback, length.out = length(missing))
    names(extra) <- missing
    pal <- c(pal, extra)
  }
  pal[levels]
}

plot_embedding <- function(df, xcol, ycol, color_col, palette, title, out_base, seed) {
  set.seed(seed)
  plot_df <- df[sample(seq_len(nrow(df))), , drop = FALSE]
  p <- ggplot(plot_df, aes(x = .data[[xcol]], y = .data[[ycol]], color = .data[[color_col]])) +
    geom_point(size = 0.2, alpha = 0.85, stroke = 0) +
    scale_color_manual(values = palette, drop = FALSE) +
    coord_equal(expand = TRUE) +
    labs(x = "UMAP 1", y = "UMAP 2", color = color_col, title = title) +
    theme_classic(base_size = 15) +
    theme(plot.title = element_text(face = "bold", size = 17), axis.text = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size = 10))
  ggsave(paste0(out_base, ".png"), p, width = 8.5, height = 6.5, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 8.5, height = 6.5)
  p
}

format_table <- function(x, max_rows = 80L) paste(capture.output(print(utils::head(x, max_rows), row.names = FALSE)), collapse = "\n")
counts_table <- function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  data.frame(level = names(tab), n = as.integer(tab), stringsAsFactors = FALSE)
}

log_msg("Reading metadata")
metadata <- read_csv(params$metadata_csv)
required_metadata <- c("row_index", "cell_barcode", "sample", "cell_type")
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata)) stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
metadata$row_index <- as.integer(metadata$row_index)
metadata$sample <- as.character(metadata$sample)
metadata$cell_type <- as.character(metadata$cell_type)

log_msg("Subsetting P56 cells and filtering cell types")
p56_info <- filter_p56_metadata(metadata)
p56_meta <- p56_info$metadata
p56_meta$batch <- factor(p56_meta$batch, levels = c("P56_batch1", "P56_batch2"))
p56_meta$cell_type <- factor(p56_meta$cell_type)
log_msg("P56 cells retained: ", nrow(p56_meta), "; batches: ", paste(names(table(p56_meta$batch)), as.integer(table(p56_meta$batch)), sep = "=", collapse = ", "))

log_msg("Reading peak list and matrix header")
peak_list <- standardize_peak_list(read_csv(peak_file))
header <- read_matrix_header(params$matrix_file)
if (nrow(peak_list) != header$n_peaks) stop("Peak list rows do not match matrix peaks")

selected_cells <- p56_meta$row_index
pass1 <- first_pass_p56(params$matrix_file, selected_cells, header, params$chunk_lines, params$progress_every)
p56_meta$full_matrix_depth <- pass1$cell_depth
nonzero_peaks <- which(pass1$peak_detection > 0)
n_select <- min(params$n_top_peaks, length(nonzero_peaks))
rank_df <- data.frame(peak_index = nonzero_peaks, peak_detection = pass1$peak_detection[nonzero_peaks])
rank_df <- rank_df[order(-rank_df$peak_detection, rank_df$peak_index), , drop = FALSE]
selected_peaks <- sort(rank_df$peak_index[seq_len(n_select)])
selected_peak_info <- peak_list[selected_peaks, , drop = FALSE]
selected_peak_info$peak_detection <- pass1$peak_detection[selected_peaks]
log_msg("Top peaks selected: ", length(selected_peaks))

pass2 <- second_pass_p56_matrix(params$matrix_file, selected_cells, selected_peaks, header, params$chunk_lines, params$progress_every)
counts <- pass2$matrix
rownames(counts) <- paste0("cell_", p56_meta$row_index)
colnames(counts) <- make.unique(as.character(selected_peak_info$name))
counts@x <- rep(1, length(counts@x))
pre_filter_dim <- dim(counts)
pre_filter_nnz <- length(counts@x)
nonempty_cells <- Matrix::rowSums(counts) > 0
nonempty_peaks <- Matrix::colSums(counts) > 0
removed_empty_cells <- sum(!nonempty_cells)
removed_empty_peaks <- sum(!nonempty_peaks)
counts <- counts[nonempty_cells, nonempty_peaks, drop = FALSE]
p56_meta <- p56_meta[nonempty_cells, , drop = FALSE]
selected_peak_info <- selected_peak_info[nonempty_peaks, , drop = FALSE]
selected_peak_info$top_peak_col <- seq_len(nrow(selected_peak_info))
log_msg("Sparse matrix built after empty filtering: ", paste(dim(counts), collapse = " x "), "; nnz=", length(counts@x))

counts_rds <- file.path(run_dir, "p56_counts_top_peaks_sparse.rds")
metadata_out <- file.path(run_dir, "p56_metadata.csv")
top_peaks_out <- file.path(run_dir, "p56_top_peak_indices.csv")
saveRDS(counts, counts_rds)
write.csv(p56_meta, metadata_out, row.names = FALSE)
write.csv(selected_peak_info, top_peaks_out, row.names = FALSE)

before <- compute_lsi_umap(counts, p56_meta, "before", params$seed)
write.csv(before$metadata, file.path(run_dir, "p56_before_lsi_umap_embedding.csv"), row.names = FALSE)
write.csv(before$lsi_embedding, file.path(run_dir, "p56_before_lsi_embedding.csv"), row.names = FALSE)
log_msg("Saved before LSI embedding: ", file.path(run_dir, "p56_before_lsi_embedding.csv"))

log_msg("Starting PACS batch-effect testing")
pacs_n <- min(params$max_pacs_peaks, ncol(counts))
pacs_rank <- selected_peak_info[order(-selected_peak_info$peak_detection, selected_peak_info$peak_index), , drop = FALSE]
pacs_peak_info <- pacs_rank[seq_len(pacs_n), , drop = FALSE]
pacs_cols <- match(pacs_peak_info$peak_index, selected_peak_info$peak_index)
if (anyNA(pacs_cols)) stop("Could not map PACS tested peaks to count matrix columns")
pacs_counts <- counts[, pacs_cols, drop = FALSE]
pacs_pic <- as(t(pacs_counts), "dgCMatrix")
rownames(pacs_pic) <- paste0("peak_", pacs_peak_info$peak_index)
colnames(pacs_pic) <- paste0("cell_", seq_len(ncol(pacs_pic)))
pacs_meta <- data.frame(cell_type = factor(p56_meta$cell_type), batch = factor(p56_meta$batch, levels = c("P56_batch1", "P56_batch2")))
cap_rates <- as.numeric(Matrix::rowSums(counts))
cap_rates <- pmax(1e-6, pmin(0.99, cap_rates / max(cap_rates)))
if (nrow(pacs_meta) != ncol(pacs_pic) || length(cap_rates) != ncol(pacs_pic)) stop("PACS input dimensions do not align")
pacs_sig <- capture.output({
  cat("PACS package path:\n")
  print(system.file(package = "PACS"))
  cat("\nPACS version:\n")
  print(packageVersion("PACS"))
  cat("\npacs_test_cumu args:\n")
  print(args(PACS::pacs_test_cumu))
  cat("\npacs_test_logit args:\n")
  print(args(PACS::pacs_test_logit))
})
writeLines(pacs_sig, file.path(run_dir, "pacs_function_signatures.txt"))
pacs_result <- pacs_test_sparse_local_fixed(covariate_meta.data = pacs_meta, formula_full = ~ cell_type + batch, formula_null = ~ cell_type, pic_matrix = pacs_pic, n_peaks_per_round = params$pacs_peaks_per_round, T_proportion_cutoff = 0.2, cap_rates = cap_rates, label = "P56 batch")
pvals <- pacs_result$pacs_p_val
pacs_df <- pacs_peak_info
pacs_df$p_value <- as.numeric(pvals[paste0("peak_", pacs_df$peak_index)])
pacs_df$fdr <- p.adjust(pacs_df$p_value, method = "BH")
pacs_df$is_batch_peak <- pacs_df$fdr <= params$fdr_cutoff
pacs_out <- file.path(run_dir, "p56_pacs_batch_peak_results.csv")
write.csv(pacs_df, pacs_out, row.names = FALSE)
batch_peak_n <- sum(pacs_df$is_batch_peak, na.rm = TRUE)
log_msg("PACS completed: tested ", nrow(pacs_df), " peaks; FDR significant batch peaks=", batch_peak_n)

tested_peak_indices <- pacs_df$peak_index
batch_peak_indices <- pacs_df$peak_index[pacs_df$is_batch_peak]
retained_tested <- setdiff(tested_peak_indices, batch_peak_indices)
if (length(retained_tested) == 0) {
  log_msg("No tested peaks retained after filtering; using original tested peaks for after UMAP with warning")
  retained_tested <- tested_peak_indices
}
retained_cols <- match(retained_tested, selected_peak_info$peak_index)
after_counts <- counts[, retained_cols, drop = FALSE]
retained_info <- selected_peak_info[retained_cols, , drop = FALSE]
retained_info$is_retained_after_pacs <- TRUE
write.csv(retained_info, file.path(run_dir, "p56_retained_peak_indices.csv"), row.names = FALSE)
summary_df <- data.frame(fdr_cutoff = params$fdr_cutoff, tested_peaks = nrow(pacs_df), significant_batch_peaks = batch_peak_n, retained_tested_peaks = length(retained_tested), stringsAsFactors = FALSE)
write.csv(summary_df, file.path(run_dir, "p56_batch_peak_summary.csv"), row.names = FALSE)

if (batch_peak_n == 0) {
  after <- before
  after$metadata <- before$metadata
  names(after$metadata)[names(after$metadata) == "before_umap_1"] <- "after_umap_1"
  names(after$metadata)[names(after$metadata) == "before_umap_2"] <- "after_umap_2"
  after$lsi_embedding <- before$lsi_embedding
  names(after$lsi_embedding)[names(after$lsi_embedding) == "before_depth"] <- "after_depth"
  log_msg("No significant batch peaks; after embedding copied from before metadata for plotting clarity")
} else {
  after <- compute_lsi_umap(after_counts, p56_meta, "after", params$seed)
}
write.csv(after$metadata, file.path(run_dir, "p56_after_lsi_umap_embedding.csv"), row.names = FALSE)
write.csv(after$lsi_embedding, file.path(run_dir, "p56_after_lsi_embedding.csv"), row.names = FALSE)
log_msg("Saved after LSI embedding: ", file.path(run_dir, "p56_after_lsi_embedding.csv"))

celltype_levels <- sort(unique(as.character(p56_meta$cell_type)))
celltype_palette <- make_palette(celltype_levels, celltype_palette_preferred, celltype_fallback)
prefix <- file.path(params$fig_dir, figure_prefix_name)
p1 <- plot_embedding(before$metadata, "before_umap_1", "before_umap_2", "batch", sample_palette, "P56 before PACS filtering by batch", paste0(prefix, "_before_by_batch"), params$seed)
p2 <- plot_embedding(before$metadata, "before_umap_1", "before_umap_2", "cell_type", celltype_palette, "P56 before PACS filtering by cell type", paste0(prefix, "_before_by_celltype"), params$seed)
p3 <- plot_embedding(after$metadata, "after_umap_1", "after_umap_2", "batch", sample_palette, "P56 after PACS batch-peak filtering by batch", paste0(prefix, "_after_by_batch"), params$seed)
p4 <- plot_embedding(after$metadata, "after_umap_1", "after_umap_2", "cell_type", celltype_palette, "P56 after PACS batch-peak filtering by cell type", paste0(prefix, "_after_by_celltype"), params$seed)
combined_status <- "combined figure skipped; patchwork/cowplot unavailable"
if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- (p1 + p2) / (p3 + p4)
  ggsave(paste0(prefix, "_four_panel.png"), combined, width = 14, height = 11, dpi = 300)
  ggsave(paste0(prefix, "_four_panel.pdf"), combined, width = 14, height = 11)
  combined_status <- "combined figure generated with patchwork"
} else if (requireNamespace("cowplot", quietly = TRUE)) {
  combined <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)
  ggsave(paste0(prefix, "_four_panel.png"), combined, width = 14, height = 11, dpi = 300)
  ggsave(paste0(prefix, "_four_panel.pdf"), combined, width = 14, height = 11)
  combined_status <- "combined figure generated with cowplot"
}
log_msg("Figures saved")

report <- c(
  "# P56 PACS Batch-Effect Peak Filtering UMAP Report",
  "",
  "## Goal",
  "Generate a first P56-specific PACS paper-style before/after UMAP demonstration.",
  "",
  "## Input Files",
  paste0("- matrix_file: `", params$matrix_file, "`"),
  paste0("- metadata_csv: `", params$metadata_csv, "`"),
  paste0("- peak_file: `", peak_file, "`"),
  "",
  "## P56 Cell Selection",
  paste0("- P56 cells before cell type filtering: ", nrow(p56_info$pre_filter)),
  paste0("- P56 cells after cell type filtering: ", nrow(p56_meta)),
  "",
  "### Batch Table",
  "```text",
  format_table(counts_table(p56_meta$batch)),
  "```",
  "",
  "### Cell Type Table",
  "```text",
  format_table(counts_table(p56_meta$cell_type)),
  "```",
  "",
  "## Cell Type Filtering",
  paste0("- min_cells_per_celltype: ", params$min_cells_per_celltype),
  paste0("- min_cells_per_batch_per_celltype: ", params$min_cells_per_batch_per_celltype),
  paste0("- removed cell types: ", if (length(p56_info$removed_cell_types)) paste(p56_info$removed_cell_types, collapse = ", ") else "none"),
  "",
  "## Matrix Dimensions",
  "```text",
  header$first_line,
  header$dims_line,
  "```",
  paste0("- P56 sparse matrix before empty filtering: ", paste(pre_filter_dim, collapse = " x "), "; nnz=", pre_filter_nnz),
  paste0("- removed empty cells: ", removed_empty_cells),
  paste0("- removed empty peaks: ", removed_empty_peaks),
  paste0("- P56 sparse matrix after empty filtering: ", paste(dim(counts), collapse = " x "), "; nnz=", length(counts@x)),
  "",
  "## Top Peak Selection",
  paste0("- requested n_top_peaks: ", params$n_top_peaks),
  paste0("- selected top peaks: ", nrow(selected_peak_info)),
  paste0("- first pass processed coordinate lines: ", pass1$processed),
  paste0("- selected P56 entries in first pass: ", pass1$retained),
  "",
  "## Before-Filtering UMAP Settings",
  paste0("- TF-IDF dim: ", paste(before$tfidf_dim, collapse = " x ")),
  paste0("- LSI dim: ", paste(before$lsi_dim, collapse = " x ")),
  paste0("- UMAP LSI dims: ", paste(before$umap_dims, collapse = ", ")),
  "",
  "## PACS Function Signatures Inspected",
  "See `pacs_function_signatures.txt`.",
  "",
  "## PACS Model Used",
  "- PACS input orientation: peaks x cells.",
  "- Original matrix orientation: cells x peaks.",
  "- full model: `~ cell_type + batch`",
  "- null model: `~ cell_type`",
  "- `cap_rates` are relative cell-depth rates from selected top peaks, scaled to max 0.99. This is a practical first-demonstration approximation because GSE157079 does not provide the Notebook 1 author `q_vec`.",
  paste0("- PACS peaks per round: ", params$pacs_peaks_per_round),
  "",
  "## PACS Results",
  paste0("- tested peaks: ", nrow(pacs_df)),
  paste0("- FDR cutoff: ", params$fdr_cutoff),
  paste0("- significant batch peaks: ", batch_peak_n),
  paste0("- retained tested peaks: ", length(retained_tested)),
  "",
  "## After-Filtering UMAP Settings",
  paste0("- TF-IDF dim: ", paste(after$tfidf_dim, collapse = " x ")),
  paste0("- LSI dim: ", paste(after$lsi_dim, collapse = " x ")),
  paste0("- UMAP LSI dims: ", paste(after$umap_dims, collapse = ", ")),
  "",
  "## Output Figures",
  paste0("- `", prefix, "_before_by_batch.png/pdf`"),
  paste0("- `", prefix, "_before_by_celltype.png/pdf`"),
  paste0("- `", prefix, "_after_by_batch.png/pdf`"),
  paste0("- `", prefix, "_after_by_celltype.png/pdf`"),
  paste0("- combined: ", combined_status),
  "",
  "## Interpretation",
  "Compare the before/after batch-colored plots to assess whether P56_batch1 and P56_batch2 separation weakens after removing PACS-significant batch peaks. Compare cell-type-colored plots to check whether biological cell type structure remains.",
  "",
  "## Limitations",
  "- P56-only analysis.",
  "- Top detected peaks only.",
  "- First author-style reconstruction, not the final full-dataset PACS paper figure.",
  "- Uses relative depth-derived cap_rates because GSE157079 lacks the author-provided q_vec used in Notebook 1."
)
writeLines(report, file.path(run_dir, "p56_pacs_batch_filter_umap_report.md"))
log_msg("Saved report: ", file.path(run_dir, "p56_pacs_batch_filter_umap_report.md"))
cat("P56 cell count: ", nrow(p56_meta), "\n", sep = "")
cat("Selected top peak count: ", nrow(selected_peak_info), "\n", sep = "")
cat("PACS tested peak count: ", nrow(pacs_df), "\n", sep = "")
cat("FDR significant batch peaks: ", batch_peak_n, "\n", sep = "")
cat("Retained peak count: ", length(retained_tested), "\n", sep = "")
cat("Before/after UMAP figures generated: TRUE\n")
