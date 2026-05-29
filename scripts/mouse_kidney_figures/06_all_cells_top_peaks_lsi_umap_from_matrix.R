#!/usr/bin/env Rscript

# All-cell, top-peak matrix-derived UMAP for GSE157079.
#
# This is the next memory-safe step toward PACS paper-style reconstruction:
# all cells are used, informative peaks are selected by detection count, and a
# new TF-IDF/LSI/UMAP embedding is computed from the MatrixMarket count file.
# It does not use precomputed GEO UMAP coordinates, does not run PACS, and does
# not remove batch-effect features.

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
  gse_dir = "/home/woodson/biostatistic/pacs/GSE157079",
  metadata_csv = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv",
  out_dir = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures",
  fig_dir = "/home/woodson/PACS_reproducing/figures/mouse_kidney",
  n_top_peaks = 20000L,
  seed = 1L,
  chunk_lines = 100000L,
  progress_every = 5000000L,
  matrix_file = ""
))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package data.table is required for streaming MatrixMarket chunks.")
}
if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Package irlba is required for LSI.")
}
if (!requireNamespace("uwot", quietly = TRUE)) {
  stop("Package uwot is required for UMAP.")
}

if (!nzchar(params$matrix_file)) {
  params$matrix_file <- file.path(params$gse_dir, "GSE157079_snATAC_cell_by_peak_matrix.txt.gz")
}
peak_file <- file.path(params$gse_dir, "GSE157079_snATAC_peak_list.csv.gz")
if (!file.exists(params$matrix_file)) stop("Missing matrix file: ", params$matrix_file)
if (!file.exists(peak_file)) stop("Missing peak list file: ", peak_file)
if (!file.exists(params$metadata_csv)) stop("Missing metadata CSV: ", params$metadata_csv)

run_dir <- file.path(params$out_dir, "gse157079_all_cells_top_peaks_lsi_umap")
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
    if (ncol(df) < length(expected)) {
      stop("Peak list has unexpected columns: ", paste(names(df), collapse = ", "))
    }
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
  if (length(first_line) != 1 || !grepl("^%%MatrixMarket", first_line)) {
    stop("Matrix file does not start with a MatrixMarket header")
  }
  if (length(dims) < 3) stop("MatrixMarket dimensions line is malformed: ", dims_line)
  list(
    first_line = first_line,
    dims_line = dims_line,
    n_cells = as.integer(dims[[1]]),
    n_peaks = as.integer(dims[[2]]),
    n_nonzero = as.numeric(dims[[3]])
  )
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
  if (length(lines) == 0) {
    return(data.table::data.table(i = integer(), j = integer(), x = numeric()))
  }
  data.table::fread(
    text = paste(lines, collapse = "\n"),
    header = FALSE,
    col.names = c("i", "j", "x"),
    showProgress = FALSE
  )
}

first_pass_counts <- function(path, dims, chunk_lines, progress_every) {
  cell_depth <- numeric(dims$n_cells)
  peak_detection <- numeric(dims$n_peaks)
  processed <- 0
  next_progress <- progress_every

  con <- open_coordinate_stream(path)
  on.exit(close(con), add = TRUE)
  log_msg("First pass start: computing cell depth and peak detection")
  repeat {
    dt <- read_coordinate_chunk(con, chunk_lines)
    if (is.null(dt)) break
    if (nrow(dt) == 0) next
    processed <- processed + nrow(dt)
    cell_depth <- cell_depth + tabulate(dt$i, nbins = dims$n_cells)
    peak_detection <- peak_detection + tabulate(dt$j, nbins = dims$n_peaks)
    if (processed >= next_progress) {
      log_msg("First pass processed ", processed, " coordinate lines")
      next_progress <- next_progress + progress_every
    }
  }
  log_msg("First pass completed: processed ", processed, " coordinate lines")
  list(cell_depth = cell_depth, peak_detection = peak_detection, processed = processed)
}

second_pass_matrix <- function(path, selected_peaks, dims, chunk_lines, progress_every) {
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
  log_msg("Second pass start: building all-cell x top-peak sparse matrix")
  repeat {
    dt <- read_coordinate_chunk(con, chunk_lines)
    if (is.null(dt)) break
    if (nrow(dt) == 0) next
    processed <- processed + nrow(dt)
    mapped_j <- peak_map[dt$j]
    keep <- mapped_j > 0L
    if (any(keep)) {
      chunk_id <- chunk_id + 1L
      i_chunks[[chunk_id]] <- dt$i[keep]
      j_chunks[[chunk_id]] <- mapped_j[keep]
      x_chunks[[chunk_id]] <- dt$x[keep]
      retained <- retained + sum(keep)
    }
    if (processed >= next_progress) {
      log_msg("Second pass processed ", processed, " coordinate lines; retained ", retained, " nonzeros")
      next_progress <- next_progress + progress_every
    }
  }
  log_msg("Second pass completed: processed ", processed, " coordinate lines; retained ", retained, " nonzeros")
  if (retained == 0) stop("No nonzero entries retained for selected top peaks")
  mat <- sparseMatrix(
    i = unlist(i_chunks, use.names = FALSE),
    j = unlist(j_chunks, use.names = FALSE),
    x = unlist(x_chunks, use.names = FALSE),
    dims = c(dims$n_cells, length(selected_peaks))
  )
  list(matrix = as(mat, "dgCMatrix"), processed = processed, retained = retained)
}

make_palette <- function(levels, preferred, fallback) {
  pal <- preferred[names(preferred) %in% levels]
  missing <- setdiff(levels, names(pal))
  if (length(missing) > 0) {
    extra <- rep(fallback, length.out = length(missing))
    names(extra) <- missing
    pal <- c(pal, extra)
  }
  pal[levels]
}

sample_palette_preferred <- c(
  "P0_batch1" = "#e90a0a",
  "P0_batch2" = "#013cbd",
  "P21_batch1" = "#088024",
  "P56_batch1" = "#fbd209",
  "P56_batch2" = "#800dbe"
)
sample_colors <- c("#e90a0a", "#013cbd", "#088024", "#fbd209", "#800dbe", "#11dfdf", "#fc8208")

celltype_palette_preferred <- c(
  "CNT" = "#66C7FF",
  "DCT" = "#D8896A",
  "Endo" = "#C69214",
  "IC" = "#B6A000",
  "immune" = "#7FB000",
  "LOH" = "#67B83F",
  "NP" = "#00A651",
  "NP_LOH" = "#00A98F",
  "PC" = "#00A6B8",
  "Podo" = "#00A7E1",
  "PT" = "#139DDF",
  "PT_out" = "#8F9BEF",
  "PT2" = "#7d4fa0",
  "stroma1" = "#f97ded",
  "stroma2" = "#dc597a"
)
celltype_colors <- c(
  "#66C7FF", "#D8896A", "#C69214", "#B6A000", "#7FB000",
  "#67B83F", "#00A651", "#00A98F", "#00A6B8", "#00A7E1",
  "#139DDF", "#8F9BEF", "#7d4fa0", "#f97ded", "#dc597a",
  "#F4A6C8", "#9ADBC5", "#E6C36A", "#BCA7FF", "#7DD3FC"
)

plot_umap <- function(df, color_col, palette, title, out_base, seed) {
  set.seed(seed)
  plot_df <- df[sample(seq_len(nrow(df))), , drop = FALSE]
  p <- ggplot(plot_df, aes(x = lsi_umap_1, y = lsi_umap_2, color = .data[[color_col]])) +
    geom_point(size = 1.0, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = palette, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_equal(expand = TRUE) +
    labs(x = "LSI UMAP 1", y = "LSI UMAP 2", color = color_col, title = title) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.title = element_text(size = 15),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      legend.key.height = grid::unit(0.42, "cm"),
      legend.key.width = grid::unit(0.42, "cm")
    )
  ggsave(paste0(out_base, ".png"), p, width = 8.5, height = 6.5, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 8.5, height = 6.5)
}

format_table <- function(x, max_rows = 80L) {
  paste(capture.output(print(utils::head(x, max_rows), row.names = FALSE)), collapse = "\n")
}

counts_table <- function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  data.frame(level = names(tab), n = as.integer(tab), stringsAsFactors = FALSE)
}

log_msg("Reading metadata")
metadata <- read_csv(params$metadata_csv)
required_metadata <- c("row_index", "cell_barcode", "sample", "cell_type")
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0) {
  stop("Metadata missing required columns: ", paste(missing_metadata, collapse = ", "))
}
metadata$row_index <- as.integer(metadata$row_index)
metadata$sample <- as.character(metadata$sample)
metadata$cell_type <- as.character(metadata$cell_type)

log_msg("Reading peak list")
peak_list <- standardize_peak_list(read_csv(peak_file))
header <- read_matrix_header(params$matrix_file)
if (nrow(metadata) != header$n_cells) stop("Metadata rows do not match matrix cells")
if (!identical(sort(metadata$row_index), seq_len(header$n_cells))) {
  stop("Metadata row_index does not cover 1:", header$n_cells)
}
if (nrow(peak_list) != header$n_peaks) stop("Peak list rows do not match matrix peaks")
if (!identical(sort(peak_list$peak_index), seq_len(header$n_peaks))) {
  stop("Peak list peak_index does not cover 1:", header$n_peaks)
}
if (params$n_top_peaks > header$n_peaks) {
  stop("Requested n_top_peaks=", params$n_top_peaks, " but matrix only has ", header$n_peaks)
}

pass1 <- first_pass_counts(params$matrix_file, header, params$chunk_lines, params$progress_every)
metadata$cell_depth <- pass1$cell_depth[metadata$row_index]

nonzero_peaks <- which(pass1$peak_detection > 0)
if (length(nonzero_peaks) == 0) stop("No peaks with nonzero detection were found")
n_select <- min(params$n_top_peaks, length(nonzero_peaks))
rank_df <- data.frame(
  peak_index = nonzero_peaks,
  peak_detection = pass1$peak_detection[nonzero_peaks]
)
rank_df <- rank_df[order(-rank_df$peak_detection, rank_df$peak_index), , drop = FALSE]
selected_peaks <- sort(rank_df$peak_index[seq_len(n_select)])
top_peak_info <- peak_list[selected_peaks, , drop = FALSE]
top_peak_info$peak_detection <- pass1$peak_detection[selected_peaks]
top_peak_info <- top_peak_info[order(-top_peak_info$peak_detection, top_peak_info$peak_index), , drop = FALSE]
log_msg("Top peaks selected: requested ", params$n_top_peaks, "; selected ", length(selected_peaks))

top_peaks_out <- file.path(run_dir, "top_peak_indices.csv")
write.csv(top_peak_info, top_peaks_out, row.names = FALSE)

pass2 <- second_pass_matrix(params$matrix_file, selected_peaks, header, params$chunk_lines, params$progress_every)
counts <- pass2$matrix
selected_peak_info <- peak_list[selected_peaks, , drop = FALSE]
rownames(counts) <- paste0("cell_", seq_len(nrow(counts)))
colnames(counts) <- make.unique(as.character(selected_peak_info$name))
log_msg("Sparse matrix built: ", paste(dim(counts), collapse = " x "), "; nnz=", length(counts@x))

pre_filter_dim <- dim(counts)
pre_filter_nnz <- length(counts@x)
nonempty_cells <- Matrix::rowSums(counts) > 0
nonempty_peaks <- Matrix::colSums(counts) > 0
removed_cells <- sum(!nonempty_cells)
removed_peaks <- sum(!nonempty_peaks)
counts <- counts[nonempty_cells, nonempty_peaks, drop = FALSE]
metadata_filtered <- metadata[nonempty_cells, , drop = FALSE]
selected_peak_info <- selected_peak_info[nonempty_peaks, , drop = FALSE]
selected_peak_info$peak_detection <- pass1$peak_detection[selected_peak_info$peak_index]
log_msg("After empty filtering: ", paste(dim(counts), collapse = " x "), "; nnz=", length(counts@x))

if (nrow(counts) < 20 || ncol(counts) < 20) {
  stop("Filtered matrix is too small for LSI/UMAP")
}

counts@x <- rep(1, length(counts@x))
depth_after_filter <- Matrix::rowSums(counts)
peak_detection_after_filter <- Matrix::colSums(counts)
metadata_filtered$top_peak_depth <- as.numeric(depth_after_filter)

counts_rds <- file.path(run_dir, "counts_sparse_top_peaks.rds")
metadata_out <- file.path(run_dir, "metadata_with_depth.csv")
top_peaks_out <- file.path(run_dir, "top_peak_indices.csv")
saveRDS(counts, counts_rds)
write.csv(metadata_filtered, metadata_out, row.names = FALSE)
write.csv(selected_peak_info, top_peaks_out, row.names = FALSE)

log_msg("TF-IDF start")
tf <- Diagonal(x = 1 / as.numeric(depth_after_filter)) %*% counts
idf <- log(1 + nrow(counts) / as.numeric(peak_detection_after_filter))
tfidf <- tf %*% Diagonal(x = idf)
tfidf@x <- log1p(tfidf@x * 1e4)

nv <- min(50L, nrow(tfidf) - 1L, ncol(tfidf) - 1L)
if (nv < 2L) stop("Not enough rows/columns for LSI")
log_msg("LSI start with irlba nv=", nv)
svd_res <- irlba::irlba(tfidf, nv = nv)
lsi <- svd_res$u %*% diag(svd_res$d, nrow = length(svd_res$d))
colnames(lsi) <- paste0("LSI_", seq_len(ncol(lsi)))
rownames(lsi) <- rownames(counts)

depth_cor <- apply(lsi, 2, function(x) suppressWarnings(cor(x, metadata_filtered$top_peak_depth, method = "spearman")))
umap_dims <- if (ncol(lsi) >= 3) 2:min(30L, ncol(lsi)) else seq_len(ncol(lsi))
log_msg("UMAP start with uwot using LSI dims: ", paste(umap_dims, collapse = ","))
set.seed(params$seed)
umap <- uwot::umap(
  lsi[, umap_dims, drop = FALSE],
  n_neighbors = 30,
  min_dist = 0.3,
  metric = "cosine",
  verbose = TRUE,
  ret_model = FALSE
)
colnames(umap) <- c("lsi_umap_1", "lsi_umap_2")
metadata_filtered$lsi_umap_1 <- umap[, 1]
metadata_filtered$lsi_umap_2 <- umap[, 2]
write.csv(metadata_filtered, metadata_out, row.names = FALSE)

sample_levels <- sort(unique(metadata_filtered$sample))
celltype_levels <- sort(unique(metadata_filtered$cell_type))
metadata_filtered$sample <- factor(metadata_filtered$sample, levels = sample_levels)
metadata_filtered$cell_type <- factor(metadata_filtered$cell_type, levels = celltype_levels)
sample_palette <- make_palette(sample_levels, sample_palette_preferred, sample_colors)
celltype_palette <- make_palette(celltype_levels, celltype_palette_preferred, celltype_colors)

sample_base <- file.path(params$fig_dir, "gse157079_all_cells_top_peaks_lsi_umap_by_sample")
celltype_base <- file.path(params$fig_dir, "gse157079_all_cells_top_peaks_lsi_umap_by_celltype")
log_msg("Saving figures")
plot_umap(
  metadata_filtered,
  "sample",
  sample_palette,
  "GSE157079 all-cell top-peak LSI UMAP by sample",
  sample_base,
  params$seed
)
plot_umap(
  metadata_filtered,
  "cell_type",
  celltype_palette,
  "GSE157079 all-cell top-peak LSI UMAP by cell type",
  celltype_base,
  params$seed
)
log_msg("Figures saved")

lsi_out <- file.path(run_dir, "lsi_embedding.csv")
umap_out <- file.path(run_dir, "lsi_umap_embedding.csv")
write.csv(data.frame(cell_id = rownames(lsi), lsi, check.names = FALSE), lsi_out, row.names = FALSE)
write.csv(
  metadata_filtered[, c("row_index", "cell_barcode", "sample", "cell_type", "cell_depth", "top_peak_depth", "lsi_umap_1", "lsi_umap_2")],
  umap_out,
  row.names = FALSE
)

report_path <- file.path(run_dir, "all_cells_top_peaks_lsi_umap_report.md")
report <- c(
  "# GSE157079 All-Cell Top-Peak LSI UMAP Report",
  "",
  "This run computes a new matrix-derived TF-IDF/LSI/UMAP embedding from all cells and top detected peaks.",
  "The precomputed GEO `umap_1`/`umap_2` columns were not used.",
  "PACS and batch-feature removal were not run in this step.",
  "",
  "## Command Arguments",
  "",
  "```text",
  paste(names(params), unlist(params), sep = " = ", collapse = "\n"),
  "```",
  "",
  "## Matrix File",
  "",
  paste0("`", params$matrix_file, "`"),
  "",
  "## MatrixMarket Header",
  "",
  "```text",
  header$first_line,
  header$dims_line,
  "```",
  "",
  "## Metadata",
  "",
  paste0("- Rows: ", nrow(metadata)),
  paste0("- Filtered rows used in UMAP: ", nrow(metadata_filtered)),
  "",
  "### Sample Table",
  "",
  "```text",
  format_table(counts_table(metadata_filtered$sample)),
  "```",
  "",
  "### Cell Type Table",
  "",
  "```text",
  format_table(counts_table(metadata_filtered$cell_type)),
  "```",
  "",
  "## Peak List",
  "",
  paste0("- Rows: ", nrow(peak_list)),
  "",
  "## First Pass Summary",
  "",
  paste0("- Coordinate lines processed: ", pass1$processed),
  paste0("- Peaks with nonzero detection: ", length(nonzero_peaks)),
  "- Cell depth summary:",
  "",
  "```text",
  paste(capture.output(print(summary(pass1$cell_depth))), collapse = "\n"),
  "```",
  "",
  "- Peak detection summary:",
  "",
  "```text",
  paste(capture.output(print(summary(pass1$peak_detection))), collapse = "\n"),
  "```",
  "",
  "## Selected Peak Summary",
  "",
  paste0("- n_top_peaks requested: ", params$n_top_peaks),
  paste0("- n_top_peaks selected: ", length(selected_peaks)),
  "- Selected peak detection summary:",
  "",
  "```text",
  paste(capture.output(print(summary(top_peak_info$peak_detection))), collapse = "\n"),
  "```",
  "",
  "## Sparse Matrix",
  "",
  paste0("- Second pass coordinate lines processed: ", pass2$processed),
  paste0("- Second pass retained nonzeros: ", pass2$retained),
  paste0("- Dimensions before empty filtering: ", paste(pre_filter_dim, collapse = " x ")),
  paste0("- Nonzeros before empty filtering: ", pre_filter_nnz),
  paste0("- Removed empty cells: ", removed_cells),
  paste0("- Removed empty peaks: ", removed_peaks),
  paste0("- Dimensions after empty filtering: ", paste(dim(counts), collapse = " x ")),
  paste0("- Nonzeros after empty filtering: ", length(counts@x)),
  "",
  "## TF-IDF / LSI",
  "",
  paste0("- TF-IDF dimensions: ", paste(dim(tfidf), collapse = " x ")),
  paste0("- LSI dimensions: ", paste(dim(lsi), collapse = " x ")),
  paste0("- LSI dimensions used for UMAP: ", paste(umap_dims, collapse = ", ")),
  "- Spearman correlations between LSI components and top-peak depth:",
  "",
  "```text",
  paste(capture.output(print(depth_cor)), collapse = "\n"),
  "```",
  "",
  "## UMAP Parameters",
  "",
  "- Package: uwot",
  "- n_neighbors: 30",
  "- min_dist: 0.3",
  "- metric: cosine",
  paste0("- UMAP rows: ", nrow(umap)),
  "",
  "## Output Files",
  "",
  paste0("- `", counts_rds, "`"),
  paste0("- `", metadata_out, "`"),
  paste0("- `", top_peaks_out, "`"),
  paste0("- `", lsi_out, "`"),
  paste0("- `", umap_out, "`"),
  paste0("- `", sample_base, ".png`"),
  paste0("- `", sample_base, ".pdf`"),
  paste0("- `", celltype_base, ".png`"),
  paste0("- `", celltype_base, ".pdf`"),
  "",
  "## Conclusion",
  "",
  "All-cell top-peak matrix-derived UMAP completed successfully if this report was written.",
  "This is the before-PACS-filtering reference candidate, not the final PACS-filtered UMAP."
)
writeLines(report, report_path)

log_msg("Saved report: ", report_path)
cat("First pass completed: TRUE\n")
cat("Second pass completed: TRUE\n")
cat("Selected top peak count: ", length(selected_peaks), "\n", sep = "")
cat("Sparse matrix dimensions: ", paste(dim(counts), collapse = " x "), "\n", sep = "")
cat("Sparse matrix nnz: ", length(counts@x), "\n", sep = "")
cat("LSI completed: TRUE\n")
cat("UMAP completed: TRUE\n")
cat("Sample figure: ", sample_base, ".png / .pdf\n", sep = "")
cat("Cell type figure: ", celltype_base, ".png / .pdf\n", sep = "")
