#!/usr/bin/env Rscript

# Pilot matrix-derived UMAP for GSE157079.
#
# This is a small, memory-conscious pilot for the PACS paper-style
# reconstruction route. It samples cells and peaks, streams the MatrixMarket
# coordinate file, builds a sparse cell x peak matrix, and computes a new
# TF-IDF/LSI/UMAP embedding. It does not use the precomputed GEO UMAP columns,
# does not run PACS, and does not perform full all-feature UMAP reconstruction.

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
  n_cells_per_sample = 500L,
  n_peaks = 20000L,
  seed = 1L,
  chunk_lines = 1000000L,
  progress_every = 5000000L
))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package data.table is required for streaming MatrixMarket chunks.")
}
if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Package irlba is required for pilot LSI.")
}
if (!requireNamespace("uwot", quietly = TRUE)) {
  stop("Package uwot is required for pilot UMAP.")
}

matrix_file <- file.path(params$gse_dir, "GSE157079_snATAC_cell_by_peak_matrix.txt.gz")
peak_file <- file.path(params$gse_dir, "GSE157079_snATAC_peak_list.csv.gz")
if (!file.exists(matrix_file)) stop("Missing matrix file: ", matrix_file)
if (!file.exists(peak_file)) stop("Missing peak list file: ", peak_file)
if (!file.exists(params$metadata_csv)) stop("Missing metadata CSV: ", params$metadata_csv)

pilot_dir <- file.path(params$out_dir, "gse157079_pilot_matrix_lsi_umap")
dir.create(pilot_dir, recursive = TRUE, showWarnings = FALSE)
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
    n_nonzero = as.integer(dims[[3]])
  )
}

stream_subset_matrix <- function(path, selected_cells, selected_peaks, dims, chunk_lines, progress_every) {
  cell_map <- integer(dims$n_cells)
  cell_map[selected_cells] <- seq_along(selected_cells)
  peak_map <- integer(dims$n_peaks)
  peak_map[selected_peaks] <- seq_along(selected_peaks)

  i_chunks <- list()
  j_chunks <- list()
  x_chunks <- list()
  chunk_id <- 0L
  processed <- 0L
  retained <- 0L
  next_progress <- progress_every

  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  invisible(readLines(con, n = 1, warn = FALSE))
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) stop("Could not find MatrixMarket dimensions line")
    if (!startsWith(line, "%")) break
  }

  log_msg("Starting matrix stream from ", path)
  repeat {
    lines <- readLines(con, n = chunk_lines, warn = FALSE)
    if (length(lines) == 0) break
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0) next
    dt <- data.table::fread(
      text = paste(lines, collapse = "\n"),
      header = FALSE,
      col.names = c("i", "j", "x"),
      showProgress = FALSE
    )
    processed <- processed + nrow(dt)
    keep <- cell_map[dt$i] > 0L & peak_map[dt$j] > 0L
    if (any(keep)) {
      chunk_id <- chunk_id + 1L
      i_chunks[[chunk_id]] <- cell_map[dt$i[keep]]
      j_chunks[[chunk_id]] <- peak_map[dt$j[keep]]
      x_chunks[[chunk_id]] <- dt$x[keep]
      retained <- retained + sum(keep)
    }
    if (processed >= next_progress) {
      log_msg("Processed ", processed, " coordinate lines; retained ", retained, " nonzeros")
      next_progress <- next_progress + progress_every
    }
  }
  log_msg("Finished matrix stream: processed ", processed, " coordinate lines; retained ", retained, " nonzeros")

  if (retained == 0L) {
    stop("No nonzero entries retained. Try increasing n_cells_per_sample or n_peaks.")
  }

  sparseMatrix(
    i = unlist(i_chunks, use.names = FALSE),
    j = unlist(j_chunks, use.names = FALSE),
    x = unlist(x_chunks, use.names = FALSE),
    dims = c(length(selected_cells), length(selected_peaks))
  )
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
  "P0_batch1" = "#FF0000",
  "P0_batch2" = "#0000FF",
  "P21_batch1" = "#00CC00",
  "P56_batch1" = "#FFD400",
  "P56_batch2" = "#A000FF"
)
sample_colors <- c("#FF0000", "#0000FF", "#00CC00", "#FFD400", "#A000FF", "#00FFFF", "#FF7F00")

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
  "PT2" = "#B879E8",
  "stroma1" = "#DD70D6",
  "stroma2" = "#F06AA8"
)
celltype_colors <- c(
  "#66C7FF", "#D8896A", "#C69214", "#B6A000", "#7FB000",
  "#67B83F", "#00A651", "#00A98F", "#00A6B8", "#00A7E1",
  "#139DDF", "#8F9BEF", "#B879E8", "#DD70D6", "#F06AA8",
  "#F4A6C8", "#9ADBC5", "#E6C36A", "#BCA7FF", "#7DD3FC"
)

plot_umap <- function(df, color_col, palette, title, out_base) {
  set.seed(1)
  plot_df <- df[sample(seq_len(nrow(df))), , drop = FALSE]
  p <- ggplot(plot_df, aes(x = lsi_umap_1, y = lsi_umap_2, color = .data[[color_col]])) +
    geom_point(size = 0.75, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = palette, drop = FALSE) +
    coord_equal(expand = TRUE) +
    labs(x = "LSI UMAP 1", y = "LSI UMAP 2", color = color_col, title = title) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.title = element_text(size = 15),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 11),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.key.width = grid::unit(0.45, "cm")
    )
  ggsave(paste0(out_base, ".png"), p, width = 8.5, height = 6.5, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 8.5, height = 6.5)
}

format_table <- function(x, max_rows = 60L) {
  paste(capture.output(print(utils::head(x, max_rows), row.names = FALSE)), collapse = "\n")
}

counts_table <- function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  data.frame(level = names(tab), n = as.integer(tab), stringsAsFactors = FALSE)
}

set.seed(params$seed)
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
header <- read_matrix_header(matrix_file)
if (nrow(metadata) != header$n_cells) {
  stop("Metadata rows ", nrow(metadata), " do not match matrix cells ", header$n_cells)
}
if (nrow(peak_list) != header$n_peaks) {
  stop("Peak list rows ", nrow(peak_list), " do not match matrix peaks ", header$n_peaks)
}
if (params$n_peaks > header$n_peaks) {
  stop("Requested n_peaks=", params$n_peaks, " but matrix only has ", header$n_peaks)
}

log_msg("Sampling cells by sample")
selected_metadata <- do.call(rbind, lapply(split(metadata, metadata$sample), function(df) {
  df[sample(seq_len(nrow(df)), min(params$n_cells_per_sample, nrow(df))), , drop = FALSE]
}))
selected_metadata <- selected_metadata[order(selected_metadata$row_index), , drop = FALSE]
selected_cells <- selected_metadata$row_index

log_msg("Sampling peaks")
selected_peaks <- sort(sample(seq_len(header$n_peaks), params$n_peaks))
selected_peak_info <- peak_list[selected_peaks, , drop = FALSE]

counts <- stream_subset_matrix(
  path = matrix_file,
  selected_cells = selected_cells,
  selected_peaks = selected_peaks,
  dims = header,
  chunk_lines = params$chunk_lines,
  progress_every = params$progress_every
)
counts <- as(counts, "dgCMatrix")
rownames(counts) <- paste0("cell_", selected_metadata$row_index)
colnames(counts) <- selected_peak_info$name
log_msg("Sparse pilot matrix built: ", paste(dim(counts), collapse = " x "), "; nnz=", length(counts@x))

pre_filter_dim <- dim(counts)
pre_filter_nnz <- length(counts@x)
nonempty_cells <- Matrix::rowSums(counts) > 0
nonempty_peaks <- Matrix::colSums(counts) > 0
removed_cells <- sum(!nonempty_cells)
removed_peaks <- sum(!nonempty_peaks)
counts <- counts[nonempty_cells, nonempty_peaks, drop = FALSE]
selected_metadata <- selected_metadata[nonempty_cells, , drop = FALSE]
selected_peak_info <- selected_peak_info[nonempty_peaks, , drop = FALSE]
log_msg("After removing empty rows/columns: ", paste(dim(counts), collapse = " x "), "; nnz=", length(counts@x))

if (nrow(counts) < 20 || ncol(counts) < 20) {
  stop("Pilot subset is too sparse after filtering. Try larger n_cells_per_sample or n_peaks.")
}

counts@x <- rep(1, length(counts@x))
depth <- Matrix::rowSums(counts)
peak_detection <- Matrix::colSums(counts)
selected_metadata$pilot_depth <- as.numeric(depth)

counts_rds <- file.path(pilot_dir, "pilot_counts_sparse.rds")
metadata_out <- file.path(pilot_dir, "pilot_metadata.csv")
peaks_out <- file.path(pilot_dir, "pilot_peak_indices.csv")
saveRDS(counts, counts_rds)
write.csv(selected_metadata, metadata_out, row.names = FALSE)
write.csv(selected_peak_info, peaks_out, row.names = FALSE)

log_msg("Starting TF-IDF")
tf <- Diagonal(x = 1 / as.numeric(depth)) %*% counts
idf <- log(1 + nrow(counts) / as.numeric(peak_detection))
tfidf <- tf %*% Diagonal(x = idf)
tfidf@x <- log1p(tfidf@x * 1e4)

nv <- min(30L, nrow(tfidf) - 1L, ncol(tfidf) - 1L)
if (nv < 2L) stop("Not enough rows/columns for LSI")
log_msg("Starting LSI with irlba nv=", nv)
svd_res <- irlba::irlba(tfidf, nv = nv)
lsi <- svd_res$u %*% diag(svd_res$d, nrow = length(svd_res$d))
colnames(lsi) <- paste0("LSI_", seq_len(ncol(lsi)))
rownames(lsi) <- rownames(counts)

depth_cor <- apply(lsi, 2, function(x) suppressWarnings(cor(x, selected_metadata$pilot_depth, method = "spearman")))
umap_dims <- if (ncol(lsi) >= 3) 2:min(30L, ncol(lsi)) else seq_len(ncol(lsi))
log_msg("Starting UMAP with uwot using LSI dims: ", paste(umap_dims, collapse = ","))
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
selected_metadata$lsi_umap_1 <- umap[, 1]
selected_metadata$lsi_umap_2 <- umap[, 2]
write.csv(selected_metadata, metadata_out, row.names = FALSE)

sample_levels <- sort(unique(selected_metadata$sample))
celltype_levels <- sort(unique(selected_metadata$cell_type))
selected_metadata$sample <- factor(selected_metadata$sample, levels = sample_levels)
selected_metadata$cell_type <- factor(selected_metadata$cell_type, levels = celltype_levels)
sample_palette <- make_palette(sample_levels, sample_palette_preferred, sample_colors)
celltype_palette <- make_palette(celltype_levels, celltype_palette_preferred, celltype_colors)

sample_base <- file.path(params$fig_dir, "gse157079_pilot_matrix_lsi_umap_by_sample")
celltype_base <- file.path(params$fig_dir, "gse157079_pilot_matrix_lsi_umap_by_celltype")
log_msg("Saving UMAP figures")
plot_umap(
  selected_metadata,
  "sample",
  sample_palette,
  "GSE157079 pilot matrix-derived LSI UMAP by sample",
  sample_base
)
plot_umap(
  selected_metadata,
  "cell_type",
  celltype_palette,
  "GSE157079 pilot matrix-derived LSI UMAP by cell type",
  celltype_base
)

lsi_out <- file.path(pilot_dir, "pilot_lsi_embedding.csv")
umap_out <- file.path(pilot_dir, "pilot_lsi_umap_embedding.csv")
write.csv(data.frame(cell_id = rownames(lsi), lsi, check.names = FALSE), lsi_out, row.names = FALSE)
write.csv(
  selected_metadata[, c("row_index", "cell_barcode", "sample", "cell_type", "pilot_depth", "lsi_umap_1", "lsi_umap_2")],
  umap_out,
  row.names = FALSE
)

report_path <- file.path(pilot_dir, "pilot_lsi_umap_report.md")
report <- c(
  "# GSE157079 Pilot Matrix-Derived LSI UMAP Report",
  "",
  "This pilot starts from the GSE157079 MatrixMarket cell-by-peak count matrix.",
  "The precomputed GEO `umap_1`/`umap_2` columns were ignored for embedding construction.",
  "",
  "## Command Arguments",
  "",
  "```text",
  paste(names(params), unlist(params), sep = " = ", collapse = "\n"),
  "```",
  "",
  "## Matrix Header",
  "",
  "```text",
  header$first_line,
  header$dims_line,
  "```",
  "",
  "## Sampling",
  "",
  paste0("- Selected cells before empty filtering: ", length(selected_cells)),
  paste0("- Selected peaks before empty filtering: ", length(selected_peaks)),
  "",
  "## Sparse Matrix",
  "",
  paste0("- Dimensions before empty filtering: ", paste(pre_filter_dim, collapse = " x ")),
  paste0("- Nonzeros before empty filtering: ", pre_filter_nnz),
  paste0("- Removed empty cells: ", removed_cells),
  paste0("- Removed empty peaks: ", removed_peaks),
  paste0("- Dimensions after empty filtering: ", paste(dim(counts), collapse = " x ")),
  paste0("- Nonzeros after empty filtering: ", length(counts@x)),
  "",
  "## Sample Table",
  "",
  "```text",
  format_table(counts_table(selected_metadata$sample)),
  "```",
  "",
  "## Cell Type Table",
  "",
  "```text",
  format_table(counts_table(selected_metadata$cell_type)),
  "```",
  "",
  "## Peak Detection Summary",
  "",
  "```text",
  paste(capture.output(print(summary(as.numeric(peak_detection)))), collapse = "\n"),
  "```",
  "",
  "## TF-IDF / LSI",
  "",
  paste0("- TF-IDF dimensions: ", paste(dim(tfidf), collapse = " x ")),
  paste0("- LSI dimensions: ", paste(dim(lsi), collapse = " x ")),
  paste0("- LSI dimensions used for UMAP: ", paste(umap_dims, collapse = ", ")),
  "- Spearman correlations between LSI components and pilot depth:",
  "",
  "```text",
  paste(capture.output(print(depth_cor)), collapse = "\n"),
  "```",
  "",
  "## UMAP",
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
  paste0("- `", peaks_out, "`"),
  paste0("- `", lsi_out, "`"),
  paste0("- `", umap_out, "`"),
  paste0("- `", sample_base, ".png`"),
  paste0("- `", sample_base, ".pdf`"),
  paste0("- `", celltype_base, ".png`"),
  paste0("- `", celltype_base, ".pdf`"),
  "",
  "## Conclusion",
  "",
  "Pilot matrix-derived UMAP succeeded. This is not the final PACS paper figure;",
  "it only verifies that the large cell-by-peak matrix can be streamed into a",
  "sparse pilot subset and used for TF-IDF/LSI/UMAP without relying on the",
  "precomputed GEO UMAP coordinates."
)
writeLines(report, report_path)

log_msg("Saved report: ", report_path)
log_msg("Finished pilot matrix-derived LSI UMAP")
cat("Pilot matrix stream completed: TRUE\n")
cat("Retained sparse matrix dimensions: ", paste(dim(counts), collapse = " x "), "\n", sep = "")
cat("Retained sparse matrix nnz: ", length(counts@x), "\n", sep = "")
cat("irlba available: TRUE\n")
cat("uwot available: TRUE\n")
cat("Sample levels: ", paste(sample_levels, collapse = ", "), "\n", sep = "")
cat("Cell type levels: ", paste(celltype_levels, collapse = ", "), "\n", sep = "")
