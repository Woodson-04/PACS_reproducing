#!/usr/bin/env Rscript

# All-cell, all-feature matrix-derived UMAP for GSE157079.
#
# This script uses all 300755 peaks from the GSE157079 cell-by-peak
# MatrixMarket file. It does not use GEO precomputed UMAP coordinates, does not
# run PACS, and does not remove batch-effect features.

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
    if (is.logical(old)) {
      defaults[[name]] <- tolower(value) %in% c("true", "t", "1", "yes", "y")
    } else if (is.integer(old)) {
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
  seed = 1L,
  chunk_lines = 100000L,
  progress_every = 5000000L,
  n_lsi = 50L,
  umap_lsi_start = 2L,
  umap_lsi_end = 30L,
  save_counts = FALSE,
  dry_run = FALSE,
  matrix_file = ""
))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package data.table is required for streaming MatrixMarket chunks.")
}

if (!params$dry_run) {
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package irlba is required for all-feature LSI.")
  }
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package uwot is required for all-feature UMAP.")
  }
}

if (!nzchar(params$matrix_file)) {
  params$matrix_file <- file.path(params$gse_dir, "GSE157079_snATAC_cell_by_peak_matrix.txt.gz")
}
peak_file <- file.path(params$gse_dir, "GSE157079_snATAC_peak_list.csv.gz")
if (!file.exists(params$matrix_file)) stop("Missing matrix file: ", params$matrix_file)
if (!file.exists(peak_file)) stop("Missing peak list file: ", peak_file)
if (!file.exists(params$metadata_csv)) stop("Missing metadata CSV: ", params$metadata_csv)

run_dir <- file.path(params$out_dir, "gse157079_all_cells_all_features_lsi_umap")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(run_dir, "all_features_lsi_umap_run.log")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
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
    df$peak_index <- seq_len(nrow(df))
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

format_table <- function(x, max_rows = 80L) {
  paste(capture.output(print(utils::head(x, max_rows), row.names = FALSE)), collapse = "\n")
}

counts_table <- function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  data.frame(level = names(tab), n = as.integer(tab), stringsAsFactors = FALSE)
}

object_size_mb <- function(x) {
  round(as.numeric(utils::object.size(x)) / 1024^2, 2)
}

estimate_dgc_memory <- function(nnz, ncol) {
  x_bytes <- 8 * nnz
  i_bytes <- 4 * nnz
  p_bytes <- 4 * (ncol + 1)
  raw_bytes <- x_bytes + i_bytes + p_bytes
  data.frame(
    component = c("x numeric", "i row indices", "p column pointers", "raw dgCMatrix core", "rough with overhead"),
    bytes = c(x_bytes, i_bytes, p_bytes, raw_bytes, raw_bytes * 1.5),
    gigabytes = round(c(x_bytes, i_bytes, p_bytes, raw_bytes, raw_bytes * 1.5) / 1024^3, 3),
    stringsAsFactors = FALSE
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
    geom_point(size = 0.75, alpha = 0.85, stroke = 0) +
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

write_feasibility_report <- function(path, params, header, metadata, peak_list, memory_estimate, aligned) {
  report <- c(
    "# GSE157079 All-Feature LSI UMAP Feasibility Report",
    "",
    "This dry run checks whether an all-cell, all-feature matrix-derived UMAP is feasible.",
    "No full matrix was streamed, no LSI/UMAP was run, and no dense matrix was created.",
    "",
    "## Conceptual Correction",
    "",
    "- Previous all-cell UMAP used top 20000 detected peaks and should be called top-peak UMAP.",
    "- This target all-feature UMAP uses all 300755 peaks from the GSE157079 cell-by-peak matrix.",
    "",
    "## Arguments",
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
    "## Alignment Checks",
    "",
    paste0("- metadata rows: ", nrow(metadata)),
    paste0("- peak list rows: ", nrow(peak_list)),
    paste0("- row_index covers matrix cells: ", aligned$row_index_ok),
    paste0("- peak_index covers matrix peaks: ", aligned$peak_index_ok),
    paste0("- metadata rows match matrix cells: ", aligned$metadata_n_ok),
    paste0("- peak rows match matrix peaks: ", aligned$peak_n_ok),
    "",
    "## Estimated dgCMatrix Memory",
    "",
    "```text",
    format_table(memory_estimate),
    "```",
    "",
    "## Feasibility Interpretation",
    "",
    "The raw dgCMatrix core is expected to require roughly the memory shown above.",
    "The full run will require additional memory for coordinate chunks, sparse operations, TF-IDF, irlba workspace, and UMAP.",
    "If available RAM is limited, the full all-feature run may fail during matrix construction, TF-IDF, or LSI.",
    "",
    "## Conclusion",
    "",
    paste0("- dry-run inputs aligned: ", all(unlist(aligned))),
    "- full all-feature UMAP should only be attempted on a Linux session with sufficient RAM.",
    "- If memory fails, keep the previous top-20000 UMAP wording and state that all-feature UMAP was attempted but was computationally limited."
  )
  writeLines(report, path)
}

write_failure_report <- function(path, params, error_msg, stage = "unknown") {
  report <- c(
    "# GSE157079 All-Feature LSI UMAP Failure Report",
    "",
    "The all-feature UMAP run did not complete.",
    "",
    "## Stage",
    "",
    stage,
    "",
    "## Error",
    "",
    "```text",
    error_msg,
    "```",
    "",
    "## Arguments",
    "",
    "```text",
    paste(names(params), unlist(params), sep = " = ", collapse = "\n"),
    "```",
    "",
    "## Interpretation",
    "",
    "Do not claim that an all-feature UMAP was generated.",
    "Keep the previous top-20000 result wording unless a later full all-feature run succeeds."
  )
  writeLines(report, path)
}

build_all_feature_matrix <- function(path, dims, chunk_lines, progress_every) {
  i_chunks <- list()
  j_chunks <- list()
  chunk_id <- 0L
  processed <- 0
  next_progress <- progress_every
  con <- open_coordinate_stream(path)
  on.exit(close(con), add = TRUE)
  log_msg("Streaming all MatrixMarket coordinate entries")
  repeat {
    dt <- read_coordinate_chunk(con, chunk_lines)
    if (is.null(dt)) break
    if (nrow(dt) == 0) next
    processed <- processed + nrow(dt)
    chunk_id <- chunk_id + 1L
    i_chunks[[chunk_id]] <- as.integer(dt$i)
    j_chunks[[chunk_id]] <- as.integer(dt$j)
    if (processed >= next_progress) {
      log_msg("Processed ", processed, " coordinate lines")
      next_progress <- next_progress + progress_every
    }
  }
  log_msg("Coordinate stream completed: processed ", processed, " entries")
  if (processed != dims$n_nonzero) {
    log_msg("WARNING: processed nonzero count ", processed, " differs from header ", dims$n_nonzero)
  }
  log_msg("Constructing sparse Matrix")
  i_all <- unlist(i_chunks, use.names = FALSE)
  j_all <- unlist(j_chunks, use.names = FALSE)
  mat <- sparseMatrix(
    i = i_all,
    j = j_all,
    x = rep(1, length(i_all)),
    dims = c(dims$n_cells, dims$n_peaks)
  )
  rm(i_chunks, j_chunks, i_all, j_all)
  gc(verbose = FALSE)
  list(matrix = as(mat, "dgCMatrix"), processed = processed, retained = length(mat@x))
}

main <- function() {
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
  aligned <- list(
    metadata_n_ok = nrow(metadata) == header$n_cells,
    peak_n_ok = nrow(peak_list) == header$n_peaks,
    row_index_ok = identical(sort(metadata$row_index), seq_len(header$n_cells)),
    peak_index_ok = identical(sort(peak_list$peak_index), seq_len(header$n_peaks))
  )
  if (!all(unlist(aligned))) {
    stop("Input alignment failed: ", paste(names(aligned), unlist(aligned), sep = "=", collapse = "; "))
  }

  memory_estimate <- estimate_dgc_memory(header$n_nonzero, header$n_peaks)
  if (params$dry_run) {
    dry_path <- file.path(run_dir, "all_features_feasibility_report.md")
    write_feasibility_report(dry_path, params, header, metadata, peak_list, memory_estimate, aligned)
    log_msg("Dry-run report saved: ", dry_path)
    print(memory_estimate)
    return(invisible(NULL))
  }

  log_msg("All-feature full run start")
  log_msg("Expected matrix dim: ", header$n_cells, " x ", header$n_peaks, "; nnz=", header$n_nonzero)
  log_msg("Estimated raw dgCMatrix core memory GB: ", memory_estimate$gigabytes[memory_estimate$component == "raw dgCMatrix core"])

  pass <- build_all_feature_matrix(params$matrix_file, header, params$chunk_lines, params$progress_every)
  counts <- pass$matrix
  log_msg("Sparse matrix built: ", paste(dim(counts), collapse = " x "), "; nnz=", length(counts@x), "; size MB=", object_size_mb(counts))
  rownames(counts) <- paste0("cell_", seq_len(nrow(counts)))
  colnames(counts) <- make.unique(as.character(peak_list$name))

  pre_filter_dim <- dim(counts)
  pre_filter_nnz <- length(counts@x)
  nonempty_cells <- Matrix::rowSums(counts) > 0
  nonempty_peaks <- Matrix::colSums(counts) > 0
  removed_cells <- sum(!nonempty_cells)
  removed_peaks <- sum(!nonempty_peaks)
  counts <- counts[nonempty_cells, nonempty_peaks, drop = FALSE]
  metadata_filtered <- metadata[nonempty_cells, , drop = FALSE]
  peak_list_filtered <- peak_list[nonempty_peaks, , drop = FALSE]
  log_msg("After empty filtering: ", paste(dim(counts), collapse = " x "), "; nnz=", length(counts@x), "; size MB=", object_size_mb(counts))

  counts@x <- rep(1, length(counts@x))
  all_feature_depth <- Matrix::rowSums(counts)
  peak_detection <- Matrix::colSums(counts)
  metadata_filtered$all_feature_depth <- as.numeric(all_feature_depth)

  if (params$save_counts) {
    counts_rds <- file.path(run_dir, "all_features_counts_sparse.rds")
    saveRDS(counts, counts_rds)
    log_msg("Saved sparse counts RDS: ", counts_rds)
  } else {
    counts_rds <- "not saved (--save_counts FALSE)"
  }
  metadata_out <- file.path(run_dir, "metadata_with_all_feature_umap.csv")
  peak_out <- file.path(run_dir, "all_feature_peak_summary.csv")
  peak_list_filtered$peak_detection <- as.numeric(peak_detection)
  write.csv(peak_list_filtered, peak_out, row.names = FALSE)

  log_msg("TF-IDF start")
  tf <- Diagonal(x = 1 / as.numeric(all_feature_depth)) %*% counts
  idf <- log(1 + nrow(counts) / as.numeric(peak_detection))
  tfidf <- tf %*% Diagonal(x = idf)
  tfidf@x <- log1p(tfidf@x * 1e4)
  log_msg("TF-IDF complete: dim=", paste(dim(tfidf), collapse = " x "), "; nnz=", length(tfidf@x), "; size MB=", object_size_mb(tfidf))

  nv <- min(params$n_lsi, nrow(tfidf) - 1L, ncol(tfidf) - 1L)
  if (nv < 2L) stop("Not enough rows/columns for LSI")
  log_msg("LSI start with irlba nv=", nv)
  svd_res <- irlba::irlba(tfidf, nv = nv)
  lsi <- svd_res$u %*% diag(svd_res$d, nrow = length(svd_res$d))
  colnames(lsi) <- paste0("LSI_", seq_len(ncol(lsi)))
  rownames(lsi) <- rownames(counts)
  log_msg("LSI complete: dim=", paste(dim(lsi), collapse = " x "))

  depth_cor <- apply(lsi, 2, function(x) suppressWarnings(cor(x, metadata_filtered$all_feature_depth, method = "spearman")))
  umap_dims <- params$umap_lsi_start:min(params$umap_lsi_end, ncol(lsi))
  if (length(umap_dims) < 2L) stop("Not enough LSI dimensions for UMAP")
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
  log_msg("UMAP complete")

  metadata_filtered$lsi_umap_1 <- umap[, 1]
  metadata_filtered$lsi_umap_2 <- umap[, 2]
  write.csv(metadata_filtered, metadata_out, row.names = FALSE)

  sample_levels <- sort(unique(metadata_filtered$sample))
  celltype_levels <- sort(unique(metadata_filtered$cell_type))
  metadata_filtered$sample <- factor(metadata_filtered$sample, levels = sample_levels)
  metadata_filtered$cell_type <- factor(metadata_filtered$cell_type, levels = celltype_levels)
  sample_palette <- make_palette(sample_levels, sample_palette_preferred, sample_colors)
  celltype_palette <- make_palette(celltype_levels, celltype_palette_preferred, celltype_colors)

  sample_base <- file.path(params$fig_dir, "gse157079_all_cells_all_features_lsi_umap_by_sample")
  celltype_base <- file.path(params$fig_dir, "gse157079_all_cells_all_features_lsi_umap_by_celltype")
  log_msg("Saving all-feature UMAP figures")
  plot_umap(metadata_filtered, "sample", sample_palette, "GSE157079 all-cell all-feature LSI UMAP by sample", sample_base, params$seed)
  plot_umap(metadata_filtered, "cell_type", celltype_palette, "GSE157079 all-cell all-feature LSI UMAP by cell type", celltype_base, params$seed)

  lsi_out <- file.path(run_dir, "all_features_lsi_embedding.csv")
  umap_out <- file.path(run_dir, "all_features_lsi_umap_embedding.csv")
  write.csv(data.frame(cell_id = rownames(lsi), lsi, check.names = FALSE), lsi_out, row.names = FALSE)
  write.csv(
    metadata_filtered[, c("row_index", "cell_barcode", "sample", "cell_type", "all_feature_depth", "lsi_umap_1", "lsi_umap_2")],
    umap_out,
    row.names = FALSE
  )

  report_path <- file.path(run_dir, "all_features_lsi_umap_report.md")
  report <- c(
    "# GSE157079 All-Cell All-Feature LSI UMAP Report",
    "",
    "This run computes a matrix-derived TF-IDF/LSI/UMAP embedding from all cells and all features.",
    "",
    "## Wording",
    "",
    "- `all-feature` means all 300755 peaks from the GSE157079 cell-by-peak matrix.",
    "- This embedding is unfiltered and does not use PACS feature removal.",
    "- This embedding does not use GEO precomputed UMAP coordinates.",
    "- This is still not the same as the PACS paper adult kidney all-feature UMAP unless sample subset, preprocessing, and parameters match.",
    "",
    "## Arguments",
    "",
    "```text",
    paste(names(params), unlist(params), sep = " = ", collapse = "\n"),
    "```",
    "",
    "## MatrixMarket Header",
    "",
    "```text",
    header$first_line,
    header$dims_line,
    "```",
    "",
    "## Input Summary",
    "",
    paste0("- metadata rows: ", nrow(metadata)),
    paste0("- peak list rows: ", nrow(peak_list)),
    paste0("- processed coordinate lines: ", pass$processed),
    paste0("- retained nonzeros: ", pass$retained),
    "",
    "## Matrix Summary",
    "",
    paste0("- dimensions before empty filtering: ", paste(pre_filter_dim, collapse = " x ")),
    paste0("- nnz before empty filtering: ", pre_filter_nnz),
    paste0("- removed empty cells: ", removed_cells),
    paste0("- removed empty peaks: ", removed_peaks),
    paste0("- dimensions after empty filtering: ", paste(dim(counts), collapse = " x ")),
    paste0("- nnz after empty filtering: ", length(counts@x)),
    paste0("- counts RDS: ", counts_rds),
    "",
    "## TF-IDF / LSI / UMAP",
    "",
    paste0("- TF-IDF dim: ", paste(dim(tfidf), collapse = " x ")),
    paste0("- LSI dim: ", paste(dim(lsi), collapse = " x ")),
    paste0("- UMAP LSI dims: ", paste(umap_dims, collapse = ", ")),
    "",
    "### LSI-depth Spearman Correlations",
    "",
    "```text",
    paste(capture.output(print(round(depth_cor, 4))), collapse = "\n"),
    "```",
    "",
    "## Sample Table",
    "",
    "```text",
    format_table(counts_table(metadata_filtered$sample)),
    "```",
    "",
    "## Cell Type Table",
    "",
    "```text",
    format_table(counts_table(metadata_filtered$cell_type)),
    "```",
    "",
    "## Output Files",
    "",
    paste0("- `", metadata_out, "`"),
    paste0("- `", lsi_out, "`"),
    paste0("- `", umap_out, "`"),
    paste0("- `", peak_out, "`"),
    paste0("- `", sample_base, ".png/pdf`"),
    paste0("- `", celltype_base, ".png/pdf`"),
    "",
    "## Conclusion",
    "",
    "The all-cell all-feature matrix-derived UMAP completed successfully if this report was written."
  )
  writeLines(report, report_path)

  comp_path <- file.path(run_dir, "top20000_vs_allfeatures_umap_comparison.md")
  comp <- c(
    "# Top-20000 vs All-Features UMAP Comparison",
    "",
    "## Top-20000 Previous Result",
    "",
    "- dimensions: 28316 x 20000",
    "- retained nonzeros: 85801336",
    "",
    "## All-Features Result",
    "",
    paste0("- dimensions: ", paste(dim(counts), collapse = " x ")),
    paste0("- retained nonzeros: ", length(counts@x)),
    "",
    "## Interpretation",
    "",
    "- Top-20000 UMAP should be called all-cell top-peak UMAP, not all-feature UMAP.",
    "- All-features UMAP uses all peaks retained after empty filtering from the 300755-peak matrix.",
    "- UMAP global shape need not match exactly between top-peak and all-feature runs.",
    "- If this all-feature run succeeds, report wording can be updated to `All-Feature UMAP`.",
    "- If visual review finds the all-feature embedding less clear, both embeddings can be shown with precise labels."
  )
  writeLines(comp, comp_path)

  log_msg("All-feature run completed")
  log_msg("Report saved: ", report_path)
}

tryCatch(
  main(),
  error = function(e) {
    msg <- conditionMessage(e)
    log_msg("ERROR: ", msg)
    failure_path <- file.path(run_dir, "all_features_failure_report.md")
    write_failure_report(failure_path, params, msg)
    stop(e)
  }
)
