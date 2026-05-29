#!/usr/bin/env Rscript

# Lightweight alignment smoke test for the GSE157079 mouse kidney matrix.
#
# This script is intentionally conservative:
# - it does not run UMAP;
# - it does not run PACS;
# - it does not call Matrix::readMM on the full 469 MB gzipped matrix;
# - it only reads the MatrixMarket header and a small prefix of coordinate rows.

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
  gse_dir = "/home/woodson/biostatistic/pacs/GSE157079",
  metadata_csv = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures/gse157079_metadata_merged.csv",
  out_dir = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures"
))

expected <- list(
  n_cells = 28316L,
  n_peaks = 300755L,
  n_nonzero = 166121193L
)

matrix_path <- file.path(params$gse_dir, "GSE157079_snATAC_cell_by_peak_matrix.txt.gz")
peak_path <- file.path(params$gse_dir, "GSE157079_snATAC_peak_list.csv.gz")
required_metadata_cols <- c("row_index", "cell_barcode", "sample", "cell_type", "umap_1", "umap_2")
required_peak_cols <- c("seqnames", "start", "end", "name")

if (!file.exists(matrix_path)) stop("Missing matrix file: ", matrix_path)
if (!file.exists(params$metadata_csv)) stop("Missing merged metadata CSV: ", params$metadata_csv)
if (!file.exists(peak_path)) stop("Missing peak list file: ", peak_path)

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
report_path <- file.path(params$out_dir, "gse157079_matrix_alignment_smoke_test.md")

read_csv_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(
      path,
      header = TRUE,
      showProgress = FALSE,
      check.names = FALSE
    )))
  }
  read.csv(gzfile(path), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

read_plain_csv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(
      path,
      header = TRUE,
      showProgress = FALSE,
      check.names = FALSE
    )))
  }
  read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

read_matrix_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  first_line <- readLines(con, n = 1, warn = FALSE)
  if (length(first_line) != 1) stop("Matrix file is empty: ", path)

  skipped_comments <- character()
  dims_line <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "%")) {
      skipped_comments <- c(skipped_comments, line)
      next
    }
    dims_line <- line
    break
  }
  if (length(dims_line) != 1) stop("Could not find MatrixMarket dimensions line")
  dims <- scan(text = dims_line, quiet = TRUE)
  if (length(dims) < 3) stop("MatrixMarket dimensions line must contain three integers: ", dims_line)
  list(
    first_line = first_line,
    comment_lines = skipped_comments,
    dims_line = dims_line,
    n_cells = as.integer(dims[[1]]),
    n_peaks = as.integer(dims[[2]]),
    n_nonzero = as.integer(dims[[3]])
  )
}

read_coordinate_prefix <- function(path, n = 10000L) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  first_line <- readLines(con, n = 1, warn = FALSE)
  if (length(first_line) != 1) stop("Matrix file is empty: ", path)

  dims_line <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) stop("Could not find MatrixMarket dimensions line before coordinates")
    if (startsWith(line, "%")) next
    dims_line <- line
    break
  }

  coord_lines <- character()
  while (length(coord_lines) < n) {
    lines <- readLines(con, n = n - length(coord_lines), warn = FALSE)
    if (length(lines) == 0) break
    lines <- lines[nzchar(lines)]
    lines <- lines[!startsWith(lines, "%")]
    coord_lines <- c(coord_lines, lines)
  }
  if (length(coord_lines) == 0) stop("No coordinate lines could be read from matrix file")

  parts <- strsplit(coord_lines, "\\s+")
  widths <- vapply(parts, length, integer(1))
  if (any(widths < 3)) {
    stop("At least one coordinate line has fewer than three fields")
  }
  coord <- data.frame(
    row = as.integer(vapply(parts, `[`, character(1), 1)),
    col = as.integer(vapply(parts, `[`, character(1), 2)),
    value = as.numeric(vapply(parts, `[`, character(1), 3))
  )
  list(dims_line = dims_line, n_lines = length(coord_lines), coord = coord)
}

standardize_peak_list <- function(df) {
  if (!all(required_peak_cols %in% names(df))) {
    expected_names <- c("peak_index", "seqnames", "start", "end", "width", "strand", "name")
    if (ncol(df) < length(expected_names)) {
      stop(
        "Peak list lacks required columns and has too few columns to standardize. Names: ",
        paste(names(df), collapse = ", ")
      )
    }
    df <- df[, seq_along(expected_names), drop = FALSE]
    names(df) <- expected_names
  } else if (!"peak_index" %in% names(df)) {
    first_name <- names(df)[[1]]
    if (is.na(first_name) || first_name == "" || grepl("^V1$|^\\.\\.\\.1$|^X$", first_name)) {
      names(df)[[1]] <- "peak_index"
    } else {
      df$peak_index <- seq_len(nrow(df))
    }
  }
  df
}

format_table <- function(x, max_rows = 50L) {
  x <- utils::head(x, max_rows)
  paste(capture.output(print(x, row.names = FALSE)), collapse = "\n")
}

format_counts <- function(x, label, max_rows = 100L) {
  tab <- sort(table(x), decreasing = TRUE)
  df <- data.frame(
    level = names(tab),
    n = as.integer(tab),
    stringsAsFactors = FALSE
  )
  paste0("\n### ", label, "\n\n```text\n", format_table(df, max_rows), "\n```\n")
}

report_store <- local({
  lines <- character()
  list(
    add = function(...) {
      text <- paste0(...)
      lines <<- c(lines, text)
      invisible(lines)
    },
    get = function() lines
  )
})
add <- report_store$add

matrix_header <- read_matrix_header(matrix_path)
metadata <- read_plain_csv(params$metadata_csv)
peak_list <- standardize_peak_list(read_csv_table(peak_path))
coord_prefix <- read_coordinate_prefix(matrix_path, n = 10000L)
coord <- coord_prefix$coord

metadata_missing <- setdiff(required_metadata_cols, names(metadata))
peak_missing <- setdiff(required_peak_cols, names(peak_list))

matrix_header_ok <- grepl("^%%MatrixMarket", matrix_header$first_line)
matrix_dims_ok <- identical(matrix_header$n_cells, expected$n_cells) &&
  identical(matrix_header$n_peaks, expected$n_peaks) &&
  identical(matrix_header$n_nonzero, expected$n_nonzero)
metadata_cols_ok <- length(metadata_missing) == 0
metadata_n_ok <- nrow(metadata) == expected$n_cells
metadata_index_ok <- metadata_cols_ok &&
  !anyNA(metadata$row_index) &&
  identical(sort(as.integer(metadata$row_index)), seq_len(expected$n_cells))
peak_cols_ok <- length(peak_missing) == 0
peak_n_ok <- nrow(peak_list) == expected$n_peaks

peak_index_ok <- TRUE
peak_index_note <- "Peak rows correspond to 300755 peaks by row order."
if ("peak_index" %in% names(peak_list)) {
  peak_index <- suppressWarnings(as.integer(peak_list$peak_index))
  peak_index_ok <- !anyNA(peak_index) && identical(sort(peak_index), seq_len(expected$n_peaks))
  peak_index_note <- if (peak_index_ok) {
    "peak_index is non-NA and covers 1:300755."
  } else {
    "peak_index is not exactly 1:300755; row order still has 300755 peaks."
  }
}

coord_bounds_ok <- all(
  !is.na(coord$row),
  !is.na(coord$col),
  coord$row >= 1,
  coord$row <= expected$n_cells,
  coord$col >= 1,
  coord$col <= expected$n_peaks
)
value_summary <- summary(coord$value)
matrix_orientation_note <- paste0(
  "Coordinate rows range from ", min(coord$row), " to ", max(coord$row),
  " and coordinate columns range from ", min(coord$col), " to ", max(coord$col),
  " in the first ", nrow(coord), " entries. Given the header dimensions ",
  matrix_header$n_cells, " x ", matrix_header$n_peaks,
  ", this is consistent with a cell x peak sparse matrix."
)

sample_levels <- if (metadata_cols_ok) sort(unique(as.character(metadata$sample))) else character()
cell_type_levels <- if (metadata_cols_ok) sort(unique(as.character(metadata$cell_type))) else character()
sample_ok <- length(sample_levels) > 0 && !anyNA(sample_levels) && all(nzchar(sample_levels))
cell_type_ok <- length(cell_type_levels) > 0 && !anyNA(cell_type_levels) && all(nzchar(cell_type_levels))

smoke_passed <- all(
  matrix_header_ok,
  matrix_dims_ok,
  metadata_cols_ok,
  metadata_n_ok,
  metadata_index_ok,
  peak_cols_ok,
  peak_n_ok,
  peak_index_ok,
  coord_bounds_ok,
  sample_ok,
  cell_type_ok
)

add("# GSE157079 Matrix Alignment Smoke Test")
add("")
add("This report checks file alignment for the PACS paper-style reconstruction route.")
add("No UMAP, PACS, or dense matrix materialization was performed.")
add("")
add("## Inputs")
add("")
add("- Matrix: `", matrix_path, "`")
add("- Merged metadata: `", params$metadata_csv, "`")
add("- Peak list: `", peak_path, "`")
add("")
add("## MatrixMarket Header")
add("")
add("```text")
add(matrix_header$first_line)
if (length(matrix_header$comment_lines) > 0) {
  add(paste(matrix_header$comment_lines, collapse = "\n"))
}
add(matrix_header$dims_line)
add("```")
add("")
add("- Header is MatrixMarket-like: `", matrix_header_ok, "`")
add("- Parsed dimensions: cells=", matrix_header$n_cells, ", peaks=", matrix_header$n_peaks, ", nonzero=", matrix_header$n_nonzero)
add("- Expected dimensions: cells=", expected$n_cells, ", peaks=", expected$n_peaks, ", nonzero=", expected$n_nonzero)
add("- Dimensions match expected: `", matrix_dims_ok, "`")
add("")
add("## Metadata")
add("")
add("- Rows: ", nrow(metadata))
add("- Columns: ", paste(names(metadata), collapse = ", "))
add("- Required columns present: `", metadata_cols_ok, "`")
if (!metadata_cols_ok) add("- Missing metadata columns: ", paste(metadata_missing, collapse = ", "))
add("- Row count matches matrix cells: `", metadata_n_ok, "`")
add("- `row_index` covers 1:28316: `", metadata_index_ok, "`")
add("")
add("First metadata rows:")
add("")
add("```text")
add(format_table(utils::head(metadata, 5)))
add("```")
add(format_counts(metadata$sample, "Sample Table"))
add(format_counts(metadata$cell_type, "Cell Type Table"))
add("")
add("## Peak List")
add("")
add("- Rows: ", nrow(peak_list))
add("- Columns: ", paste(names(peak_list), collapse = ", "))
add("- Required columns present: `", peak_cols_ok, "`")
if (!peak_cols_ok) add("- Missing peak columns: ", paste(peak_missing, collapse = ", "))
add("- Row count matches matrix peaks: `", peak_n_ok, "`")
add("- Peak index check: `", peak_index_ok, "`")
add("- Peak index note: ", peak_index_note)
add("")
add("First peak rows:")
add("")
add("```text")
add(format_table(utils::head(peak_list[, intersect(c("peak_index", "seqnames", "start", "end", "width", "strand", "name"), names(peak_list)), drop = FALSE], 5)))
add("```")
add("")
add("## Coordinate-Line Sanity Check")
add("")
add("- Coordinate lines parsed: ", nrow(coord))
add("- Row index min/max: ", min(coord$row), " / ", max(coord$row))
add("- Column index min/max: ", min(coord$col), " / ", max(coord$col))
add("- Coordinate bounds within expected dimensions: `", coord_bounds_ok, "`")
add("- Value summary:")
add("")
add("```text")
add(paste(capture.output(print(value_summary)), collapse = "\n"))
add("```")
add("")
add(matrix_orientation_note)
add("")
add("## Conclusion")
add("")
add("- Smoke test passed: `", smoke_passed, "`")
if (smoke_passed) {
  add("- The GSE157079 matrix, merged metadata, and peak list appear sufficient and aligned for the next PACS paper-style pipeline.")
  add("- The matrix orientation appears to be cell x peak, matching the file name and header dimensions.")
} else {
  add("- The files need review before full PACS paper-style reconstruction.")
}
add("")
add("Recommended next step: write a small sparse-loading prototype that creates a sparse Matrix object from the MatrixMarket file in a controlled output directory, verifies row/column names from metadata and peak list, and then plans all-feature versus PACS-filtered UMAP without dense materialization.")

report_lines <- report_store$get()
writeLines(report_lines, report_path)

cat("Saved smoke-test report: ", normalizePath(report_path, mustWork = FALSE), "\n", sep = "")
cat("Smoke test passed: ", smoke_passed, "\n", sep = "")
cat("Matrix dimensions: ", matrix_header$n_cells, " x ", matrix_header$n_peaks, " nnz=", matrix_header$n_nonzero, "\n", sep = "")
cat("Metadata rows: ", nrow(metadata), "\n", sep = "")
cat("Peak list rows: ", nrow(peak_list), "\n", sep = "")
cat("Sample levels: ", paste(sample_levels, collapse = ", "), "\n", sep = "")
cat("Cell type levels: ", paste(cell_type_levels, collapse = ", "), "\n", sep = "")
cat("Matrix appears cell x peak: ", coord_bounds_ok && matrix_header$n_cells < matrix_header$n_peaks, "\n", sep = "")

if (!smoke_passed) {
  quit(status = 1)
}
