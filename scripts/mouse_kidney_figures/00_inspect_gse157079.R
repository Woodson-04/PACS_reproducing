#!/usr/bin/env Rscript

# Inspect downloaded GSE157079 files without extracting or modifying the
# read-only source directory. The large cell-by-peak matrix is sampled only from
# its header/first lines.

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
  out_dir = "/home/woodson/PACS_reproducing/results/mouse_kidney_figures",
  preview_n = "6"
))
preview_n <- as.integer(params$preview_n)

files <- c(
  umap = "GSE157079_snATAC_UMAP_coordinates.csv.gz",
  metadata = "GSE157079_snATAC_metadata.csv.gz",
  peak_list = "GSE157079_snATAC_peak_list.csv.gz",
  cell_by_peak_matrix = "GSE157079_snATAC_cell_by_peak_matrix.txt.gz"
)

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(params$out_dir, "gse157079_inspection_report.md")

read_preview <- function(path, n = 6) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

split_header <- function(line) {
  if (length(line) == 0 || is.na(line) || !nzchar(line)) return(character())
  if (grepl("\t", line, fixed = TRUE)) {
    strsplit(line, "\t", fixed = TRUE)[[1]]
  } else {
    strsplit(line, ",", fixed = TRUE)[[1]]
  }
}

safe_fread <- function(path, nrows = Inf) {
  if (!requireNamespace("data.table", quietly = TRUE)) return(NULL)
  tryCatch(
    as.data.frame(data.table::fread(
      path,
      header = TRUE,
      nrows = nrows,
      showProgress = FALSE,
      check.names = FALSE
    )),
    error = function(e) NULL
  )
}

cell_value <- function(df, row, col) {
  if (is.null(df) || nrow(df) < row || ncol(df) < col) return(NA_character_)
  tolower(as.character(df[[col]][row]))
}

header_like_first_row <- function(df, nm) {
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  if (nm == "umap") {
    return(cell_value(df, 1, 2) %in% c("umap-1", "umap.1", "umap_1") ||
      cell_value(df, 1, 3) %in% c("umap-2", "umap.2", "umap_2"))
  }
  if (nm == "metadata") {
    return(cell_value(df, 1, 2) %in% c("barcodes", "barcode") ||
      cell_value(df, 1, 3) %in% c("samples", "sample") ||
      cell_value(df, 1, 4) %in% c("clusters", "cluster"))
  }
  if (nm == "peak_list") {
    return(cell_value(df, 1, 2) %in% c("seqnames", "seqname", "chrom", "chr") ||
      cell_value(df, 1, 3) == "start" ||
      cell_value(df, 1, 4) == "end" ||
      cell_value(df, 1, 7) == "name")
  }
  FALSE
}

drop_header_like_first_row <- function(df, nm) {
  detected <- header_like_first_row(df, nm)
  if (detected) {
    df <- df[-1, , drop = FALSE]
  }
  list(data = df, detected = detected)
}

guess_matrix_format <- function(lines) {
  if (length(lines) == 0) return("unknown: empty preview")
  first <- lines[[1]]
  if (startsWith(first, "%%MatrixMarket")) return("MatrixMarket-like")
  delim <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
  fields <- strsplit(first, delim, fixed = TRUE)[[1]]
  if (length(fields) <= 3) return("sparse/triplet-like or coordinate table")
  if (length(lines) >= 2) {
    second_fields <- strsplit(lines[[2]], delim, fixed = TRUE)[[1]]
    if (length(second_fields) == length(fields)) return("dense table-like")
  }
  "other delimited text"
}

find_matching_cols <- function(cols, patterns) {
  cols[grepl(paste(patterns, collapse = "|"), cols, ignore.case = TRUE)]
}

report <- character()
add <- function(...) {
  report <<- c(report, paste0(...))
}

add("# GSE157079 Inspection Report")
add("")
add("Source directory: `", params$gse_dir, "`")
add("")
add("This report was generated without extracting or modifying source `.gz` files.")
add("")

for (nm in names(files)) {
  path <- file.path(params$gse_dir, files[[nm]])
  exists_file <- file.exists(path)
  info <- if (exists_file) file.info(path) else NULL
  add("## ", nm)
  add("")
  add("- File: `", path, "`")
  add("- Exists: ", exists_file)
  if (exists_file) {
    add("- Size bytes: ", info$size)
    lines <- read_preview(path, preview_n)
    cols <- split_header(lines[[1]])
    add("- Preview lines read: ", length(lines))
    add("- Delimited column count from first line: ", length(cols))
    if (length(cols) > 0) {
      add("- First-line fields: `", paste(head(cols, 20), collapse = "`, `"), if (length(cols) > 20) "`, ..." else "`")
    }
    add("")
    add("Preview:")
    add("")
    add("```text")
    add(paste(lines, collapse = "\n"))
    add("```")
    add("")

    if (nm != "cell_by_peak_matrix") {
      dt0 <- safe_fread(path, nrows = 0)
      if (!is.null(dt0)) {
        add("- Parsed column names: `", paste(names(dt0), collapse = "`, `"), "`")
      }
      dt <- safe_fread(path)
      if (!is.null(dt)) {
        clean <- drop_header_like_first_row(dt, nm)
        dt_clean <- clean$data
        add("- Parsed dimensions: ", nrow(dt), " rows x ", ncol(dt), " columns")
        if (clean$detected) {
          add("- Header-like first data row: detected and removed for clean row count")
        } else {
          add("- Header-like first data row: not detected")
        }
        add("- Clean dimensions: ", nrow(dt_clean), " rows x ", ncol(dt_clean), " columns")
        barcode_cols <- find_matching_cols(names(dt), c("barcode", "cell", "cellid", "cell_id"))
        umap_cols <- find_matching_cols(names(dt), c("umap", "UMAP"))
        celltype_cols <- find_matching_cols(names(dt), c("cell.?type", "annotation", "cluster", "identity"))
        peak_coord_cols <- find_matching_cols(names(dt), c("^chr$", "chrom", "start", "end", "peak", "seqnames", "name"))
        add("- Candidate barcode columns: ", if (length(barcode_cols)) paste(barcode_cols, collapse = ", ") else "none")
        add("- Candidate UMAP columns: ", if (length(umap_cols)) paste(umap_cols, collapse = ", ") else "none")
        add("- Candidate cell type/annotation columns: ", if (length(celltype_cols)) paste(celltype_cols, collapse = ", ") else "none")
        add("- Candidate peak coordinate columns: ", if (length(peak_coord_cols)) paste(peak_coord_cols, collapse = ", ") else "none")
      } else {
        add("- Full parse skipped or failed; preview above is still available.")
      }
    } else {
      add("- Matrix format guess: ", guess_matrix_format(lines))
      add("- Full matrix read: skipped intentionally because this file is large.")
    }
  }
  add("")
}

writeLines(report, out_file)
cat("Saved GSE157079 inspection report: ", normalizePath(out_file, mustWork = FALSE), "\n", sep = "")
