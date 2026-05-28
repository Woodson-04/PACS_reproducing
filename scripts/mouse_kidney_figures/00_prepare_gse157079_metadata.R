#!/usr/bin/env Rscript

# Prepare lightweight GSE157079 metadata/UMAP tables. This script deliberately
# avoids reading the large cell-by-peak matrix.
#
# GEO CSV files have an unnamed first column. For this dataset:
# - UMAP:     row_index, umap_1, umap_2
# - metadata: row_index, cell_barcode, sample, cell_type
# The UMAP table has no barcode, so row_index is the correct merge key.

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
  expected_cells = "28316"
))
expected_cells <- as.integer(params$expected_cells)

paths <- list(
  umap = file.path(params$gse_dir, "GSE157079_snATAC_UMAP_coordinates.csv.gz"),
  metadata = file.path(params$gse_dir, "GSE157079_snATAC_metadata.csv.gz"),
  peak_list = file.path(params$gse_dir, "GSE157079_snATAC_peak_list.csv.gz")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) {
  stop("Missing required GSE157079 files: ", paste(missing, collapse = ", "))
}

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
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

cell_value <- function(df, row, col) {
  if (nrow(df) < row || ncol(df) < col) return(NA_character_)
  tolower(as.character(df[[col]][row]))
}

diagnose_index_problem <- function(df, index, label, header_like = FALSE) {
  na_rows <- which(is.na(index))
  dup_values <- unique(index[duplicated(index) & !is.na(index)])
  first_rows <- paste(capture.output(print(utils::head(df, 5))), collapse = " | ")
  paste0(
    label, " row_index contains NA or duplicate values. ",
    "nrow=", nrow(df), "; names=", paste(names(df), collapse = ", "),
    "; NA row positions=", if (length(na_rows)) paste(head(na_rows, 20), collapse = ", ") else "none",
    "; duplicated row_index values=", if (length(dup_values)) paste(head(dup_values, 20), collapse = ", ") else "none",
    "; first row looks header-like=", header_like,
    "; first rows: ", first_rows
  )
}

standardize_umap <- function(df) {
  if (ncol(df) < 3) {
    stop("UMAP file must have at least 3 columns: row_index, umap_1, umap_2")
  }
  header_like <- cell_value(df, 1, 2) %in% c("umap-1", "umap.1", "umap_1") ||
    cell_value(df, 1, 3) %in% c("umap-2", "umap.2", "umap_2")
  if (header_like) {
    df <- df[-1, , drop = FALSE]
  }
  df <- df[, seq_len(3), drop = FALSE]
  names(df) <- c("row_index", "umap_1", "umap_2")
  df$row_index <- as.integer(df$row_index)
  df$umap_1 <- as.numeric(df$umap_1)
  df$umap_2 <- as.numeric(df$umap_2)
  if (anyNA(df$row_index) || anyDuplicated(df$row_index)) {
    stop(diagnose_index_problem(df, df$row_index, "UMAP", header_like))
  }
  df
}

standardize_metadata <- function(df) {
  if (ncol(df) < 4) {
    stop("Metadata file must have at least 4 columns: row_index, barcodes, samples, clusters")
  }
  header_like <- cell_value(df, 1, 2) %in% c("barcodes", "barcode") ||
    cell_value(df, 1, 3) %in% c("samples", "sample") ||
    cell_value(df, 1, 4) %in% c("clusters", "cluster")
  if (header_like) {
    df <- df[-1, , drop = FALSE]
  }
  df <- df[, seq_len(4), drop = FALSE]
  names(df) <- c("row_index", "cell_barcode", "sample", "cell_type")
  df$row_index <- as.integer(df$row_index)
  df$cell_barcode <- as.character(df$cell_barcode)
  df$sample <- as.character(df$sample)
  df$cell_type <- as.character(df$cell_type)
  if (anyNA(df$row_index) || anyDuplicated(df$row_index)) {
    stop(diagnose_index_problem(df, df$row_index, "Metadata", header_like))
  }
  if (anyNA(df$cell_barcode) || any(df$cell_barcode == "")) {
    stop("Metadata cell_barcode contains NA or empty values")
  }

  sample_barcode_id <- paste(df$sample, df$cell_barcode, sep = "__")
  if (anyDuplicated(sample_barcode_id)) {
    dup_pairs <- unique(sample_barcode_id[duplicated(sample_barcode_id)])
    stop(
      "Metadata sample + cell_barcode contains duplicate values. Example duplicated pairs: ",
      paste(head(dup_pairs, 20), collapse = ", ")
    )
  }

  dup_barcode_n <- sum(duplicated(df$cell_barcode))
  if (dup_barcode_n > 0) {
    message(
      "Note: cell_barcode is not globally unique across samples; ",
      dup_barcode_n,
      " duplicated barcode rows found. This is expected for 10x-style barcodes. ",
      "Using row_index as merge key and sample + cell_barcode as uniqueness check."
    )
  }
  df
}

standardize_peak_list <- function(df) {
  expected <- c("peak_index", "seqnames", "start", "end", "width", "strand", "name")
  if (ncol(df) < length(expected)) {
    stop("Peak list file must have at least columns: ", paste(expected, collapse = ", "))
  }
  header_like <- cell_value(df, 1, 2) %in% c("seqnames", "seqname", "chrom", "chr") ||
    cell_value(df, 1, 3) == "start" ||
    cell_value(df, 1, 4) == "end" ||
    cell_value(df, 1, 7) == "name"
  if (header_like) {
    df <- df[-1, , drop = FALSE]
  }
  df <- df[, seq_along(expected), drop = FALSE]
  names(df) <- expected
  df$peak_index <- as.integer(df$peak_index)
  df
}

umap_raw <- read_table(paths$umap)
metadata_raw <- read_table(paths$metadata)
peak_raw <- read_table(paths$peak_list)

umap <- standardize_umap(umap_raw)
metadata <- standardize_metadata(metadata_raw)
peak_list <- standardize_peak_list(peak_raw)

if (nrow(umap) != nrow(metadata)) {
  stop(
    "UMAP and metadata row counts differ: UMAP=", nrow(umap),
    ", metadata=", nrow(metadata),
    ". Expected both to describe the same cells."
  )
}
if (!is.na(expected_cells) && expected_cells > 0L && nrow(umap) != expected_cells) {
  stop(
    "Unexpected GSE157079 cell count: UMAP=", nrow(umap),
    ", metadata=", nrow(metadata),
    ", expected=", expected_cells
  )
}
if (!identical(sort(umap$row_index), sort(metadata$row_index))) {
  stop("UMAP and metadata row_index values do not match")
}

merged <- merge(metadata, umap, by = "row_index", all = FALSE, sort = TRUE)
merged <- merged[, c("row_index", "cell_barcode", "sample", "cell_type", "umap_1", "umap_2")]

if (nrow(merged) != nrow(metadata)) {
  stop(
    "Merged metadata row count ", nrow(merged),
    " does not equal metadata row count ", nrow(metadata)
  )
}

metadata_out <- file.path(params$out_dir, "gse157079_metadata_merged.csv")
peak_out <- file.path(params$out_dir, "gse157079_peak_list_preview.csv")
summary_out <- file.path(params$out_dir, "gse157079_metadata_summary.csv")

write.csv(merged, metadata_out, row.names = FALSE)
write.csv(head(peak_list, 5000), peak_out, row.names = FALSE)

summary_df <- rbind(
  data.frame(
    summary_type = "merge",
    field = "row_index",
    value = "metadata_umap_overlap",
    n = nrow(merged),
    stringsAsFactors = FALSE
  ),
  data.frame(
    summary_type = "cell_type_counts",
    field = "cell_type",
    value = names(sort(table(merged$cell_type), decreasing = TRUE)),
    n = as.integer(sort(table(merged$cell_type), decreasing = TRUE)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    summary_type = "sample_counts",
    field = "sample",
    value = names(sort(table(merged$sample), decreasing = TRUE)),
    n = as.integer(sort(table(merged$sample), decreasing = TRUE)),
    stringsAsFactors = FALSE
  )
)
write.csv(summary_df, summary_out, row.names = FALSE)

cat("Merged metadata rows: ", nrow(merged), "\n", sep = "")
cat("Merged metadata columns: ", paste(names(merged), collapse = ", "), "\n", sep = "")
cat("Peak preview columns: ", paste(names(peak_list), collapse = ", "), "\n", sep = "")
cat("Saved: ", normalizePath(metadata_out, mustWork = FALSE), "\n", sep = "")
cat("Saved: ", normalizePath(peak_out, mustWork = FALSE), "\n", sep = "")
cat("Saved: ", normalizePath(summary_out, mustWork = FALSE), "\n", sep = "")
