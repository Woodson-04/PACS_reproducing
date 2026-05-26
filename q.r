#!/usr/bin/env Rscript

# PACS Notebook 1 reproduction: Type I error and power on real kidney data.
# This script keeps the author's core workflow:
# - original pmats: cells x peaks
# - PACS input: peaks x cells
# - PT-only data for permuted-label Type I error
# - PT+LOH data for actual-label power
# Baseline methods are optional because the current PACS package does not export
# the original notebook's other_methods_for_differential_updated.R functions.

suppressPackageStartupMessages({
  library(Matrix)
  library(PACS)
})

parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    return(defaults)
  }

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    name <- sub("^--", "", key)
    if (!name %in% names(defaults)) {
      stop("Unknown argument: --", name)
    }
    if (i == length(args)) {
      stop("Missing value for argument: --", name)
    }
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
  data_dir = "/home/woodson/biostatistic/pacs/PACS_data",
  output_dir = "results",
  seed = 3384L,
  n_repeat = 5L,
  n_cell_sample = 500L,
  n_features_sample = 10000L,
  count_outlier_cap = 10,
  feature_total_count_cut = 18,
  t_prop_cutoff = 0.2,
  snap_bcv = 0.4,
  run_baselines = FALSE,
  overwrite = FALSE
))

options(error = function() {
  traceback(2)
  q(status = 1)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE))
} else {
  getwd()
}
if (!grepl("^(/|[A-Za-z]:)", params$output_dir)) {
  params$output_dir <- file.path(script_dir, params$output_dir)
}

set.seed(params$seed)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(params$output_dir, paste0("kidney_notebook1_", timestamp))
if (dir.exists(run_dir) && !params$overwrite) {
  stop("Output directory already exists: ", run_dir)
}
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(run_dir, "run.log")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

get_fun <- function(fun_name) {
  if (exists(fun_name, where = asNamespace("PACS"), inherits = FALSE)) {
    return(get(fun_name, envir = asNamespace("PACS")))
  }
  if (exists(fun_name, mode = "function", inherits = TRUE)) {
    return(get(fun_name, mode = "function", inherits = TRUE))
  }
  NULL
}

check_sample_size <- function(available, requested, label) {
  if (requested > available) {
    stop(label, ": requested ", requested, " but only ", available, " available")
  }
}

check_pacs_input <- function(pic_matrix, group_info, cap_rates, label) {
  if (!inherits(pic_matrix, "sparseMatrix")) {
    stop(label, ": pic_matrix must inherit from sparseMatrix")
  }
  if (is.null(rownames(pic_matrix)) || length(rownames(pic_matrix)) != nrow(pic_matrix)) {
    stop(label, ": rownames length must equal nrow(pic_matrix)")
  }
  if (is.null(colnames(pic_matrix)) || length(colnames(pic_matrix)) != ncol(pic_matrix)) {
    stop(label, ": colnames length must equal ncol(pic_matrix)")
  }
  if (nrow(group_info) != ncol(pic_matrix)) {
    stop(label, ": nrow(group.info) must equal ncol(pic_matrix)")
  }
  if (length(cap_rates) != ncol(pic_matrix)) {
    stop(label, ": length(cap_rates) must equal ncol(pic_matrix)")
  }
}

store_pvalues <- function(target, values, expected_len, method, scenario, repeat_id) {
  values <- as.numeric(values)
  if (length(values) != expected_len) {
    msg <- paste0(
      method, " returned ", length(values), " p-values for ", scenario,
      " repeat ", repeat_id, "; expected ", expected_len,
      ". Refusing to store mismatched p-values."
    )
    log_msg("ERROR: ", msg)
    stop(msg)
  }
  target[, repeat_id] <- values
  target
}

pacs_test_sparse_local_fixed <- function(covariate_meta.data,
                                         formula_full,
                                         formula_null,
                                         pic_matrix,
                                         n_peaks_per_round = NULL,
                                         T_proportion_cutoff = 0.25,
                                         cap_rates,
                                         par_initial_null = NULL,
                                         par_initial_full = NULL,
                                         n_cores = 1,
                                         verbose = TRUE,
                                         label = "") {
  log_msg(
    "Calling local fixed PACS sparse wrapper ", label,
    ": dim=", paste(dim(pic_matrix), collapse = " x "),
    ", rownames=", length(rownames(pic_matrix)),
    ", colnames=", length(colnames(pic_matrix)),
    ", meta_rows=", nrow(covariate_meta.data),
    ", cap_rates=", length(cap_rates)
  )

  if (is.null(fn$pacs_test_cumu) || is.null(fn$pacs_test_logit)) {
    stop("pacs_test_sparse_local_fixed requires pacs_test_cumu and pacs_test_logit from library(PACS)")
  }

  n_cell <- ncol(pic_matrix)
  n_peaks <- nrow(pic_matrix)
  p_names <- rownames(pic_matrix)

  if (nrow(covariate_meta.data) != n_cell) {
    stop("number of cells do not match between meta.data and data matrix")
  }
  if (length(cap_rates) != n_cell) {
    stop("number of cells do not match between cap_rates and data matrix")
  }
  if (is.null(p_names)) {
    p_names <- paste("f", seq_len(n_peaks), sep = "_")
    rownames(pic_matrix) <- p_names
  }
  if (is.null(colnames(pic_matrix))) {
    colnames(pic_matrix) <- paste0("cell_", seq_len(n_cell))
  }
  if (is.null(n_peaks_per_round)) {
    n_peaks_per_round <- min(floor(2^29 / n_cell), n_peaks)
  }

  if (inherits(pic_matrix, "Matrix")) {
    pic_matrix_2 <- pic_matrix
    pic_matrix_2@x[pic_matrix_2@x == 1] <- 0
    pic_matrix_2 <- Matrix::drop0(pic_matrix_2)
    pic_matrix_2@x <- rep(1, length = length(pic_matrix_2@x))
    pic_matrixbin <- pic_matrix
    pic_matrixbin@x <- rep(1, length = length(pic_matrixbin@x))
  } else {
    pic_matrix_2 <- Matrix::Matrix(pic_matrix, sparse = TRUE)
    pic_matrix_2@x[pic_matrix_2@x == 1] <- 0
    pic_matrix_2 <- Matrix::drop0(pic_matrix_2)
    pic_matrix_2@x <- rep(1, length = length(pic_matrix_2@x))
    pic_matrixbin <- Matrix::Matrix(pic_matrix, sparse = TRUE)
    pic_matrixbin@x <- rep(1, length = length(pic_matrixbin@x))
  }

  rs <- Matrix::rowSums(pic_matrixbin)
  rs2 <- Matrix::rowSums(pic_matrix_2)
  p_2 <- rs2 / rs
  p_2[is.na(p_2)] <- 0
  n_p_2 <- sum(p_2 >= T_proportion_cutoff)
  n_p_b <- sum(p_2 < T_proportion_cutoff)

  if (verbose) {
    log_msg(label, ": ", n_p_2, " peaks consider cumulative logit models")
    log_msg(label, ": ", n_p_b, " peaks consider logit models")
  }

  f_sel <- names(p_2)[p_2 >= T_proportion_cutoff]
  f_b_sel <- names(p_2)[p_2 < T_proportion_cutoff]
  rm(pic_matrix_2)
  gc(verbose = FALSE)

  p_cumu <- list()
  p_logit <- list()

  if (n_p_2 >= 1) {
    pic_matrix_cumu <- pic_matrix[f_sel, , drop = FALSE]
    n_iters <- ceiling(n_p_2 / n_peaks_per_round)
    p_cumu <- vector("list", n_iters)
    for (jj in seq_len(n_iters)) {
      peak_start <- (jj - 1) * n_peaks_per_round + 1
      peak_end <- min(n_p_2, jj * n_peaks_per_round)
      pic_dense <- as.matrix(pic_matrix_cumu[peak_start:peak_end, , drop = FALSE])
      log_msg(label, ": pacs_test_cumu block ", jj, "/", n_iters, " dim=", paste(dim(pic_dense), collapse = " x "))
      p_cumu[[jj]] <- fn$pacs_test_cumu(
        covariate_meta.data = covariate_meta.data,
        max_T = 2,
        formula_full = formula_full,
        formula_null = formula_null,
        pic_matrix = pic_dense,
        cap_rates = cap_rates,
        n_cores = n_cores,
        par_initial_null = par_initial_null,
        par_initial_full = par_initial_full
      )
    }
    rm(pic_matrix_cumu, pic_dense)
    gc(verbose = FALSE)
  }

  if (n_p_b >= 1) {
    pic_matrixbin_logit <- pic_matrixbin[f_b_sel, , drop = FALSE]
    n_iters_b <- ceiling(n_p_b / n_peaks_per_round)
    p_logit <- vector("list", n_iters_b)
    for (jj in seq_len(n_iters_b)) {
      peak_start <- (jj - 1) * n_peaks_per_round + 1
      peak_end <- min(n_p_b, jj * n_peaks_per_round)
      pic_dense <- as.matrix(pic_matrixbin_logit[peak_start:peak_end, , drop = FALSE])
      log_msg(label, ": pacs_test_logit block ", jj, "/", n_iters_b, " dim=", paste(dim(pic_dense), collapse = " x "))
      p_logit[[jj]] <- fn$pacs_test_logit(
        covariate_meta.data = covariate_meta.data,
        formula_full = formula_full,
        formula_null = formula_null,
        pic_matrix = pic_dense,
        cap_rates = cap_rates,
        n_cores = n_cores,
        par_initial_null = par_initial_null,
        par_initial_full = par_initial_full
      )
    }
    rm(pic_matrixbin_logit, pic_dense)
    gc(verbose = FALSE)
  }

  p_val_cumu <- if (length(p_cumu) > 0) {
    unlist(lapply(p_cumu, function(x) x$pacs_p_val), use.names = TRUE)
  } else {
    numeric(0)
  }
  p_val_logit <- if (length(p_logit) > 0) {
    unlist(lapply(p_logit, function(x) x$pacs_p_val), use.names = TRUE)
  } else {
    numeric(0)
  }
  p_val <- c(p_val_cumu, p_val_logit)[p_names]

  if (length(p_val) != nrow(pic_matrix)) {
    log_msg("ERROR: ", label, " p-value length ", length(p_val), " != nrow(pic_matrix) ", nrow(pic_matrix))
    stop(label, ": merged p-value length does not match nrow(pic_matrix)")
  }
  if (!identical(names(p_val), p_names)) {
    log_msg("ERROR: ", label, " merged p-value names do not match input peak names")
    stop(label, ": merged p-value names do not match rownames(pic_matrix)")
  }
  if (anyNA(p_val)) {
    missing_names <- p_names[is.na(p_val)]
    log_msg("ERROR: ", label, " merged p-values contain NA for first missing peaks: ", paste(head(missing_names, 10), collapse = ", "))
    stop(label, ": merged p-values contain NA")
  }

  convergence <- matrix(
    NA,
    nrow = nrow(pic_matrix),
    ncol = 2,
    dimnames = list(p_names, c("null", "full"))
  )

  list(pacs_converged = convergence, pacs_p_val = p_val)
}

log_msg("Starting PACS kidney Notebook 1 reproduction")
log_msg("Output directory: ", normalizePath(run_dir, mustWork = FALSE))
log_msg("Data directory: ", params$data_dir)

if (!dir.exists(params$data_dir)) {
  stop("Data directory does not exist: ", params$data_dir)
}

kidney_features_to_keep <- readRDS(file.path(params$data_dir, "kidney_features_to_keep.rds"))
load(file.path(params$data_dir, "data_for_test_for_t1e_power.rdata"))
r_by_ct_est <- readRDS(file.path(params$data_dir, "r_by_ct_est_kidney_adult.rds"))

stopifnot(exists("pmats"), exists("x.sp_cluster2"))
stopifnot("q_vec_new" %in% names(r_by_ct_est))

log_msg("Original pmats dim (cells x peaks): ", paste(dim(pmats), collapse = " x "))

# PT-only data for permuted labels / Type I error.
pmatpt <- pmats[x.sp_cluster2 == "PT", kidney_features_to_keep]
q_val <- r_by_ct_est$q_vec_new[x.sp_cluster2 == "PT"]

# PT+LOH data for actual labels / power.
actual_pmat <- pmats[x.sp_cluster2 %in% c("PT", "LOH"), kidney_features_to_keep]
actual_q_val <- r_by_ct_est$q_vec_new[x.sp_cluster2 %in% c("PT", "LOH")]
actual_cell_type_labels <- x.sp_cluster2[x.sp_cluster2 %in% c("PT", "LOH")]

# Notebook orientation: peaks x cells for PACS input.
p_by_c <- t(pmatpt)
actual_p_by_c <- t(actual_pmat)

# Remove outlier counts and drop explicit zeros, matching the notebook.
p_by_c@x[p_by_c@x >= params$count_outlier_cap] <- 0
actual_p_by_c@x[actual_p_by_c@x >= params$count_outlier_cap] <- 0
p_by_c <- drop0(p_by_c)
actual_p_by_c <- drop0(actual_p_by_c)

# Filter peaks separately for the two notebook scenarios.
p_by_c <- p_by_c[rowSums(p_by_c) > params$feature_total_count_cut, , drop = FALSE]
actual_p_by_c <- actual_p_by_c[rowSums(actual_p_by_c) > params$feature_total_count_cut, , drop = FALSE]

p_by_c1 <- actual_p_by_c[, actual_cell_type_labels == "PT", drop = FALSE]
q_val1 <- actual_q_val[actual_cell_type_labels == "PT"]
p_by_c2 <- actual_p_by_c[, actual_cell_type_labels == "LOH", drop = FALSE]
q_val2 <- actual_q_val[actual_cell_type_labels == "LOH"]

log_msg("PT-only p_by_c dim after filtering: ", paste(dim(p_by_c), collapse = " x "))
log_msg("Actual PT p_by_c1 dim after filtering: ", paste(dim(p_by_c1), collapse = " x "))
log_msg("Actual LOH p_by_c2 dim after filtering: ", paste(dim(p_by_c2), collapse = " x "))

n_repeat <- params$n_repeat
n_cell_sample <- params$n_cell_sample
n_features_sample <- params$n_features_sample

check_sample_size(nrow(p_by_c), n_features_sample, "n_features_sample against PT-only peaks")
check_sample_size(nrow(actual_p_by_c), n_features_sample, "n_features_sample against PT+LOH peaks")
check_sample_size(ncol(p_by_c), n_cell_sample * 2, "permutation PT cells")
check_sample_size(ncol(p_by_c1), n_cell_sample, "actual PT cells")
check_sample_size(ncol(p_by_c2), n_cell_sample, "actual LOH cells")

fn <- list(pacs_test_sparse = get_fun("pacs_test_sparse"))
fn$pacs_test_cumu <- get_fun("pacs_test_cumu")
fn$pacs_test_logit <- get_fun("pacs_test_logit")
if (is.null(fn$pacs_test_cumu) || is.null(fn$pacs_test_logit)) {
  stop("pacs_test_cumu and pacs_test_logit must be available from library(PACS)")
}

capture.output({
  cat("PACS package path:\n")
  print(system.file(package = "PACS"))
  cat("\nPACS package version:\n")
  print(packageVersion("PACS"))
  cat("\nPACS namespace exports:\n")
  print(getNamespaceExports("PACS"))
  for (nm in c("pacs_test_sparse", "pacs_test_cumu", "pacs_test_logit")) {
    cat("\n\n==== ", nm, " formals ====\n", sep = "")
    print(formals(fn[[nm]]))
    cat("\n==== ", nm, " source ====\n", sep = "")
    cat(paste(deparse(fn[[nm]]), collapse = "\n"))
    cat("\n")
  }
}, file = file.path(run_dir, "pacs_function_sources.txt"))

baseline_names <- c("seurat", "archR", "snapATAC", "fisher")
baseline_fun_names <- c(
  seurat = "seurat_method2_subsample",
  archR = "archR_method",
  snapATAC = "snapATAC_method",
  fisher = "fisher_method"
)

baseline_source <- list(
  type = "none",
  path = NA_character_,
  note = "run_baselines is FALSE; no baseline methods loaded"
)
available_baselines <- character()
if (params$run_baselines) {
  original_baseline_file <- file.path(script_dir, "other_methods_for_differential_updated.R")
  cleanroom_baseline_file <- file.path(script_dir, "baseline_methods_notebook1.R")
  if (file.exists(original_baseline_file)) {
    source(original_baseline_file)
    baseline_source <- list(
      type = "original_notebook",
      path = original_baseline_file,
      note = "Using original notebook baseline methods"
    )
    log_msg("Using original notebook baseline methods: ", original_baseline_file)
  } else if (file.exists(cleanroom_baseline_file)) {
    source(cleanroom_baseline_file)
    baseline_source <- list(
      type = "clean_room_reimplemented",
      path = cleanroom_baseline_file,
      note = "Using clean-room reimplemented baseline methods; not original author baseline"
    )
    log_msg("Using clean-room reimplemented baseline methods; not original author baseline: ", cleanroom_baseline_file)
  } else {
    stop(
      "run_baselines is TRUE, but no baseline method file was found in the project. ",
      "Expected other_methods_for_differential_updated.R or baseline_methods_notebook1.R"
    )
  }

  for (method in baseline_names) {
    f <- get_fun(baseline_fun_names[[method]])
    if (is.null(f)) {
      log_msg(
        "Baseline unavailable: ", method, " (", baseline_fun_names[[method]],
        "). This is expected with current library(PACS); original baseline source is missing."
      )
    } else {
      fn[[baseline_fun_names[[method]]]] <- f
      available_baselines <- c(available_baselines, method)
    }
  }
  missing_baselines <- setdiff(baseline_names, available_baselines)
  if (length(missing_baselines) > 0) {
    stop(
      "Baseline source was loaded, but these required methods are missing: ",
      paste(missing_baselines, collapse = ", ")
    )
  }
  log_msg(
    "Running available baseline methods: ", paste(available_baselines, collapse = ", "),
    ". Treat clean-room methods as approximate unless original notebook source is restored."
  )
} else {
  log_msg("run_baselines is FALSE; running PACS only.")
}

methods_all <- c("our", available_baselines)

p_value_permuted_label <- setNames(vector("list", length(methods_all)), methods_all)
p_value_actual_label <- setNames(vector("list", length(methods_all)), methods_all)
for (m in methods_all) {
  p_value_permuted_label[[m]] <- matrix(NA_real_, nrow = n_features_sample, ncol = n_repeat)
  p_value_actual_label[[m]] <- matrix(NA_real_, nrow = n_features_sample, ncol = n_repeat)
}

cells_sampled1_mat <- matrix(NA_integer_, nrow = n_cell_sample, ncol = n_repeat)
cells_sampled2_mat <- matrix(NA_integer_, nrow = n_cell_sample, ncol = n_repeat)
act_cells_sampled1_mat <- matrix(NA_integer_, nrow = n_cell_sample, ncol = n_repeat)
act_cells_sampled2_mat <- matrix(NA_integer_, nrow = n_cell_sample, ncol = n_repeat)
features_sampled_mat <- matrix(NA_integer_, nrow = n_features_sample, ncol = n_repeat)

for (iii in seq_len(n_repeat)) {
  act_cells_sampled1_mat[, iii] <- sample.int(ncol(p_by_c1), n_cell_sample, replace = FALSE)
  act_cells_sampled2_mat[, iii] <- sample.int(ncol(p_by_c2), n_cell_sample, replace = FALSE)

  all_cells <- sample.int(ncol(p_by_c), n_cell_sample * 2, replace = FALSE)
  cells_sampled1_mat[, iii] <- all_cells[seq_len(n_cell_sample)]
  cells_sampled2_mat[, iii] <- all_cells[(n_cell_sample + 1):(2 * n_cell_sample)]

  features_sampled_mat[, iii] <- sample.int(nrow(p_by_c), n_features_sample, replace = FALSE)
}

for (iii in seq_len(n_repeat)) {
  log_msg("Repeat ", iii, " / ", n_repeat)

  data_matrix_pos <- p_by_c[, cells_sampled1_mat[, iii], drop = FALSE]
  data_matrix_neg <- p_by_c[, cells_sampled2_mat[, iii], drop = FALSE]
  true_q_pos <- q_val[cells_sampled1_mat[, iii]]
  true_q_neg <- q_val[cells_sampled2_mat[, iii]]

  act_data_matrix_pos <- p_by_c1[, act_cells_sampled1_mat[, iii], drop = FALSE]
  act_data_matrix_neg <- p_by_c2[, act_cells_sampled2_mat[, iii], drop = FALSE]
  act_true_q_pos <- q_val1[act_cells_sampled1_mat[, iii]]
  act_true_q_neg <- q_val2[act_cells_sampled2_mat[, iii]]

  # Notebook behavior: compute read depths before feature subsetting.
  n_reads_cell <- c(colSums(data_matrix_pos), colSums(data_matrix_neg))
  act_n_reads_cell <- c(colSums(act_data_matrix_pos), colSums(act_data_matrix_neg))

  feature_idx <- features_sampled_mat[, iii]
  data_matrix_pos <- data_matrix_pos[feature_idx, , drop = FALSE]
  data_matrix_neg <- data_matrix_neg[feature_idx, , drop = FALSE]
  act_data_matrix_pos <- act_data_matrix_pos[feature_idx, , drop = FALSE]
  act_data_matrix_neg <- act_data_matrix_neg[feature_idx, , drop = FALSE]

  group.info <- data.frame(group = c(rep(0, ncol(data_matrix_pos)), rep(1, ncol(data_matrix_neg))))
  data_mat <- Matrix(cbind(data_matrix_pos, data_matrix_neg), sparse = TRUE)
  rownames(data_mat) <- paste0("f_", seq_len(nrow(data_mat)))
  colnames(data_mat) <- paste0("c_", seq_len(ncol(data_mat)))
  cap_rates <- c(true_q_pos, true_q_neg)
  check_pacs_input(data_mat, group.info, cap_rates, paste0("permuted repeat ", iii))

  act_group.info <- data.frame(group = c(rep(0, ncol(act_data_matrix_pos)), rep(1, ncol(act_data_matrix_neg))))
  act_data_mat <- Matrix(cbind(act_data_matrix_pos, act_data_matrix_neg), sparse = TRUE)
  rownames(act_data_mat) <- paste0("f_", seq_len(nrow(act_data_mat)))
  colnames(act_data_mat) <- paste0("c_", seq_len(ncol(act_data_mat)))
  act_cap_rates <- c(act_true_q_pos, act_true_q_neg)
  check_pacs_input(act_data_mat, act_group.info, act_cap_rates, paste0("actual repeat ", iii))

  our_p <- pacs_test_sparse_local_fixed(
    covariate_meta.data = group.info,
    formula_full = ~ factor(group),
    formula_null = ~ 1,
    pic_matrix = data_mat,
    n_peaks_per_round = NULL,
    T_proportion_cutoff = params$t_prop_cutoff,
    cap_rates = cap_rates,
    label = paste0("permuted repeat ", iii)
  )$pacs_p_val
  p_value_permuted_label[["our"]] <- store_pvalues(
    p_value_permuted_label[["our"]], our_p, n_features_sample, "our", "permuted", iii
  )

  act_our_p <- pacs_test_sparse_local_fixed(
    covariate_meta.data = act_group.info,
    formula_full = ~ factor(group),
    formula_null = ~ 1,
    pic_matrix = act_data_mat,
    n_peaks_per_round = NULL,
    T_proportion_cutoff = params$t_prop_cutoff,
    cap_rates = act_cap_rates,
    label = paste0("actual repeat ", iii)
  )$pacs_p_val
  p_value_actual_label[["our"]] <- store_pvalues(
    p_value_actual_label[["our"]], act_our_p, n_features_sample, "our", "actual", iii
  )

  if ("seurat" %in% available_baselines) {
    p_value_permuted_label[["seurat"]] <- store_pvalues(
      p_value_permuted_label[["seurat"]],
      fn$seurat_method2_subsample(data_matrix_pos, data_matrix_neg, n_reads_cell),
      n_features_sample, "seurat", "permuted", iii
    )
    p_value_actual_label[["seurat"]] <- store_pvalues(
      p_value_actual_label[["seurat"]],
      fn$seurat_method2_subsample(act_data_matrix_pos, act_data_matrix_neg, act_n_reads_cell),
      n_features_sample, "seurat", "actual", iii
    )
  }

  if (any(c("snapATAC", "fisher") %in% available_baselines)) {
    data_matrix_pos_bin <- data_matrix_pos
    data_matrix_pos_bin@x[data_matrix_pos_bin@x != 0] <- 1
    data_matrix_neg_bin <- data_matrix_neg
    data_matrix_neg_bin@x[data_matrix_neg_bin@x != 0] <- 1
    act_data_matrix_pos_bin <- act_data_matrix_pos
    act_data_matrix_pos_bin@x[act_data_matrix_pos_bin@x != 0] <- 1
    act_data_matrix_neg_bin <- act_data_matrix_neg
    act_data_matrix_neg_bin@x[act_data_matrix_neg_bin@x != 0] <- 1
  }

  if ("snapATAC" %in% available_baselines) {
    p_value_permuted_label[["snapATAC"]] <- store_pvalues(
      p_value_permuted_label[["snapATAC"]],
      fn$snapATAC_method(data_matrix_pos_bin, data_matrix_neg_bin, bcv = params$snap_bcv),
      n_features_sample, "snapATAC", "permuted", iii
    )
    p_value_actual_label[["snapATAC"]] <- store_pvalues(
      p_value_actual_label[["snapATAC"]],
      fn$snapATAC_method(act_data_matrix_pos_bin, act_data_matrix_neg_bin, bcv = params$snap_bcv),
      n_features_sample, "snapATAC", "actual", iii
    )
  }

  if ("fisher" %in% available_baselines) {
    p_value_permuted_label[["fisher"]] <- store_pvalues(
      p_value_permuted_label[["fisher"]],
      fn$fisher_method(data_matrix_pos_bin, data_matrix_neg_bin),
      n_features_sample, "fisher", "permuted", iii
    )
    p_value_actual_label[["fisher"]] <- store_pvalues(
      p_value_actual_label[["fisher"]],
      fn$fisher_method(act_data_matrix_pos_bin, act_data_matrix_neg_bin),
      n_features_sample, "fisher", "actual", iii
    )
  }

  if ("archR" %in% available_baselines) {
    p_value_permuted_label[["archR"]] <- store_pvalues(
      p_value_permuted_label[["archR"]],
      fn$archR_method(data_matrix_pos, data_matrix_neg),
      n_features_sample, "archR", "permuted", iii
    )
    p_value_actual_label[["archR"]] <- store_pvalues(
      p_value_actual_label[["archR"]],
      fn$archR_method(act_data_matrix_pos, act_data_matrix_neg),
      n_features_sample, "archR", "actual", iii
    )
  }
}

t1power_mat <- matrix(
  NA_real_,
  nrow = length(methods_all),
  ncol = 5,
  dimnames = list(methods_all, c("t1e", "t1e_sd", "power", "power_sd", "scenario"))
)
t1power_mat[, "scenario"] <- 1

for (jj in methods_all) {
  t1e_each_repeat <- colMeans(p_value_permuted_label[[jj]] < 0.05, na.rm = TRUE)
  t1power_mat[jj, "t1e"] <- mean(t1e_each_repeat, na.rm = TRUE)
  t1power_mat[jj, "t1e_sd"] <- sd(t1e_each_repeat, na.rm = TRUE)
}

cutoffs <- setNames(rep(NA_real_, length(methods_all)), methods_all)
for (jj in methods_all) {
  cutoffs[[jj]] <- min(0.05, quantile(p_value_permuted_label[[jj]], 0.05, na.rm = TRUE))
}

power_mat <- matrix(
  NA_real_,
  nrow = length(methods_all),
  ncol = n_repeat,
  dimnames = list(methods_all, paste0("rep", seq_len(n_repeat)))
)

union_methods <- intersect(c("our", "seurat", "archR", "snapATAC"), methods_all)
log_msg("Pseudo-true union methods: ", paste(union_methods, collapse = ", "))
if (identical(union_methods, "our")) {
  log_msg("PACS-only run: power for 'our' is defined against an 'our'-only union and is not the full notebook pseudo-true power comparison.")
}

for (ii in seq_len(n_repeat)) {
  union_true <- Reduce(
    "|",
    lapply(union_methods, function(m) p_value_actual_label[[m]][, ii] < cutoffs[[m]])
  )

  denom <- sum(union_true, na.rm = TRUE)
  if (denom == 0) {
    log_msg("Repeat ", ii, " has an empty pseudo-true union; power left as NA.")
    next
  }

  for (jj in methods_all) {
    power_mat[jj, ii] <- sum(
      p_value_actual_label[[jj]][, ii] < cutoffs[[jj]] & union_true,
      na.rm = TRUE
    ) / denom
  }
}

t1power_mat[, "power"] <- rowMeans(power_mat, na.rm = TRUE)
t1power_mat[, "power_sd"] <- apply(power_mat, 1, sd, na.rm = TRUE)

print(t1power_mat)

result <- list(
  params = params,
  run_dir = run_dir,
  dimensions = list(
    pmats = dim(pmats),
    p_by_c = dim(p_by_c),
    actual_p_by_c = dim(actual_p_by_c),
    p_by_c1 = dim(p_by_c1),
    p_by_c2 = dim(p_by_c2)
  ),
  cells_sampled1_mat = cells_sampled1_mat,
  cells_sampled2_mat = cells_sampled2_mat,
  act_cells_sampled1_mat = act_cells_sampled1_mat,
  act_cells_sampled2_mat = act_cells_sampled2_mat,
  features_sampled_mat = features_sampled_mat,
  p_value_permuted_label = p_value_permuted_label,
  p_value_actual_label = p_value_actual_label,
  cutoffs = cutoffs,
  union_methods = union_methods,
  power_mat = power_mat,
  t1power_mat = t1power_mat,
  baseline_source = baseline_source
)

saveRDS(result, file.path(run_dir, "pacs_kidney_notebook1_result.rds"))
write.csv(
  data.frame(
    method = rownames(t1power_mat),
    t1power_mat,
    baseline_source = baseline_source$type,
    row.names = NULL
  ),
  file.path(run_dir, "summary.csv"),
  row.names = FALSE
)

capture.output(sessionInfo(), file = file.path(run_dir, "session_info.txt"))

log_msg("Saved result RDS: ", file.path(run_dir, "pacs_kidney_notebook1_result.rds"))
log_msg("Saved summary CSV: ", file.path(run_dir, "summary.csv"))
log_msg("Saved session info: ", file.path(run_dir, "session_info.txt"))
log_msg("Finished")
