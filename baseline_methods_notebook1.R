# Clean-room baseline methods for PACS Notebook 1.
#
# This is a clean-room reimplementation for Notebook 1 baseline comparison.
# It is not the original other_methods_for_differential_updated.R.
# Results should be reported as reimplemented/approximate baselines.
#
# Design choices:
# - seurat_method2_subsample: per-peak logistic regression LRT with cell-group
#   label as outcome, log sequencing depth as an adjustment covariate, and raw
#   peak count as the tested predictor. This follows the Seurat LR baseline
#   description in spirit while keeping dependencies to base R + Matrix.
# - snapATAC_method: edgeR exact test on binarized peak accessibility counts,
#   with fixed dispersion bcv^2. This matches the Notebook 1 comment
#   "snapATAC-- edgeR method" more closely than a logistic-regression surrogate.
# - fisher_method: exact 2x2 Fisher test on binarized accessibility by group.
# - archR_method: approximate bias-aware baseline using nearest-neighbor depth
#   matching followed by Wilcoxon rank-sum tests on log1p depth-normalized
#   counts. This is not the official ArchR implementation, which includes
#   additional bias matching/model details.

suppressPackageStartupMessages({
  library(Matrix)
})

.as_numeric_row <- function(mat, i) {
  as.numeric(mat[i, ])
}

.safe_glm_lrt <- function(full_formula, null_formula, data) {
  full_fit <- tryCatch(
    suppressWarnings(glm(full_formula, data = data, family = binomial())),
    error = function(e) NULL
  )
  null_fit <- tryCatch(
    suppressWarnings(glm(null_formula, data = data, family = binomial())),
    error = function(e) NULL
  )
  if (is.null(full_fit) || is.null(null_fit)) {
    return(NA_real_)
  }
  dev_diff <- null_fit$deviance - full_fit$deviance
  df_diff <- null_fit$df.residual - full_fit$df.residual
  if (!is.finite(dev_diff) || !is.finite(df_diff) || df_diff <= 0) {
    return(NA_real_)
  }
  stats::pchisq(dev_diff, df = df_diff, lower.tail = FALSE)
}

seurat_method2_subsample <- function(mat_pos, mat_neg, n_reads_cell) {
  counts <- cbind(mat_pos, mat_neg)
  n_peaks <- nrow(counts)
  group <- c(rep(0L, ncol(mat_pos)), rep(1L, ncol(mat_neg)))

  if (length(n_reads_cell) != ncol(counts)) {
    stop("seurat_method2_subsample: length(n_reads_cell) must equal number of cells")
  }

  base_df <- data.frame(
    group = group,
    log_depth = log1p(as.numeric(n_reads_cell))
  )

  pvals <- vapply(seq_len(n_peaks), function(i) {
    peak_value <- .as_numeric_row(counts, i)
    if (length(unique(peak_value)) < 2L || length(unique(group)) < 2L) {
      return(NA_real_)
    }
    dat <- transform(base_df, peak_value = peak_value)
    .safe_glm_lrt(group ~ log_depth + peak_value, group ~ log_depth, dat)
  }, numeric(1))

  names(pvals) <- rownames(counts)
  pvals
}

snapATAC_method <- function(mat_pos_bin, mat_neg_bin, bcv = 0.4) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("snapATAC_method requires the edgeR package for the clean-room edgeR baseline")
  }

  counts <- as.matrix(cbind(mat_pos_bin, mat_neg_bin))
  group <- factor(
    c(rep("pos", ncol(mat_pos_bin)), rep("neg", ncol(mat_neg_bin))),
    levels = c("pos", "neg")
  )

  y <- edgeR::DGEList(counts = counts, group = group)
  y <- edgeR::calcNormFactors(y)
  test <- edgeR::exactTest(y, dispersion = bcv^2)
  pvals <- test$table$PValue
  names(pvals) <- rownames(counts)
  pvals
}

fisher_method <- function(mat_pos_bin, mat_neg_bin) {
  pos_access <- Matrix::rowSums(mat_pos_bin != 0)
  neg_access <- Matrix::rowSums(mat_neg_bin != 0)
  pos_inaccess <- ncol(mat_pos_bin) - pos_access
  neg_inaccess <- ncol(mat_neg_bin) - neg_access

  pvals <- vapply(seq_along(pos_access), function(i) {
    # Rows are accessibility status (accessible, inaccessible);
    # columns are groups (pos, neg). Orientation does not affect the two-sided
    # Fisher p-value, but keeping this layout makes the contrast explicit.
    tab <- matrix(
      c(pos_access[[i]], neg_access[[i]], pos_inaccess[[i]], neg_inaccess[[i]]),
      nrow = 2,
      byrow = FALSE
    )
    tryCatch(stats::fisher.test(tab)$p.value, error = function(e) NA_real_)
  }, numeric(1))

  names(pvals) <- rownames(mat_pos_bin)
  pvals
}

.depth_match_indices <- function(depth_pos, depth_neg) {
  if (length(depth_pos) <= length(depth_neg)) {
    target <- depth_pos
    pool <- depth_neg
    matched_neg <- vapply(target, function(d) which.min(abs(pool - d)), integer(1))
    list(pos = seq_along(depth_pos), neg = unique(matched_neg))
  } else {
    target <- depth_neg
    pool <- depth_pos
    matched_pos <- vapply(target, function(d) which.min(abs(pool - d)), integer(1))
    list(pos = unique(matched_pos), neg = seq_along(depth_neg))
  }
}

archR_method <- function(mat_pos, mat_neg) {
  counts <- cbind(mat_pos, mat_neg)
  n_pos <- ncol(mat_pos)
  depth_pos <- as.numeric(Matrix::colSums(mat_pos))
  depth_neg <- as.numeric(Matrix::colSums(mat_neg))
  matched <- .depth_match_indices(depth_pos, depth_neg)

  pos_cols <- matched$pos
  neg_cols <- n_pos + matched$neg
  selected_cols <- c(pos_cols, neg_cols)
  selected_counts <- counts[, selected_cols, drop = FALSE]
  selected_depth <- as.numeric(Matrix::colSums(selected_counts))
  selected_depth[selected_depth <= 0] <- 1
  scale_factor <- 10000 / selected_depth
  n_pos_matched <- length(pos_cols)
  n_peaks <- nrow(counts)

  pvals <- vapply(seq_len(n_peaks), function(i) {
    normalized <- log1p(.as_numeric_row(selected_counts, i) * scale_factor)
    x <- normalized[seq_len(n_pos_matched)]
    y <- normalized[(n_pos_matched + 1L):length(normalized)]
    if (length(unique(normalized)) < 2L) {
      return(NA_real_)
    }
    tryCatch(
      suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value),
      error = function(e) NA_real_
    )
  }, numeric(1))

  names(pvals) <- rownames(counts)
  pvals
}
