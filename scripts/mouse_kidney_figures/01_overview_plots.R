#!/usr/bin/env Rscript

# Mouse kidney figure reproduction: overview plots from existing PACS_data and
# Notebook 1 benchmark outputs. This script is read-only with respect to source
# data and writes figures under figures/mouse_kidney by default.

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
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
    defaults[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  defaults
}

params <- parse_args(list(
  data_dir = "/home/woodson/biostatistic/pacs/PACS_data",
  notebook1_result_dir = "/home/woodson/PACS_reproducing/results/20260526_2318_large_baseline",
  output_dir = "/home/woodson/PACS_reproducing/figures/mouse_kidney"
))

required_files <- c(
  file.path(params$data_dir, "data_for_test_for_t1e_power.rdata"),
  file.path(params$data_dir, "kidney_features_to_keep.rds"),
  file.path(params$data_dir, "r_by_ct_est_kidney_adult.rds"),
  file.path(params$notebook1_result_dir, "pacs_kidney_notebook1_result.rds"),
  file.path(params$notebook1_result_dir, "summary.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(params$data_dir, "data_for_test_for_t1e_power.rdata"))
kidney_features_to_keep <- readRDS(file.path(params$data_dir, "kidney_features_to_keep.rds"))
r_by_ct_est <- readRDS(file.path(params$data_dir, "r_by_ct_est_kidney_adult.rds"))
result <- readRDS(file.path(params$notebook1_result_dir, "pacs_kidney_notebook1_result.rds"))
summary_df <- read.csv(file.path(params$notebook1_result_dir, "summary.csv"), check.names = FALSE)

if (!exists("pmats") || !exists("x.sp_cluster2")) {
  stop("data_for_test_for_t1e_power.rdata must contain pmats and x.sp_cluster2")
}
if (!inherits(pmats, "Matrix")) {
  stop("pmats is expected to be a Matrix object")
}
if (!"q_vec_new" %in% names(r_by_ct_est)) {
  stop("r_by_ct_est_kidney_adult.rds must contain q_vec_new")
}
if (!"our" %in% names(result$p_value_permuted_label)) {
  stop("Notebook 1 result must contain p_value_permuted_label[['our']]")
}

theme_set(theme_bw(base_size = 12))
pal <- c(
  our = "#2f6f9f",
  seurat = "#c85a54",
  archR = "#5b8f5a",
  snapATAC = "#8b6bb1",
  fisher = "#d08c2d"
)

# 1. Cell type composition.
cell_counts <- as.data.frame(table(cell_type = x.sp_cluster2), stringsAsFactors = FALSE)
cell_counts <- cell_counts[order(cell_counts$Freq, decreasing = TRUE), ]
cell_counts$cell_type <- factor(cell_counts$cell_type, levels = cell_counts$cell_type)

p_cell_counts <- ggplot(cell_counts, aes(x = cell_type, y = Freq)) +
  geom_col(fill = "#4f7cac", width = 0.72) +
  coord_flip() +
  labs(x = NULL, y = "Cells", title = "Kidney Cell Type Composition") +
  theme(panel.grid.major.y = element_blank())

ggsave(
  file.path(params$output_dir, "cell_type_counts_barplot.png"),
  p_cell_counts,
  width = 7,
  height = max(4, 0.28 * nrow(cell_counts)),
  dpi = 300
)

# 2. PT vs LOH depth and capture-rate distribution.
pt_loh_cells <- x.sp_cluster2 %in% c("PT", "LOH")
if (!any(pt_loh_cells)) {
  stop("No PT or LOH cells found in x.sp_cluster2")
}
depth_df <- data.frame(
  cell_type = factor(x.sp_cluster2[pt_loh_cells], levels = c("PT", "LOH")),
  total_counts = as.numeric(Matrix::rowSums(pmats[pt_loh_cells, kidney_features_to_keep, drop = FALSE])),
  q_est = as.numeric(r_by_ct_est$q_vec_new[pt_loh_cells])
)

p_depth <- ggplot(depth_df, aes(x = total_counts, fill = cell_type, color = cell_type)) +
  geom_density(alpha = 0.24, linewidth = 0.6) +
  scale_x_log10() +
  scale_fill_manual(values = c(PT = "#2f6f9f", LOH = "#c85a54")) +
  scale_color_manual(values = c(PT = "#2f6f9f", LOH = "#c85a54")) +
  labs(x = "Total accessibility counts across kept kidney peaks (log10)", y = "Density", title = "PT vs LOH Depth Distribution") +
  theme(legend.title = element_blank())

ggsave(
  file.path(params$output_dir, "pt_loh_depth_distribution.png"),
  p_depth,
  width = 7,
  height = 4.5,
  dpi = 300
)

# 3. Notebook 1 benchmark Type I error and power.
summary_long <- rbind(
  data.frame(method = summary_df$method, metric = "Type I error", value = summary_df$t1e),
  data.frame(method = summary_df$method, metric = "Power", value = summary_df$power)
)
summary_long$method <- factor(summary_long$method, levels = summary_df$method)
summary_long$metric <- factor(summary_long$metric, levels = c("Type I error", "Power"))

p_benchmark <- ggplot(summary_long, aes(x = method, y = value, fill = method)) +
  geom_col(width = 0.72) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = pal[levels(summary_long$method)], na.value = "#777777") +
  labs(x = NULL, y = "Value", title = "Notebook 1 PACS Benchmark") +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

ggsave(
  file.path(params$output_dir, "pacs_benchmark_t1e_power_barplot.png"),
  p_benchmark,
  width = 8,
  height = 4.5,
  dpi = 300
)

# 4. PACS permuted-label QQ plot.
our_p <- as.numeric(result$p_value_permuted_label[["our"]])
our_p <- our_p[is.finite(our_p) & !is.na(our_p) & our_p > 0 & our_p <= 1]
if (length(our_p) < 10) {
  stop("Need at least 10 valid PACS permuted p-values for QQ plot")
}
our_p <- sort(our_p)
qq_df <- data.frame(
  expected = -log10(ppoints(length(our_p))),
  observed = -log10(our_p)
)
max_axis <- max(qq_df$expected, qq_df$observed, na.rm = TRUE)

p_qq <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_abline(slope = 1, intercept = 0, color = "grey55", linetype = "dashed") +
  geom_point(color = "#2f6f9f", alpha = 0.55, size = 0.65) +
  coord_equal(xlim = c(0, max_axis), ylim = c(0, max_axis)) +
  labs(x = "Expected -log10(p)", y = "Observed -log10(p)", title = "PACS Permuted-Label QQ Plot")

ggsave(
  file.path(params$output_dir, "pacs_permuted_qq_plot.png"),
  p_qq,
  width = 5,
  height = 5,
  dpi = 300
)

cat("Saved overview figures to: ", normalizePath(params$output_dir, mustWork = FALSE), "\n", sep = "")
