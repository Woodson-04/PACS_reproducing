# PACS 下一阶段交接说明

本文档用于开启新的 ChatGPT 对话或交给 Codex 继续工作。它不是对全部聊天记录的流水账总结，而是面向下一阶段复现工作的结构化背景说明。

当前仓库：

```text
/home/woodson/PACS_reproducing
```

参考旧项目与原始数据目录：

```text
/home/woodson/biostatistic/pacs
```

重要原则：

- `/home/woodson/biostatistic` 只读，不修改、不移动、不删除。
- 不修改已安装 PACS R 包。
- 不解压、覆盖或编辑原始 GSE157079 `.gz` 文件。
- 大型 MatrixMarket 文件只通过流式方式读取。
- 新代码、结果、图、报告只写入 `/home/woodson/PACS_reproducing`。

---

## 1. 项目总目标

本项目目标是复现 PACS 论文中的关键分析流程，并形成可展示、可追踪、可继续扩展的复现项目。

当前已经完成从 **Notebook 1 benchmark 复现** 到 **GSE157079 mouse kidney P56 batch-effect feature filtering** 的阶段性推进。当前主线是：

```text
PACS paper-style reconstruction
= 检测 batch-effect features
→ 移除 batch-associated features
→ 重新构建 UMAP
→ 定量比较 batch mixing 与 cell-type structure preservation
```

当前结果应表述为：

```text
P56 two-batch top10000 setting 下的 Fig.3a–d author-style reconstruction
```

不能表述为：

```text
exact full Fig.3 reproduction
```

---

## 2. 当前已经完成的工作

### 2.1 Notebook 1 benchmark

主脚本：

```text
q.r
```

large benchmark run 参数：

```text
n_repeat = 5
n_cell_sample = 500
n_features_sample = 10000
run_baselines = TRUE
output_dir = results/20260526_2318_large_baseline
```

核心结果：

| method | Type I error | power | 说明 |
|---|---:|---:|---|
| our / PACS | 0.04008 | 0.83337 | PACS 主方法 |
| seurat | 0.06342 | 0.82344 | clean-room baseline |
| archR | 0.04096 | 0.67437 | clean-room approximate baseline |
| snapATAC | 0.01810 | 0.76094 | clean-room edgeR-style baseline |
| fisher | 0.02208 | 0.76630 | binary Fisher exact test |

说明：

- PACS 主方法结果与作者 Notebook 1 数值接近。
- baseline 方法不能声称完全复现作者原始 baseline，因为作者 helper 文件 `other_methods_for_differential_updated.R` 未找到。
- 当前 baseline 是 clean-room reimplementation。
- PACS 0.2.2 中 `pacs_test_sparse()` 存在 mixed cumulative/logit branch 的 `rownames` 合并 bug，`q.r` 中实现了 local fixed wrapper。

---

### 2.2 GSE157079 数据接入与矩阵对齐

GSE157079 原始数据文件位于：

```text
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_UMAP_coordinates.csv.gz
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_metadata.csv.gz
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_peak_list.csv.gz
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz
```

已完成 metadata、UMAP、peak list 和 MatrixMarket sparse matrix 对齐。

核心结果：

| object | value |
|---|---:|
| cells | 28,316 |
| peaks | 300,755 |
| nonzero entries | 166,121,193 |
| metadata rows | 28,316 |
| peak list rows | 300,755 |
| matrix orientation | cell × peak |

标准 metadata 输出：

```text
results/mouse_kidney_figures/gse157079_metadata_merged.csv
```

重要修复：

- 原始 CSV 第一列空列名已标准化为 `row_index`。
- metadata 与 GEO UMAP coordinates 使用 `row_index` 合并。
- 10x-style `cell_barcode` 不能作为全局唯一键，因为 barcode 可跨 sample 重复。

---

### 2.3 All-cell top-peak matrix-derived UMAP

已从 GSE157079 原始 cell-by-peak MatrixMarket 矩阵出发，完成 all-cell top-peak TF-IDF / LSI / UMAP。

主脚本：

```text
scripts/mouse_kidney_figures/06_all_cells_top_peaks_lsi_umap_from_matrix.R
```

当前结果：

| item | value |
|---|---:|
| cells | 28,316 |
| selected top peaks | 20,000 |
| retained nonzeros | 85,801,336 |
| sparse matrix | 28,316 × 20,000 |
| LSI dimension | 28,316 × 50 |
| UMAP input | LSI_2:LSI_30 |

核心图：

```text
figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png
figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png
```

注意：

- 这不是 all-features UMAP。
- 这不是 PACS-filtered UMAP。
- 它只是 all-cell top-peak matrix-derived reference embedding。
- 完整 300,755 peaks 的 all-features UMAP 已尝试运行，但在 sparse Matrix 构建阶段被系统终止，推测主要原因是内存压力。

---

### 2.4 P56 PACS batch-effect peak filtering

当前主设置：

```text
n_top_peaks = 10000
max_pacs_peaks = 10000
fdr_cutoff = 0.05
```

结果目录：

```text
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved
```

P56 子集：

| item | value |
|---|---:|
| P56 cells | 13,526 |
| P56_batch1 | 7,129 |
| P56_batch2 | 6,397 |
| cell types | 15 |

PACS model：

```text
full model: ~ cell_type + batch
null model: ~ cell_type
```

PACS peak results：

| item | value |
|---|---:|
| tested peaks | 10,000 |
| significant batch peaks | 6,305 |
| retained peaks | 3,695 |
| before matrix | 13,526 × 10,000 |
| after matrix | 13,526 × 3,695 |

核心图：

```text
figures/mouse_kidney/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png
```

解释：

- before by batch 显示 P56_batch1 与 P56_batch2 有明显 batch separation。
- after by batch 显示移除 PACS-significant batch peaks 后 batch separation 减弱。
- before/after by cell type 用于评估 biological cell-type structure 是否保留。

---

### 2.5 Batch mixing quantification

本项目已经在 PCA-logNorm、LSI 和 UMAP 三类坐标空间中计算 paper-style normalized batch mixing score。

主结果：

| coordinate space | before | after | role |
|---|---:|---:|---|
| PCA-logNorm PC1:30 | 0.0510 | 0.3097 | 主指标；最接近原文 normalized PCA mixing |
| PCA-logNorm PC1:50 | 0.0685 | 0.4239 | PCA 维度敏感性 |
| LSI_2:30 | 0.0520 | 0.3343 | scATAC-aware supporting metric |
| UMAP-space | 0.02649 | 0.65073 | visualization-level auxiliary metric |

核心图：

```text
figures/mouse_kidney/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png
```

结论：

- PACS filtering 后三类空间中的 batch mixing score 均提高。
- 因此当前结果不能只解释为 UMAP visualization artifact。

---

### 2.6 LSI_1-depth QC

LSI_1 与 cell depth 高度相关，因此 LSI-space metric 使用 LSI_2:30。

| stage | component | Spearman(depth) | Pearson(depth) |
|---|---|---:|---:|
| before | LSI_1 | -0.9983 | -0.9739 |
| after | LSI_1 | 0.9976 | 0.9830 |

核心图：

```text
figures/mouse_kidney/gse157079_p56_lsi_depth_correlation_before_after.png
```

结论：

```text
LSI_1 is essentially a depth/QC axis.
Excluding LSI_1 and using LSI_2:30 is justified.
```

---

### 2.7 PCA-logNorm parameter stability

已重新计算 PCA-logNorm space 中的 parameter stability，不再只依赖 UMAP-space auxiliary metrics。

结果文件：

```text
results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space_parameter_stability/p56_pca_parameter_stability_scores.csv
```

核心结果：

| setting | PC1:30 before | PC1:30 after | PC1:50 before | PC1:50 after | retained peaks |
|---|---:|---:|---:|---:|---:|
| 5000 tested | 0.0510 | 0.3679 | 0.0685 | 0.5032 | 1,860 |
| 10000 tested | 0.0510 | 0.3097 | 0.0685 | 0.4239 | 3,695 |
| 20000 top / 10000 tested | 0.0398 | 0.2928 | 0.0556 | 0.4097 | 3,792 |

结论：

- 三组参数均显示 after > before。
- 5000 tested after score 最高，但 retained peaks 最少，可能更激进。
- 当前仍推荐 10000 top / 10000 tested / FDR0.05 作为主展示设置，因为它在 batch mixing improvement、retained features 与 cell-type preservation 之间更平衡。

---

## 3. 当前仓库和报告状态

README 已更新为对外展示版本，包含最重要图表和结果摘要。

阶段性 PDF/HTML 报告已基本完成，核心定位是：

```text
PACS 复现项目阶段报告：GSE157079 UMAP 重建与 P56 批次效应特征去除
P56 two-batch top10000 setting 下的 Fig.3a–d author-style reconstruction
```

报告相关目录：

```text
report/pdf_report/
```

---

## 4. 必须避免的错误表述

后续写作、汇报和 README 中必须避免以下表述：

1. 不能说当前结果是 exact full Fig.3 reproduction。
2. 不能说 all-cell top-20,000 peak UMAP 是 all-features UMAP。
3. 不能说 GEO precomputed UMAP 是 PACS-filtered ground truth。
4. 不能说 Notebook 1 baseline 完全复现作者 baseline；只能说 clean-room reimplementation。
5. 不能忽略 `cap_rates` 是 depth-derived approximation，因为公开 GSE157079 文件没有作者 Notebook 1 的 `q_vec`。
6. 不能把 P56-only two-batch subset 说成完整 adult kidney atlas。
7. 不能把 UMAP-space score 当作主指标；PCA-logNorm 更接近 PACS 原文 normalized PCA mixing。
8. 不能把 LSI_1 直接纳入主 LSI metric；LSI_1 已证明与 depth 高度相关。

---

## 5. 下一阶段核心问题：scATAC-seq 测序深度 / capture-rate bias

下一阶段建议从“batch-effect feature filtering”转向“测序深度 / capture-rate bias 的建模与校正”。这是 PACS 的方法学核心之一，也最能体现数学和统计理解。

### 5.1 为什么这个方向值得做

scATAC-seq 中每个细胞的总 fragments / accessible counts 差异很大。高深度细胞更容易在更多 peaks 上观测到 accessibility，低深度细胞则容易出现 dropout 或 sparsity。这会导致两个问题：

1. 技术深度差异被误判为 biological accessibility difference。
2. 降维结果中的第一主轴或第一 LSI component 可能主要反映 depth，而非细胞类型或真实调控状态。

当前项目已经发现：

```text
before LSI_1 Spearman(depth) = -0.9983
after  LSI_1 Spearman(depth) =  0.9976
```

这说明 depth 是当前数据中非常强的技术轴。

---

### 5.2 下一阶段优先选择哪张图

优先建议：

```text
优先选择 PACS 论文中与 sequencing depth / capture rate correction / Type I error calibration 直接相关的 benchmark 或方法验证图。
```

不要立刻进入 Fig.3e/f downstream biology。理由是：

- Fig.3e/f 通常更偏 downstream biological interpretation，例如 cell-type-specific peaks、motif、gene activity 或轨迹展示；
- 这些图需要额外 annotation、motif database、gene linkage 或 scRNA reference，工作量大且容易偏离 PACS 的统计核心；
- 你目前最有竞争力的方向是用数学和统计解释 PACS 如何处理 scATAC 的 depth/capture bias。

更合适的下一阶段图应满足以下条件：

| 优先级 | 图类型 | 是否建议 | 原因 |
|---|---|---|---|
| 1 | depth/capture rate correction 相关 benchmark 图 | 强烈建议 | 直接体现 PACS 方法学核心 |
| 2 | p-value calibration / Type I error under depth confounding | 强烈建议 | 可延续 Notebook 1 benchmark，数学含量高 |
| 3 | LSI/PCA component-depth correlation diagnostic | 建议 | 可基于现有结果扩展，成本低 |
| 4 | Fig.3e/f downstream peak-to-gene / motif / heatmap | 暂缓 | 生物注释复杂，和 depth 主题不完全一致 |
| 5 | full all-features UMAP | 暂缓 | 主要受内存限制，方法学收益不如 depth 分析 |

如果论文中有专门展示 PACS 对测序深度、capture rate、cell-specific detection probability 或 p-value calibration 的图，应优先复现那张图。若不确定图号，下一阶段第一步应让 Codex 读取论文 PDF、supplement、作者 notebook 和 README，定位所有与以下关键词相关的图：

```text
sequencing depth
read depth
capture rate
q_i
cell-specific effect
cell-specific depth
calibration
Type I error
power
p-value
permuted label
coverage
library size
```

然后再决定具体复现哪一张。

---

### 5.3 如果论文没有直接的 depth figure，建议自建一张“方法学补充图”

即使论文没有单独的 depth figure，也可以基于当前结果构建一个新的、很有价值的 reproduction-extension figure：

```text
Depth-aware QC and calibration panel for GSE157079 P56
```

建议四联图：

1. cell depth distribution by batch；
2. LSI_1 vs log(depth) scatter；
3. depth-stratified batch mixing / cell-type mixing；
4. PACS p-value 或 significant peak fraction 与 peak depth / detection frequency 的关系。

这张图的价值：

- 直接回应 scATAC-seq 测序深度问题；
- 与 PACS 的 cell-specific correction 思想一致；
- 能体现数学统计分析能力；
- 不依赖复杂 motif/gene annotation；
- 可作为当前 Fig.3a–d reconstruction 的方法学延伸。

---

## 6. 给 Codex 的下一阶段任务建议

下一阶段不要一开始就大规模重跑所有 pipeline。先做 paper/source inspection，再决定图。

### 任务 A：定位 PACS 论文中 depth/capture-rate 相关图

让 Codex 执行：

```text
Search the repository and local reference directories for PACS paper PDF, supplementary PDF, notebooks, and scripts.
Identify all figures, panels, notebooks, or result sections related to sequencing depth, capture rate, q_i, cell-specific depth effect, p-value calibration, Type I error under null/permuted labels, and depth confounding.
Create a report listing candidate figures/panels for next-stage reproduction.
Do not run heavy analysis yet.
```

输出建议：

```text
results/next_stage_depth_bias/depth_related_figure_candidates.md
```

报告内容应包括：

- 候选图编号；
- 图的科学问题；
- 所需数据；
- 当前仓库是否已有数据；
- 复现难度；
- 推荐优先级。

---

### 任务 B：基于现有 P56 数据做 depth diagnostics

在不重跑 PACS 的前提下，优先生成：

```text
results/next_stage_depth_bias/
figures/mouse_kidney/
```

建议脚本：

```text
scripts/mouse_kidney_figures/10_p56_depth_bias_diagnostics.R
```

输入：

```text
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_counts_top_peaks_sparse.rds
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_metadata.csv
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_before_lsi_embedding.csv
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_after_lsi_embedding.csv
results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_pacs_batch_peak_results.csv
```

输出图建议：

```text
figures/mouse_kidney/gse157079_p56_depth_distribution_by_batch.png
figures/mouse_kidney/gse157079_p56_lsi1_vs_depth_scatter.png
figures/mouse_kidney/gse157079_p56_peak_detection_vs_pacs_significance.png
figures/mouse_kidney/gse157079_p56_depth_bias_diagnostic_four_panel.png
```

输出报告：

```text
results/next_stage_depth_bias/p56_depth_bias_diagnostics_report.md
```

---

### 任务 C：如果要进一步贴近 PACS 方法，重构 depth/capture-rate 模型解释

需要整理 PACS 中 `q_i`、`cap_rates` 或 cell-specific capture probability 的数学含义。

建议报告：

```text
results/next_stage_depth_bias/pacs_depth_correction_math_note.md
```

内容包括：

- scATAC counts 与 depth/capture rate 的关系；
- 为什么简单 logistic / Fisher / naive differential test 会受 depth 影响；
- PACS 如何用 cell-specific correction 控制 depth；
- 当前 GSE157079 public data 缺少作者 q_vec，因此本项目用 depth-derived approximation；
- 这个 approximation 对 p-value calibration 和 downstream filtering 的潜在影响。

---

## 7. 下一阶段推荐工作顺序

推荐顺序：

1. **冻结当前阶段**：README、PDF 报告、核心图表、结果路径全部确认。
2. **定位论文 depth 相关图**：先不要盲目做 Fig.3e/f，先找论文中是否有 depth/capture-rate correction 图。
3. **做 P56 depth diagnostics**：基于现有 saved matrix 和 LSI embedding，先生成 depth distribution、LSI_1-depth scatter、peak detection vs PACS p-value 等轻量结果。
4. **决定是否复现论文对应图**：如果论文有明确 depth panel，就复现；如果没有，就把 depth diagnostics 作为 reproduction-extension。
5. **再考虑 Fig.3e/f**：只有在 depth 方向完成后，再进入 downstream biological interpretation。

---

## 8. 新对话启动 prompt

可以把下面这段复制到新的 ChatGPT 对话中：

```text
我正在继续 PACS 论文复现项目。当前 GitHub 仓库是 Woodson-04/PACS_reproducing，本地路径是 /home/woodson/PACS_reproducing。旧项目和原始数据在 /home/woodson/biostatistic/pacs，只读。

当前已经完成：Notebook 1 benchmark；GSE157079 metadata / matrix / peak list 对齐；all-cell top-20000 peak TF-IDF/LSI/UMAP；P56 two-batch top10000 PACS batch-effect peak filtering；before/after UMAP 四联图；PCA-logNorm / LSI / UMAP 三层 normalized batch mixing score；LSI_1-depth QC；README 和阶段性 PDF 报告。

当前主结果：P56 cells = 13526；P56_batch1 = 7129；P56_batch2 = 6397；tested peaks = 10000；significant batch peaks = 6305；retained peaks = 3695；PCA-logNorm PC1:30 score 0.0510 -> 0.3097；PCA-logNorm PC1:50 score 0.0685 -> 0.4239；LSI_2:30 score 0.0520 -> 0.3343；UMAP-space score 0.02649 -> 0.65073。

注意：当前结果是 P56 two-batch top10000 setting 下的 Fig.3a-d author-style reconstruction，不是 exact full Fig.3 reproduction。all-cell top-20000 UMAP 不是 all-features UMAP。GEO UMAP 不是 PACS-filtered ground truth。baseline 是 clean-room reimplementation。cap_rates 是 depth-derived approximation。

下一阶段我想重点解决 scATAC-seq 的测序深度 / capture-rate bias 问题。请先帮我判断 PACS 论文中最适合复现哪一张与 sequencing depth / capture rate / q_i / Type I error calibration / p-value calibration 相关的图。如果论文没有直接 depth figure，则请设计一个基于当前 P56 数据的 depth diagnostics reproduction-extension panel。
```

---

## 9. 一句话总结

当前阶段已经完成 PACS 从 benchmark 到 P56 Fig.3a–d author-style UMAP reconstruction 的核心闭环；下一阶段最值得推进的方向不是马上做 downstream motif/gene 图，而是系统处理 scATAC-seq 测序深度 / capture-rate bias，因为这更贴近 PACS 的方法学核心，也更能体现数学统计能力。
