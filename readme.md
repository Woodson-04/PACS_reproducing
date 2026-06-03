# PACS 论文复现项目

## 1. 项目目标

本项目目标是复现 PACS 论文中的关键分析流程，并逐步形成可展示、可追踪、可继续扩展的复现项目。

目前项目分为两个主要阶段：

1. **Notebook 1 benchmark 复现**：复现作者 Notebook 1 中基于 real kidney data 的 Type I error 和 power workflow。
2. **mouse kidney / GSE157079 图复现**：从 GSE157079 mouse kidney snATAC 数据出发，转向 PACS 论文中的 batch-effect feature removal UMAP 逻辑。

当前项目主线是 **PACS paper-style reconstruction**，也就是围绕 PACS 论文中“检测 batch-effect features、移除这些 features、重新构建 UMAP”的逻辑推进。项目不是完整复现原始 `dev-kidney-snATAC` atlas；该仓库只作为理解原始 mouse kidney snATAC 数据处理背景的参考。

## 2. 项目路径和数据约束

当前项目目录：

```text
/home/woodson/PACS_reproducing
```

旧项目、参考代码和原始数据目录：

```text
/home/woodson/biostatistic/pacs
```

重要约束：

- `/home/woodson/biostatistic` 只读，不应修改、移动、删除或覆盖其中任何文件。
- 不修改已安装的 PACS R 包。
- 不解压、覆盖或编辑原始 GSE157079 `.gz` 文件。
- 大型 GSE157079 matrix 文件只通过流式方式读取，不解压到磁盘。
- 当前项目中的新代码、报告、图和中间结果均保存在 `/home/woodson/PACS_reproducing` 内。

GSE157079 关键数据文件：

```text
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_UMAP_coordinates.csv.gz
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_metadata.csv.gz
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_peak_list.csv.gz
/home/woodson/biostatistic/pacs/GSE157079/GSE157079_snATAC_cell_by_peak_matrix.txt.gz
```

## 3. AI 辅助与个人工作说明

本项目采用人机协作方式推进，较多使用了 ChatGPT 和 Codex。为了保证学术沟通中的透明性，这里明确说明各方角色。

ChatGPT 和 Codex 在项目中承担了大量辅助工作，包括复现路线规划、代码设计、R 脚本撰写、错误分析、调试建议、文件整理、报告初稿生成、结果解释框架整理和后续任务拆解。Codex 主要参与代码实现、脚本修改、项目结构整理、报告文件生成和命令执行建议；ChatGPT 主要参与研究路线讨论、统计和生物信息学解释、结果总结和汇报材料组织。

本人主要负责研究方向选择、复现目标确定、科学主线判断、数据路径和运行环境提供、关键 Linux/R 命令的实际运行、输出结果检查、是否调整路线的决策，以及面向老师的汇报和科学解释。也就是说，AI 工具承担了较多实现和组织层面的工作，但研究问题选择、结果核验、路线取舍和最终解释由本人负责。

因此，本项目应被理解为一个 **human-AI collaborative reproduction workflow**。后续汇报中会明确说明 AI 辅助的范围，既不夸大个人手工编码贡献，也不把研究工作表述为 AI 独立完成。

## 4. 第一阶段：Notebook 1 benchmark 复现

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

结果摘要：

| method | Type I error | power | 说明 |
|---|---:|---:|---|
| our / PACS | 0.04008 | 0.83337 | PACS 主方法，接近作者 Notebook |
| seurat | 0.06342 | 0.82344 | clean-room baseline |
| archR | 0.04096 | 0.67437 | clean-room approximate baseline |
| snapATAC | 0.01810 | 0.76094 | clean-room edgeR-style baseline |
| fisher | 0.02208 | 0.76630 | binary Fisher exact test |

解释：

- PACS 主方法结果与作者 Notebook 1 数值接近，可以说明 PACS benchmark workflow 已基本复现。
- baseline 方法不能声称完全复现作者原始 baseline，因为作者 helper 文件 `other_methods_for_differential_updated.R` 未找到。
- 当前 baseline 是 clean-room reimplementation，应在汇报中标注为近似/重实现版本。
- PACS 0.2.2 中 `pacs_test_sparse()` 存在 mixed cumulative/logit branch 的 `rownames` 合并 bug。`q.r` 中实现了本地 fixed wrapper，直接调用 `pacs_test_cumu()` 和 `pacs_test_logit()` 后安全合并 p value。

## 5. 第二阶段：GSE157079 mouse kidney 图复现

GSE157079 metadata、UMAP、peak list 和 matrix 文件已经接入项目。

标准化后的 metadata 输出：

```text
results/mouse_kidney_figures/gse157079_metadata_merged.csv
```

标准列：

```text
row_index
cell_barcode
sample
cell_type
umap_1
umap_2
```

重要修复：

- 修复了 `data.table` 下用 `df[row, col]` 导致列索引解释错误的问题。
- 原始 CSV 第一列为空列名，已标准化为 `row_index`。
- 10x-style `cell_barcode` 在不同 sample 间可能重复，因此不把 barcode 作为全局唯一键；metadata 和 UMAP 用 `row_index` 合并，`sample + cell_barcode` 只作为 sanity check。

matrix alignment smoke test 结果：

```text
matrix = 28316 cells x 300755 peaks
nonzero entries = 166121193
metadata rows = 28316
peak list rows = 300755
orientation = cell x peak
```

报告：

```text
results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md
```

GEO precomputed UMAP 已绘制为 overview 图，但它不等同于 PACS 论文中 all-feature vs PACS-filtered UMAP。后续主线使用 matrix-derived embedding，而不是直接使用 GEO 的 `umap_1` / `umap_2`。

## 6. Matrix-derived LSI/UMAP 当前进展

pilot 脚本：

```text
scripts/mouse_kidney_figures/05_pilot_gse157079_lsi_umap_from_matrix.R
```

该脚本使用抽样 cells 和随机 peaks，从 MatrixMarket 文件流式构建 sparse matrix，并完成 TF-IDF / LSI / UMAP。

all-cell top-peak 脚本：

```text
scripts/mouse_kidney_figures/06_all_cells_top_peaks_lsi_umap_from_matrix.R
```

当前最佳结果：

```text
all cells = 28316
top detected peaks = 20000
retained nonzeros = 85801336
sparse matrix = 28316 x 20000
LSI = 28316 x 50
UMAP uses LSI dims 2:30
```

输出图：

```text
figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png
figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png
```

解释：

- sample-colored UMAP 显示明显 sample-associated structure，说明 batch/sample effect 仍然存在。
- cell-type-colored UMAP 仍保留一定生物学结构。
- 这组图可以作为当前 **before-filtering baseline**，为下一步 PACS batch-effect feature removal 做参照。

## 7. 当前可展示成果

建议会议展示材料：

1. `results/project_progress_summary_for_teacher.md`
2. `results/20260526_2318_large_baseline/summary.csv`
3. `figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png`
4. `figures/mouse_kidney/pacs_permuted_qq_plot.png`
5. `results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md`
6. `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png`
7. `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png`
8. `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md`

## 8. 下一步计划

下一步技术路线：

1. 使用 PACS 在 GSE157079 中检测 batch-associated peaks。
2. 移除显著 batch-effect peaks。
3. 用过滤后的 feature set 重新计算 TF-IDF / LSI / UMAP。
4. 构建 before/after 四联图：
   - before by sample；
   - after by sample；
   - before by cell type；
   - after by cell type。
5. 比较移除 batch-effect peaks 后 sample-associated structure 是否减弱，同时 cell-type biological structure 是否保留。

## 9. 目录说明

```text
scripts/
```

项目脚本目录。

```text
scripts/mouse_kidney_figures/
```

GSE157079 intake、matrix alignment、matrix-derived UMAP 和后续 mouse kidney figure reconstruction 脚本。

```text
results/
```

Notebook 1 benchmark 结果、GSE157079 中间结果、报告和会议用总结。

```text
figures/mouse_kidney/
```

mouse kidney 相关图，包括 GEO UMAP、matrix-derived pilot UMAP 和 all-cell top-peak UMAP。

```text
meeting_materials_20260530/
```

会议用材料目录，包含精简报告、图和表格。

```text
archive_pre_meeting_cleanup_20260530/
```

会前整理归档目录。当前环境下未强制删除文件；归档操作应在确认后执行。
