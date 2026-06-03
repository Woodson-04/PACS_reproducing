# 会前汇报指南

## 3 分钟版本

1. 项目目标是复现 PACS 论文，先完成 Notebook 1 benchmark，再转向 mouse kidney / GSE157079 图复现。
2. Notebook 1 benchmark 已经跑通。PACS 主方法 Type I error = 0.04008，power = 0.83337，和作者结果接近。
3. GSE157079 的 metadata、peak list 和大矩阵已经完成对齐检查：28316 cells x 300755 peaks，166121193 nonzero entries，矩阵方向是 cell x peak。
4. 已经从矩阵本身而不是 GEO 预计算 UMAP 出发，完成 all-cell top-20000 peak TF-IDF / LSI / UMAP。
5. 当前 UMAP 显示 sample-associated structure 明显，同时 cell type structure 仍存在。下一步就是用 PACS 找 batch-effect peaks，移除后重算 UMAP，做 before/after 图。

## 8 分钟版本

1. **项目定位**：这是 PACS 论文复现项目，不是完整复现原始 dev-kidney-snATAC atlas。
2. **Notebook 1 benchmark**：严格按作者 Notebook 1 的 PT-only permutation 和 PT+LOH actual label 逻辑复现；PACS 主方法结果接近作者。
3. **baseline 说明**：作者原始 baseline helper 文件缺失，所以 Seurat、ArchR、snapATAC、Fisher 是 clean-room reimplementation，不声称精确复现作者 baseline。
4. **PACS 包问题**：PACS 0.2.2 `pacs_test_sparse()` 有 mixed branch rownames bug，已在 `q.r` 中用本地 fixed wrapper 绕过。
5. **GSE157079 intake**：修复 metadata/UMAP 第一列、barcode 不全局唯一、data.table 索引等问题，得到标准 metadata。
6. **matrix alignment**：确认矩阵维度、metadata 行数、peak list 行数和 cell x peak 方向。
7. **matrix-derived UMAP**：先做 sampled pilot，再做 all-cell top-peak UMAP，不使用 GEO 预计算坐标。
8. **下一步**：PACS batch-effect peak detection -> remove peaks -> recompute UMAP -> before/after four-panel figure。

## 推荐展示顺序

1. `results/project_progress_summary_for_teacher.md`
2. `results/20260526_2318_large_baseline/summary.csv`
3. `figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png`
4. `figures/mouse_kidney/pacs_permuted_qq_plot.png`
5. `results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md`
6. `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png`
7. `figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png`
8. `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md`

## 可能问题与回答

### 1. PACS 和 baseline 差异为什么存在？

PACS 主方法是按作者 Notebook 1 的核心逻辑复现的，结果接近作者。baseline 差异主要来自作者原始 helper 文件 `other_methods_for_differential_updated.R` 未找到，因此当前 baseline 是 clean-room reimplementation。它们可用于流程比较，但不应声称完全复现作者 baseline 数字。

### 2. 为什么 GEO UMAP 和论文 UMAP 不一样？

GEO 文件中的 UMAP 是预计算坐标，只适合作为数据 overview。PACS 论文的关键逻辑是 all-feature UMAP 和移除 batch-effect features 后的 filtered-feature UMAP。后者需要从 cell-by-peak matrix 重新计算 embedding，而不是直接画 GEO 坐标。

### 3. 为什么现在 sample UMAP 分离明显？

当前 all-cell top-peak UMAP 是 before-filtering baseline，使用 top detected peaks，没有移除 batch-associated peaks。因此 sample/batch structure 明显是合理的，也正好说明下一步需要 PACS batch-effect feature removal。

### 4. 下一步如何证明 PACS 去除了 batch effect？

计划用 PACS 检测 batch-associated peaks，移除显著 peaks 后重新计算 TF-IDF / LSI / UMAP。然后比较 before/after 四联图：sample 分离应减弱，而 cell type 生物学结构应尽量保留。

### 5. AI 在项目中具体做了什么？

本项目是 human-AI collaborative reproduction workflow。AI 工具主要辅助路线规划、代码撰写、调试建议、报告组织和结果解释框架。本人负责复现目标选择、科学路线判断、关键 Linux/R 命令运行、结果检查、是否调整路线的决策，以及最终向老师汇报和解释。
