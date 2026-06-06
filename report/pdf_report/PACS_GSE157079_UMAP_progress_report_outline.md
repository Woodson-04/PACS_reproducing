# PACS 复现项目：GSE157079 UMAP 重建与 P56 批次效应去除报告

## 摘要

- 本报告总结 PACS 复现项目从 Notebook 1 benchmark 到 GSE157079 mouse kidney UMAP 重建的当前进展。
- Notebook 1 已证明 PACS 主 workflow 可以在 real kidney benchmark 上跑通。
- 当前重点结果是 P56 子集上的 PACS batch-effect feature filtering：在 13,526 个 P56 细胞和 top 10,000 peaks 中，PACS 检出 6,305 个 FDR 0.05 显著 batch-associated peaks，过滤后保留 3,695 个 peaks。
- 过滤后 batch mixing 在 UMAP、LSI 和 PCA-logNorm 空间均明显改善，同时 cell-type 结构仍然保留。
- 当前结果是 paper-style reconstruction 的阶段性版本，不声明为作者原图逐像素复现。

## 1. 任务背景与复现目标

- 项目目标是复现 PACS 论文中的关键分析流程，尤其是 mouse kidney 数据中的 batch-effect feature detection 和过滤前后 UMAP 展示。
- 前期 Notebook 1 benchmark 已完成，用于确认 PACS 主方法 Type I error / power workflow 能跑通。
- 当前阶段目标是从 GSE157079 原始 cell-by-peak matrix 出发，构建 matrix-derived UMAP，并在 P56 两个 batch 上演示 PACS 去除 batch-associated peaks 的效果。

## 2. 数据基础与子集选择

- 数据来源：GSE157079 snATAC。
- 已确认 cell-by-peak MatrixMarket 文件维度为 28,316 cells x 300,755 peaks，非零项 166,121,193。
- metadata 标准列：`row_index`, `cell_barcode`, `sample`, `cell_type`, `umap_1`, `umap_2`。
- P56 子集：
  - P56_batch1: 7,129 cells
  - P56_batch2: 6,397 cells
  - total P56 cells: 13,526
  - cell types: 15 类。

## 3. 从 GEO UMAP 到原始矩阵重建

- GEO 预计算 UMAP 只作为参考，不作为 PACS filtering 前后的 embedding 输入。
- 已完成从 MatrixMarket count matrix 出发的 TF-IDF / LSI / UMAP pipeline。
- all-cell top-peaks UMAP 证明可以从原始矩阵构建出具有生物学 cell-type 结构的 embedding。
- 后续 P56 before/after UMAP 均基于重新计算的 LSI/UMAP，不使用 GEO `umap_1` / `umap_2`。

## 4. P56 PACS batch peak filtering 流程

- 输入矩阵方向：原始文件为 cells x peaks；PACS 输入为 peaks x cells。
- 当前主设置：
  - `n_top_peaks = 10000`
  - `max_pacs_peaks = 10000`
  - `fdr_cutoff = 0.05`
  - full model: `~ cell_type + batch`
  - null model: `~ cell_type`
- 因 GSE157079 未提供 Notebook 1 中作者的 `q_vec`，当前使用 selected top peaks 的相对 depth-derived cap rates 作为 practical approximation。
- PACS 结果：
  - tested peaks: 10,000
  - significant batch peaks: 6,305
  - retained peaks: 3,695

## 5. UMAP 四联图结果

- 四联图展示：
  - before filtering colored by batch
  - before filtering colored by cell type
  - after filtering colored by batch
  - after filtering colored by cell type
- 视觉结论：
  - before 图中 P56_batch1 与 P56_batch2 有强烈分离。
  - after 图中 batch 分离明显减弱。
  - cell-type 结构仍然可见，说明过滤没有简单地把生物学结构全部抹掉。
- 推荐图：`fig_umap_four_panel_p56_10000.png`。

## 6. Paper-style normalized batch mixing score 数学定义

设第 `i` 个细胞的 k 近邻集合为 `N_k(i)`，batch label 为 `b_i`。

单细胞 batch mixing score:

```text
s_i = (1/k) sum_{j in N_k(i)} 1{b_j != b_i}
```

整体 observed score:

```text
S_obs = (1/n) sum_i s_i
```

设 cell type by batch 计数矩阵为：

```text
M = (m_ab)
```

其中 `a` 表示 cell type，`b` 表示 batch。考虑 cell-type composition 后的 expected score:

```text
S_exp = (1/n) sum_a sum_b m_ab [sum_{c != b} m_ac / sum_c m_ac]
```

归一化 score:

```text
S_norm = S_obs / S_exp
```

- `S_norm` 越高，表示 batch mixing 越接近期望随机混合。
- 当前报告同时计算 UMAP-space、LSI-space 和 PCA-logNorm space 的 score。

## 7. PCA-logNorm mixing score：最接近原文 normalized PCA mixing

- PACS 论文报告的是 normalized PCA mixing，因此 PCA-logNorm 分数是本项目中最接近原文 metric 的补充。
- 当前 PCA-logNorm 结果：
  - PC1:20: 0.0398 -> 0.2576
  - PC1:30: 0.0510 -> 0.3097
  - PC1:50: 0.0685 -> 0.4239
- 解释：
  - before 分数极低，说明原始 top-peak 空间 batch separation 很强。
  - after 分数明显提高，说明 PACS filtering 在 PCA-like 空间中减少了 batch signal。
  - PC 数越多，after 分数越高，提示部分 batch correction effect 分布在较高维成分中。

## 8. LSI-space 补充指标与 LSI_1-depth QC

- 对 snATAC 数据，TF-IDF/LSI 是比普通 PCA 更常用的表示空间，因此 LSI-space score 是 scATAC-aware analogue。
- 当前 LSI sensitivity 结果：
  - LSI_1:30: 0.0576 -> 0.2803
  - LSI_2:30: 0.0520 -> 0.3343
  - LSI_2:50: 0.0793 -> 0.3730
  - LSI_1:50: 0.0766 -> 0.3088
- LSI_1-depth QC：
  - before LSI_1 Spearman depth correlation = -0.9983
  - after LSI_1 Spearman depth correlation = 0.9976
- 因此 LSI_1 主要是 depth/QC component；使用 LSI_2:30 或 LSI_2:50 作为 batch mixing 评估更合理。

## 9. UMAP-space 与 GEO reference

- UMAP-space 是可视化空间指标，直观但不应作为唯一结论。
- 当前 UMAP-space paper-style score:
  - P56 10000/10000 before: 0.0265
  - P56 10000/10000 after: 0.6507
- GEO precomputed P56 UMAP reference:
  - normalized score = 0.9390
- 解释：
  - GEO UMAP 是作者/数据集提供的预计算坐标，可能包含完整 atlas 处理、不同特征选择或 batch handling。
  - 当前 PACS P56 result 不应直接要求达到 GEO reference score，但 after 明显优于 before。

## 10. 参数稳定性分析

- 已比较三组设置：
  - 5000/5000/FDR0.05
  - 10000/10000/FDR0.05
  - 20000/10000/FDR0.05
- 关键观察：
  - 5000 设置：significant peaks = 3,140，retained peaks = 1,860。
  - 10000 设置：significant peaks = 6,305，retained peaks = 3,695。
  - 20000/10000 设置：tested peaks = 10,000，significant peaks = 6,208，retained peaks = 3,792。
- 10000/10000 是当前主结果，因为它在 batch mixing 改善和 cell-type structure preservation 之间较平衡。

## 11. 当前亮点

- 已经从原始 MatrixMarket count matrix 走通 metadata alignment、sparse loading、TF-IDF/LSI/UMAP。
- PACS filtering 后 batch mixing 在多个空间中均改善：
  - PCA-logNorm PC1:50: 0.0685 -> 0.4239
  - LSI_2:30: 0.0520 -> 0.3343
  - UMAP-space: 0.0265 -> 0.6507
- LSI_1-depth QC 已验证，说明后续使用 LSI_2 开始的维度有明确依据。
- 参数比较显示 current main setting 不是偶然单点结果。

## 12. 局限性

- 当前是 P56-only 分析，不是完整 adult kidney all-batch figure。
- 当前只使用 top detected peaks，不是 all 300,755 peaks。
- cap rates 使用 depth-derived approximation，而不是作者原始 `q_vec`。
- 原始作者完整 mouse kidney figure code 尚未完全定位，因此当前是 PACS paper-style reconstruction，而非逐行复现。
- UMAP 图形状和长宽比不应与作者图做机械比较。

## 13. 下一步计划

- 将当前材料整理成正式中文 PDF 报告。
- 若时间允许，补充：
  - 更完整的 20000/20000 或 cell-type stratified PACS filtering。
  - 更接近作者设置的 all-feature 或 higher-top-peaks UMAP。
  - 对 retained/removed peaks 做 genomic annotation 或 motif enrichment。
- 汇报时强调：当前结果已证明 PACS workflow 在 GSE157079 P56 subset 上有效减少 batch-associated feature signal。

## 附录

- A. 主要文件路径清单：见 `materials/result_snapshot_manifest.md`。
- B. 图表来源说明：见 `materials/source_notes/report_source_notes.md`。
- C. UMAP 形状与长宽比说明：见 `materials/source_notes/umap_aspect_ratio_note.md`。
- D. 关键 CSV 表格：见 `materials/tables/`。
