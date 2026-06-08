# PACS GSE157079 / P56 UMAP 复现阶段报告（Draft）

**项目主题：** PACS 论文复现与 GSE157079 小鼠肾脏 snATAC 数据 batch-effect filtering 图形重建  
**当前定位：** P56 two-batch top10000 setting 下的 Fig.3a–d author-style reconstruction  
**重要说明：** 本报告不是 exact full Fig.3 reproduction；当前结果是基于 P56 two-batch subset、top10000 peaks 与 depth-derived `cap_rates` 的 PACS paper-style proof of concept。  
**版本：** Markdown draft，仅供审阅；尚未生成 PDF。

---

## 1. 标题页

本阶段工作围绕两个目标展开：

1. 复现 PACS 作者 Notebook 1 中 real kidney data 的 Type I error / power benchmark workflow；
2. 在公开 GSE157079 mouse kidney snATAC 数据上，构建从原始 cell-by-peak matrix 出发的 UMAP，并完成 P56 two-batch top10000 setting 下的 PACS batch-effect peak filtering 与 before/after UMAP 对照。

当前最重要结论：

> PACS filtering 在 P56 two-batch setting 中显著降低 batch-associated accessibility structure，并在 UMAP-space、LSI-space、PCA-logNorm-space 三类指标中均显示 batch mixing 改善，同时保留有意义的 cell-type structure。

---

## 2. 图表总览页

| 编号 | 图 / 表 | 内容 | 位置 |
|---|---|---|---|
| Fig. 1 | Notebook 1 benchmark | PACS Type I error / power | `materials/figures/pacs_benchmark_t1e_power_barplot.png` |
| Fig. 2 | Permuted-label QQ plot | Type I error 下 PACS p-value 行为 | `materials/figures/pacs_permuted_qq_plot.png` |
| Fig. 3 | all-cell UMAP by sample | 从原始矩阵重建的 all-cell UMAP | `materials/figures/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png` |
| Fig. 4 | all-cell UMAP by cell type | 同一 UMAP 的细胞类型结构 | `materials/figures/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png` |
| Fig. 5 | P56 PACS before/after UMAP | Fig.3a–d author-style 四联图 | `materials/figures/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png` |
| Fig. 6 | PCA/LSI/UMAP mixing score | batch mixing 定量比较 | `materials/figures/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png` |
| Fig. 7 | LSI-depth QC | LSI_1 与 depth 相关性 | `materials/figures/gse157079_p56_lsi_depth_correlation_before_after.png` |
| Fig. 8 | LSI dimension sensitivity | LSI 维度选择敏感性 | `materials/figures/gse157079_p56_lsi_dimension_sensitivity_batch_mixing.png` |
| Fig. 9 | GEO vs PACS metrics | GEO reference 与 PACS-filtered UMAP 对比 | `materials/figures/gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png` |

> 注：`gse157079_p56_pacs_setting_comparison_metrics.png` 当前未在项目图目录中找到，已记录于 `MISSING_REPORT_FIGURES.md`。参数稳定性结果以表格形式展示。

---

## 3. 摘要

Notebook 1 benchmark 复现显示，PACS 主方法在 real kidney data 上达到 Type I error = **0.04008**、power = **0.83337**，与作者 notebook 结果接近。由于作者 baseline helper 文件未找到，Seurat、ArchR、snapATAC、Fisher baseline 结果均标注为 clean-room reimplemented baselines。

GSE157079 数据接入阶段确认公开矩阵为 **28316 cells x 300755 peaks**，包含 **166121193** 个 nonzero entries，metadata 与 peak list 均与矩阵维度对齐。随后从原始矩阵直接完成 all-cell TF-IDF / LSI / UMAP 重建，使用 top **20000** detected peaks，保留 **85801336** 个 nonzero entries。

在 P56 subset 中，P56_batch1 = **7129** cells，P56_batch2 = **6397** cells，总计 **13526** cells。PACS 在 top10000 peaks 中测试 batch effect，使用 full model `~ cell_type + batch` 与 null model `~ cell_type`，识别 **6305** 个 FDR-significant batch peaks，保留 **3695** peaks 用于 after-filtering UMAP。过滤后 batch mixing 在 PCA-logNorm、LSI-space 与 UMAP-space 中均显著改善。

---

## 4. 复现目标与边界

本报告关注的是 **PACS paper-style reconstruction**，不是完整重建原始 dev-kidney-snATAC atlas pipeline。

当前复现目标：

- Notebook 1 PACS benchmark workflow；
- GSE157079 原始矩阵接入与 matrix-derived UMAP；
- P56 two-batch top10000 setting 下的 PACS batch-effect peak filtering；
- Fig.3a–d author-style before/after UMAP 四联图；
- normalized batch mixing score 的 PCA/LSI/UMAP 三空间定量。

当前边界：

- 不声称 exact full Fig.3 reproduction；
- 不声称使用作者完整 preprocessing 参数；
- 不声称 GEO precomputed UMAP 是 PACS-filtered UMAP；
- 不声称 baseline methods 完全复现作者原始 baseline 代码；
- GSE157079 缺少作者 Notebook 1 使用的 `q_vec`，因此本阶段使用 depth-derived `cap_rates` 作为 practical approximation。

---

## 5. Notebook 1 Benchmark

Notebook 1 large benchmark 配置：

- `n_repeat = 5`
- `n_cell_sample = 500`
- `n_features_sample = 10000`
- `run_baselines = TRUE`

| 方法 | Type I error | Power | 说明 |
|---|---:|---:|---|
| PACS / our | 0.04008 | 0.83337 | PACS 主方法 |
| Seurat | 0.06342 | 0.82344 | clean-room baseline |
| ArchR | 0.04096 | 0.67437 | clean-room approximation |
| snapATAC | 0.01810 | 0.76094 | clean-room edgeR-style baseline |
| Fisher | 0.02208 | 0.76630 | binary Fisher exact test |

![Notebook 1 benchmark](materials/figures/pacs_benchmark_t1e_power_barplot.png)

![Permuted-label QQ plot](materials/figures/pacs_permuted_qq_plot.png)

结论：

- PACS Type I error 接近 0.05；
- PACS power 接近作者 notebook；
- baseline 由于作者 helper 文件缺失，只能作为 clean-room reimplemented comparison；
- 主方法 PACS workflow 已足以支持后续 mouse kidney figure reconstruction。

---

## 6. GSE157079 数据接入与矩阵对齐

GSE157079 数据包括 metadata、UMAP coordinate、peak list 与 cell-by-peak MatrixMarket sparse matrix。矩阵对齐 smoke test 确认：

| 数据对象 | 数量 / 维度 |
|---|---:|
| cells | 28316 |
| peaks | 300755 |
| nonzero entries | 166121193 |
| metadata rows | 28316 |
| peak list rows | 300755 |
| matrix orientation | cell x peak |

标准化后的 metadata 字段：

- `row_index`
- `cell_barcode`
- `sample`
- `cell_type`
- `umap_1`
- `umap_2`

关键样本规模：

| sample | cells |
|---|---:|
| P56_batch1 | 7129 |
| P56_batch2 | 6397 |
| P0_batch1 | 5993 |
| P0_batch2 | 5436 |
| P21_batch1 | 3361 |

结论：

- GSE157079 四个公开文件可支持 core matrix-derived UMAP reconstruction；
- 精确复现作者图仍需要作者 preprocessing、feature filtering、embedding 参数等细节。

---

## 7. 从原始矩阵重建 All-Cell UMAP

为避免依赖 GEO precomputed UMAP，本项目从原始 cell-by-peak sparse matrix 出发，构建 TF-IDF / LSI / UMAP embedding。

| 项目 | 数值 |
|---|---:|
| cells | 28316 |
| selected top peaks | 20000 |
| retained nonzeros | 85801336 |
| sparse matrix dimension | 28316 x 20000 |
| LSI dimension | 28316 x 50 |
| UMAP input | LSI 2:30 |

![All-cell matrix-derived UMAP by sample](materials/figures/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png)

![All-cell matrix-derived UMAP by cell type](materials/figures/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png)

结果解释：

- sample-colored UMAP 显示明显 sample-associated structure；
- cell-type-colored UMAP 保留生物学细胞类型结构；
- 该结果作为 before-PACS-filtering 的全细胞参考，而不是最终 PACS-filtered 图。

---

## 8. P56 子集与 PACS Batch-Effect Peak Filtering

选择 P56 subset 的原因是 P56 有两个 batch，且两个 batch cell 数量均充足，适合进行 same-age two-batch batch-effect test。

| 项目 | 数值 |
|---|---:|
| P56 cells | 13526 |
| P56_batch1 cells | 7129 |
| P56_batch2 cells | 6397 |
| tested peaks | 10000 |
| significant batch peaks | 6305 |
| retained peaks | 3695 |
| full model | `~ cell_type + batch` |
| null model | `~ cell_type` |

PACS 模型解释：

- 原始矩阵方向为 cells x peaks；
- PACS 输入方向为 peaks x cells；
- full model 检验 cell type 调整后的 batch effect；
- null model 保留 cell type 结构；
- 显著 batch peaks 被移除后，用 retained peaks 重新构建 LSI / UMAP。

当前设置：

```text
n_top_peaks = 10000
max_pacs_peaks = 10000
fdr_cutoff = 0.05
```

---

## 9. Fig.3a–d Author-Style UMAP 四联图

下图是当前最核心的 P56 two-batch top10000 setting 下的 Fig.3a–d author-style reconstruction：

![P56 PACS before/after UMAP four-panel](materials/figures/gse157079_p56_pacs_batch_filter_p56_top10000_test10000_fdr005_lsi_saved_four_panel.png)

图像阅读方式：

- before by batch：P56_batch1 与 P56_batch2 在 UMAP 中存在强 batch separation；
- after by batch：移除 PACS-significant batch peaks 后，batch separation 明显减弱；
- before / after by cell type：用于判断生物学 cell-type structure 是否被过度破坏。

当前图应表述为：

> P56 two-batch top10000 setting 下的 Fig.3a–d author-style reconstruction。

不应表述为：

> exact full Fig.3 reproduction。

---

## 10. Normalized Batch Mixing Score 数学定义

对每个 cell \(i\)，设其 batch label 为 \(b_i\)，其 \(k\)-nearest neighbors 为 \(N_k(i)\)。cell-level batch mixing score 定义为：

$$
s_i = \frac{1}{k}\sum_{j \in N_k(i)} I(b_j \ne b_i)
$$

全局 observed score 为：

$$
S_{\mathrm{obs}} = \frac{1}{n}\sum_i s_i
$$

考虑 cell type composition 后，expected score 定义为：

$$
S_{\mathrm{exp}} =
\frac{1}{n}
\sum_a \sum_b
m_{ab}
\left(
\frac{\sum_{d \ne b} m_{ad}}{\sum_d m_{ad}}
\right)
$$

其中 \(m_{ab}\) 表示 cell type \(a\)、batch \(b\) 中的 cell 数。

normalized batch mixing score 定义为：

$$
S_{\mathrm{norm}} = \frac{S_{\mathrm{obs}}}{S_{\mathrm{exp}}}
$$

解释：

- \(S_{\mathrm{norm}}\) 越高，batch mixing 越接近期望混合状态；
- \(S_{\mathrm{norm}}\) 不应单独解释，应同时考虑 cell type preservation；
- PCA-space score 最接近 PACS paper 的 normalized PCA mixing 思路；
- LSI-space 是 scATAC-aware analogue；
- UMAP-space 是 visualization-level approximation。

---

## 11. PCA-LogNorm 主指标

PCA-logNorm 是当前最接近 PACS paper normalized PCA mixing 的 author-like metric。

| PCA dimensions | before | after |
|---|---:|---:|
| PC1:20 | 0.0398 | 0.2576 |
| PC1:30 | 0.0510 | 0.3097 |
| PC1:50 | 0.0685 | 0.4239 |

![PCA/LSI/UMAP normalized batch mixing comparison](materials/figures/gse157079_p56_10000_pca_lsi_umap_paper_style_score_comparison.png)

结论：

- PC1:20、PC1:30、PC1:50 三个读数均显示 after-filtering score 明显高于 before；
- 这支持 PACS filtering 的改善不是 UMAP visualization artifact；
- PC1:30 可作为主读数，PC1:50 可作为高维敏感性支持。

---

## 12. LSI-Space 与 LSI_1-Depth QC

scATAC 数据通常使用 TF-IDF / LSI 表示稀疏 chromatin accessibility matrix。因此，LSI-space score 是 PCA-logNorm 之外的重要支持指标。

| 坐标空间 | 维度 / 设置 | before | after |
|---|---:|---:|---:|
| LSI-space | LSI 2:30 | 0.0520 | 0.3343 |

LSI_1-depth QC：

| 阶段 | component | Spearman(depth) | Pearson(depth) |
|---|---|---:|---:|
| before | LSI_1 | -0.9983 | -0.9739 |
| after | LSI_1 | 0.9976 | 0.9830 |

![LSI depth correlation QC](materials/figures/gse157079_p56_lsi_depth_correlation_before_after.png)

解释：

- LSI_1 与 depth 几乎完全相关；
- 使用 LSI 2:30 计算 batch mixing score 是合理的；
- LSI-space 的 before -> after 改善支持 PACS filtering 在 scATAC-aware high-dimensional space 中有效。

LSI 维度敏感性：

![LSI dimension sensitivity](materials/figures/gse157079_p56_lsi_dimension_sensitivity_batch_mixing.png)

---

## 13. UMAP-Space 与 GEO Reference

UMAP-space score 是可视化层面的支持指标：

| setting | normalized batch mixing score |
|---|---:|
| P56 10000/10000 before | 0.02649 |
| P56 10000/10000 after | 0.65073 |
| GEO precomputed UMAP P56 | 0.93903 |

![GEO reference vs P56 PACS UMAP metrics](materials/figures/gse157079_geo_vs_p56_pacs_umap_metrics_comparison.png)

解释：

- P56 PACS-filtered UMAP 的 UMAP-space score 明显高于 before；
- GEO P56 precomputed UMAP score 更高，但该 UMAP 可能来自不同 atlas pipeline 或 integration/correction 流程；
- GEO reference 可作为 public reference embedding，不应等同于 PACS-filtered UMAP。

---

## 14. 参数稳定性分析

当前比较了 5000/5000/FDR0.05 与 10000/10000/FDR0.05 两个设置。

| setting | significant peaks | retained peaks | normalized batch entropy | same-batch fraction | cell type silhouette | same-celltype fraction |
|---|---:|---:|---:|---:|---:|---:|
| 5000/5000/FDR0.05 | 3140 | 1860 | 0.69564 | 0.69693 | 0.15933 | 0.66504 |
| 10000/10000/FDR0.05 | 6305 | 3695 | 0.71111 | 0.68101 | 0.30198 | 0.77592 |

结论：

- 10000/10000/FDR0.05 相比 5000/5000/FDR0.05 更优；
- 它改善 batch mixing，同时更好保留 cell-type structure；
- 因此当前主结果采用 top10000 / test10000 / FDR0.05。

---

## 15. 当前复现亮点

1. Notebook 1 PACS 主方法已成功复现，Type I error 与 power 接近作者结果。
2. GSE157079 原始矩阵、metadata、peak list 已验证对齐。
3. 已从原始 MatrixMarket sparse matrix 出发重建 all-cell TF-IDF / LSI / UMAP。
4. P56 two-batch setting 中完成 PACS batch-effect peak filtering。
5. before/after UMAP 四联图已形成清晰 author-style reconstruction。
6. PCA-logNorm、LSI-space、UMAP-space 三类 normalized batch mixing score 均支持 PACS filtering 改善 batch mixing。
7. LSI_1-depth QC 解释了为什么使用 LSI 2:30。

---

## 16. 局限性

当前结果仍有以下限制：

- P56-only two-batch subset，不是 full all-age dataset；
- top10000 peaks setting，不是 all 300755 peaks；
- depth-derived `cap_rates` 是 practical approximation，不是作者原始 `q_vec`；
- Notebook 1 baseline methods 是 clean-room reimplementation；
- GEO precomputed UMAP 不是 PACS-filtered UMAP；
- 当前报告展示的是 Fig.3a–d author-style reconstruction，不是 exact full Fig.3 reproduction。

---

## 17. 下一步计划

建议下一阶段按优先级推进：

1. 审阅本 Markdown draft，确认术语与结论表述；
2. 补齐缺失或未复制的 report material figures；
3. 整理并统一 figure captions；
4. 若需要更强结果，可扩展至更多 top peaks 或 cell-type stratified PACS filtering；
5. 在 Markdown 审阅通过后，再生成 PDF；
6. 若用于论文级展示，进一步补充 exact author preprocessing 的证据或说明。

---

## 18. 附录 A：核心文件索引

核心结果文件：

- `results/current_important_results_summary.md`
- `results/20260526_2318_large_baseline/summary.csv`
- `results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md`
- `results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md`
- `results/mouse_kidney_figures/gse157079_p56_pacs_batch_filter_umap_p56_top10000_test10000_fdr005_lsi_saved/p56_pacs_batch_filter_umap_report.md`
- `results/mouse_kidney_figures/paper_style_batch_mixing_score_pca_space/p56_10000_pca_space_score_report.md`
- `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_paper_style_metric_final_interpretation.md`
- `results/mouse_kidney_figures/paper_style_batch_mixing_score_lsi_space/p56_lsi_depth_correlation_report.md`
- `results/mouse_kidney_figures/geo_precomputed_umap_quantification/geo_precomputed_umap_quantification_report.md`

报告材料：

- `materials/tables/`
- `materials/figures/`
- `materials/source_notes/report_source_notes.md`
- `MISSING_REPORT_FIGURES.md`

---

## 19. 附录 B：AI 协作说明

本项目中 AI 协作主要用于：

- 阅读和整理已有 R 脚本、结果表和 Markdown 报告；
- 根据作者 notebook 逻辑修复 PACS 主流程中的输入方向、p-value 合并和输出记录问题；
- 在不修改原始数据和 PACS 安装包的前提下，构建本地 wrapper 与分析脚本；
- 整理 GSE157079 数据接入、matrix-derived UMAP、P56 PACS filtering 和 batch mixing quantification 的结果报告；
- 生成当前 Markdown draft。

所有关键数值均来自当前项目已有输出文件。本 draft 未重新运行 PACS、UMAP 或 MatrixMarket streaming，也未修改源数据。

