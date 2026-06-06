# UMAP 形状与长宽比说明

当前 UMAP 图的整体形状和长宽比例不应被过度解读。UMAP 坐标本身没有固定方向、尺度或长宽比；同一批点可以因为绘图面板宽高、坐标范围、padding、是否使用 `coord_equal()`、是否裁剪边界、UMAP 初始化和参数不同，而呈现出更扁、更宽或更接近正方形的视觉效果。

作者论文中的 UMAP 面板可能采用了更接近平方面板的排版，也可能使用了不同的 UMAP 初始化、特征集合、细胞范围或 batch-correction 后坐标。因此，我们当前的 P56 子集、top 10000 peaks、PACS retained peaks UMAP 不需要也不应该强行匹配作者完整 adult kidney 多批次 UMAP 的外形。

如果用于汇报展示，可以在不改变坐标和定量指标的前提下，用更方正的画布、`coord_equal()`、统一 padding 或裁剪方式改善视觉呈现。但仅为了让图“更像正方形”而重新运行 UMAP 没有必要；真正关键的是 PCA/LSI/UMAP 空间中的 batch mixing 定量指标和 cell-type structure 是否被保留。
