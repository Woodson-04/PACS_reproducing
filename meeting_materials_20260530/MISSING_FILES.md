# Missing Or Not Copied Files

Source-file existence was checked earlier, and the requested meeting source
files are present in the project. However, this Codex/PowerShell SSHFS session
was denied permission when copying binary and text files into
`meeting_materials_20260530`.

The directory structure and this guide were created. If the copied meeting
files are not present, run the copy commands from a Linux terminal inside:

```bash
cd /home/woodson/PACS_reproducing
mkdir -p meeting_materials_20260530/reports meeting_materials_20260530/figures meeting_materials_20260530/tables

cp results/project_progress_summary_for_teacher.md meeting_materials_20260530/reports/
cp results/mouse_kidney_figures/gse157079_matrix_alignment_smoke_test.md meeting_materials_20260530/reports/
cp results/mouse_kidney_figures/gse157079_all_cells_top_peaks_lsi_umap/all_cells_top_peaks_lsi_umap_report.md meeting_materials_20260530/reports/
cp notebook1_reproduction_report.md meeting_materials_20260530/reports/
cp baseline_official_methods_review.md meeting_materials_20260530/reports/

cp figures/mouse_kidney/pacs_benchmark_t1e_power_barplot.png meeting_materials_20260530/figures/
cp figures/mouse_kidney/pacs_permuted_qq_plot.png meeting_materials_20260530/figures/
cp figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_sample.png meeting_materials_20260530/figures/
cp figures/mouse_kidney/gse157079_all_cells_top_peaks_lsi_umap_by_celltype.png meeting_materials_20260530/figures/
cp figures/mouse_kidney/gse157079_umap_by_sample.png meeting_materials_20260530/figures/
cp figures/mouse_kidney/gse157079_umap_by_celltype.png meeting_materials_20260530/figures/

cp results/20260526_2318_large_baseline/summary.csv meeting_materials_20260530/tables/
cp results/mouse_kidney_figures/gse157079_metadata_summary.csv meeting_materials_20260530/tables/
```
