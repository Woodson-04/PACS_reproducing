# Generate Final PDF

The final Markdown is ready:

```text
report/pdf_report/PACS_GSE157079_UMAP_progress_report_FINAL.md
```

Run these commands in a Linux terminal after copying figures:

```bash
cd /home/woodson/PACS_reproducing
bash report/pdf_report/COPY_REPORT_FIGURES.sh
```

Preferred option with Quarto:

```bash
cd /home/woodson/PACS_reproducing/report/pdf_report
quarto render PACS_GSE157079_UMAP_progress_report_FINAL.md --to pdf
```

Alternative with Pandoc and XeLaTeX:

```bash
cd /home/woodson/PACS_reproducing/report/pdf_report
pandoc PACS_GSE157079_UMAP_progress_report_FINAL.md \
  --pdf-engine=xelatex \
  -V mainfont="Noto Sans CJK SC" \
  -V CJKmainfont="Noto Sans CJK SC" \
  -o PACS_GSE157079_UMAP_progress_report_FINAL.pdf
```

If the Chinese font is unavailable, try:

```bash
fc-list :lang=zh | head
```

Then replace `Noto Sans CJK SC` with an available Chinese font.

Verify:

```bash
ls -lh PACS_GSE157079_UMAP_progress_report_FINAL.pdf
```
