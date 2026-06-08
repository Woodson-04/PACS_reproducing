# Final Report Manifest

## Final Outputs

| item | status | path |
|---|---|---|
| Final Markdown | generated | `report/pdf_report/PACS_GSE157079_UMAP_progress_report_FINAL.md` |
| Final HTML | generated | `report/pdf_report/PACS_GSE157079_UMAP_progress_report_FINAL.html` |
| Final PDF | not generated in current environment | `report/pdf_report/PACS_GSE157079_UMAP_progress_report_FINAL.pdf` |
| CSS | generated | `report/pdf_report/report_style.css` |
| Figure copy script | generated | `report/pdf_report/COPY_REPORT_FIGURES.sh` |
| Source notes | generated | `report/pdf_report/materials/source_notes/report_source_notes.md` |

## Embedded Materials

- Required figures present in `materials/figures/`: 9 / 9
- Core Markdown tables in `materials/tables/`: 7

## HTML Generation Method

HTML was generated as a standalone local HTML file with:

- `report_style.css`
- MathJax support for LaTeX formulas
- local image references under `materials/figures/`

Because `pandoc` was unavailable, HTML was generated via a controlled local
HTML rendering step rather than a Pandoc pipeline.

## PDF Generation Method Attempted

Attempted method:

- Microsoft Edge headless

Outcome:

- Failed due to Windows PowerShell / SSHFS access-denied errors in Edge
  crashpad / mojo / profile startup.
- No final PDF was produced.

## Cleanup Actions

Cleanup was **not performed**, because PDF generation did not succeed. The
following files were therefore intentionally retained:

- `PACS_GSE157079_UMAP_progress_report_DRAFT.md`
- `COPY_REPORT_FIGURES_INSTRUCTIONS.md`
- `GENERATE_FINAL_PDF_INSTRUCTIONS.md`
- `MISSING_REPORT_FIGURES.md`
- `render_markdown_to_html_temp.ps1`

## Warnings

- Chinese text exists in the HTML source.
- MathJax support exists in the HTML source.
- Formula rendering was not visually verified in a browser by this run.
- PDF Chinese rendering and formula rendering were not verified because PDF
  generation failed in the current environment.
