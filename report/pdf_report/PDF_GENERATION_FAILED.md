# PDF Generation Failed

Final Markdown and HTML were generated successfully, but PDF generation failed
in the current Windows PowerShell / SSHFS environment.

## Generated Files

- `report/pdf_report/PACS_GSE157079_UMAP_progress_report_FINAL.md`
- `report/pdf_report/PACS_GSE157079_UMAP_progress_report_FINAL.html`
- `report/pdf_report/report_style.css`

## Figure Status

All 9 required report figures are present under:

```text
report/pdf_report/materials/figures/
```

## Attempted PDF Method

Microsoft Edge headless was attempted:

```text
C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe
```

The browser failed before PDF output was written.

## Error Summary

The observed errors included:

```text
CreateFile: access denied
crash server failed to launch
mojo platform channel check failed
Failed to create a ProcessSingleton for profile directory
```

Specifying a temporary user-data-dir did not resolve the issue in this
environment.

## Recommended Linux Commands

Run PDF generation in a Linux terminal:

```bash
cd /home/woodson/PACS_reproducing
bash report/pdf_report/COPY_REPORT_FIGURES.sh
```

If Chromium is installed:

```bash
cd /home/woodson/PACS_reproducing/report/pdf_report
chromium-browser \
  --headless \
  --disable-gpu \
  --no-sandbox \
  --print-to-pdf=/home/woodson/PACS_reproducing/report/pdf_report/PACS_GSE157079_UMAP_progress_report_FINAL.pdf \
  --print-to-pdf-no-header \
  file:///home/woodson/PACS_reproducing/report/pdf_report/PACS_GSE157079_UMAP_progress_report_FINAL.html
```

If Chromium is unavailable, install a PDF route:

```bash
sudo apt update
sudo apt install -y chromium-browser
```

or:

```bash
sudo apt install -y pandoc wkhtmltopdf
```

Then verify:

```bash
ls -lh /home/woodson/PACS_reproducing/report/pdf_report/PACS_GSE157079_UMAP_progress_report_FINAL.pdf
```
