# Report Figure Status

All required source figures for the final draft were found under
`figures/mouse_kidney/`.

The figure files have **not** been copied into
`report/pdf_report/materials/figures/` in the current Codex/PowerShell/SSHFS
environment, because binary file copying previously failed with access-denied
errors.

Please run the Linux copy script manually:

```bash
cd /home/woodson/PACS_reproducing
bash report/pdf_report/COPY_REPORT_FIGURES.sh
```

After that, `report/pdf_report/materials/figures/` should contain the required
PNG files and the PDF can be generated safely.
