# Archive Manifest

This manifest was prepared before moving files. In the current Codex/PowerShell
SSHFS session, move operations requiring new archived file creation were denied
by the sandbox/permission layer, so no files were moved by this assistant.

Recommended archive candidates after confirmation:

| original path | archived path | reason |
|---|---|---|
| `results/kidney_notebook1_20260526_172547/` | `archive_pre_meeting_cleanup_20260530/results/kidney_notebook1_20260526_172547/` | early Notebook 1 small/debug result, not meeting-facing |
| `results/kidney_notebook1_20260526_205633/` | `archive_pre_meeting_cleanup_20260530/results/kidney_notebook1_20260526_205633/` | early baseline/debug result |
| `results/kidney_notebook1_20260526_221904/` | `archive_pre_meeting_cleanup_20260530/results/kidney_notebook1_20260526_221904/` | medium baseline development history |
| `cleanup_manifest.md` | `archive_pre_meeting_cleanup_20260530/cleanup_manifest.md` | old cleanup planning document replaced by current retention recommendation |

Suggested Linux commands if confirmed:

```bash
cd /home/woodson/PACS_reproducing
mkdir -p archive_pre_meeting_cleanup_20260530/results
mv results/kidney_notebook1_20260526_172547 archive_pre_meeting_cleanup_20260530/results/
mv results/kidney_notebook1_20260526_205633 archive_pre_meeting_cleanup_20260530/results/
mv results/kidney_notebook1_20260526_221904 archive_pre_meeting_cleanup_20260530/results/
mv cleanup_manifest.md archive_pre_meeting_cleanup_20260530/
```
