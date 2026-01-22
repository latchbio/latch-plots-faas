<goal>
Download study metadata and supplementary files from GEO. Collect paper context for downstream steps.
</goal>

<context_files>
Write large context to accession-specific folder (read by downstream steps):
```
/opt/latch/plots-faas/runtime/mount/agent_config/context/curation/tmp/{accession}/
├── study_metadata.txt      # GSE + SRA metadata
├── paper_text.txt          # User-provided paper abstract/methods/results
└── downloaded_files.txt    # List of downloaded file paths (one per line)
```
</context_files>

<method>
**Setup:**
```python
import sys
from pathlib import Path

sys.path.insert(0, "/opt/latch/plots-faas/runtime/mount/agent_config/context/curation/lib")
from curate.geo import gsm_to_gse

accession = "GSE12345"  # from user

# Convert GSM to GSE if needed (all functions expect GSE)
if accession.upper().startswith("GSM"):
    gse_id = gsm_to_gse(accession)
    if gse_id:
        print(f"Converted {accession} → {gse_id}")
        accession = gse_id
    else:
        print(f"Could not find parent GSE for {accession}")

CONTEXT_DIR = Path(f"/opt/latch/plots-faas/runtime/mount/agent_config/context/curation/tmp/{accession}")
CONTEXT_DIR.mkdir(parents=True, exist_ok=True)
```

**From GSE ID:**
1. Get study metadata and write to file:
   ```python
   from curate.geo import construct_study_metadata
   study_metadata = construct_study_metadata(accession)
   (CONTEXT_DIR / "study_metadata.txt").write_text(study_metadata)
   ```

2. List/download supplementary files:
   ```python
   from curate.geo import list_gse_supplementary_files, download_gse_supplementary_files
   files = list_gse_supplementary_files(accession)  # preview first
   downloaded_paths = download_gse_supplementary_files(accession, target_dir)
   (CONTEXT_DIR / "downloaded_files.txt").write_text("\n".join(str(p) for p in downloaded_paths))
   ```

3. Collect paper text via widget (NOT in chat - avoids history bloat):
   ```python
   from lplots.widgets.text import w_text_input

   paper_input = w_text_input(
       label="Paper Text",
       placeholder="Paste abstract, methods, and results here...",
       multiline=True,
   )
   ```
   Then in a separate cell after user fills it:
   ```python
   paper_text = paper_input.value
   if paper_text and paper_text.strip():
       (CONTEXT_DIR / "paper_text.txt").write_text(paper_text)
   ```
   Critical for cell count estimation, cell typing, metadata harmonization.
   **Never ask user to paste paper text in chat** - always use widget.

**User has data uploaded:**
- Skip GEO download
- Select files via w_ldata_picker, write paths to downloaded_files.txt
- Still collect paper text and any available metadata
</method>

<library>
curate.geo
</library>

<self_eval_criteria>
- Context directory exists: `/opt/latch/.../curation/tmp/{accession}/`
- `study_metadata.txt` exists in context dir (or confirmed unavailable)
- `downloaded_files.txt` exists with valid paths
- `paper_text.txt` exists (strongly recommended)
</self_eval_criteria>
