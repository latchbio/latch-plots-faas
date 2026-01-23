<goal>
Download study metadata and supplementary files from GEO. Collect paper context for downstream steps.
</goal>

<method>
**1. Create paper text widget first (do not skip):**
```python
from lplots.widgets.text import w_text_input

paper_input = w_text_input(
    label="Paper Text",
    placeholder="Paste abstract, methods, and results here...",
    multiline=True,
)
```
Tell user: "Paste paper text (abstract/methods/results) in the widget while I download the GEO data."

**2. Setup and convert GSM→GSE if needed:**
```python
import sys
from pathlib import Path
sys.path.insert(0, "/opt/latch/plots-faas/runtime/mount/agent_config/context/curation/lib")
from curate.geo import gsm_to_gse, construct_study_metadata, download_gse_supplementary_files

if accession.upper().startswith("GSM"):
    accession = gsm_to_gse(accession) or accession

CONTEXT_DIR = Path(f"/opt/latch/plots-faas/runtime/mount/agent_config/context/curation/tmp/{accession}")
CONTEXT_DIR.mkdir(parents=True, exist_ok=True)
```

**3. Download metadata and files:**
```python
(CONTEXT_DIR / "study_metadata.txt").write_text(construct_study_metadata(accession))
downloaded_paths = download_gse_supplementary_files(accession, target_dir)
(CONTEXT_DIR / "downloaded_files.txt").write_text("\n".join(str(p) for p in downloaded_paths))
```

**4. Save paper text (after user fills widget):**
```python
if paper_input.value and paper_input.value.strip():
    (CONTEXT_DIR / "paper_text.txt").write_text(paper_input.value)
```
</method>

<workflows>
</workflows>

<library>
curate.geo
</library>

<self_eval_criteria>
- Paper text widget displayed and user informed to fill it
- `study_metadata.txt` written with platform info (check for platform_id field)
- `downloaded_files.txt` written with valid file paths
- `paper_text.txt` written (confirm user provided input)
</self_eval_criteria>
