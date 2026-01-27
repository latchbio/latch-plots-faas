# SpatialBench Eval Rubric (Shared)

This section is **authoritative** and is intended to be reused verbatim in:
- the human-facing curriculum (as the canonical rubric), and
- the judge system prompt (as the grading rubric).

## Core Rule: Specify WHAT, not HOW
- For all evals: tasks should specify **what to return** (exact JSON fields) and the **scientific goal**, not a step-by-step method or parameters.
- Exception: **PROCEDURAL** evals must name the method/technique to apply (but still should not fix exact parameters unless necessary).

## "Notes" are not leakage
- The `notes` field is internal and **not shown to the agent**.
- Never treat notes content as an anti-shortcut leak.

---

## Principle 1 — Verifiability (objective grading possible)

**PASS** when all are true:
- Task specifies an **exact JSON output format** (required field names + expected value types).
- Grader type matches the output structure (e.g., single number → numeric tolerance; single letter → multiple choice).
- Success condition is automatically checkable; no subjective interpretation required.

**FAIL** when any are true:
- Output format is missing, ambiguous, or does not list required fields.
- Grader type does not match the expected output structure.
- Task relies on subjective language ("interesting", "meaningful", "significant", "interpret") **without** an objective operational definition.

**Do NOT fail solely because**:
- thresholds are not specified (that can be good for anti-shortcut),
- algorithm choice is not specified (SCIENTIFIC evals intentionally allow this),
- the task does not explain "how to identify" features (that can preserve anti-shortcut).

---

## Principle 2 — Scientific Durability (robust to reasonable choices OR tolerance covers variation)

Durability is judged relative to `metadata.eval_type`:

### SCIENTIFIC evals
**Goal:** Agent chooses methods + parameters to achieve a scientific goal.
**PASS** if:
- reasonable methods should converge on the same answer **OR**
- tolerances are wide enough to accept valid scientific variation.
**FAIL** if:
- reasonable methodological variation yields materially different answers and tolerance is too tight, or
- task depends on brittle artifacts (random seed, library version quirks, arbitrary coordinate systems).

**Red flag:** scientific + open-ended phrasing + tight numeric tolerance.

### PROCEDURAL evals
**Goal:** Agent applies a named method correctly; parameters not specified.
**PASS** if:
- method is explicitly named,
- tolerance accommodates reasonable parameter choices within that method.
**FAIL** if:
- method is not named (then it's not procedural),
- task over-specifies exact parameters/steps in a way that tests implementation trivia,
- output depends on unstable artifacts (seed, exact coordinate values).

### OBSERVATIONAL evals
**Goal:** Exposure / trajectory generation / loose testing.
Durability is relaxed. Primary focus becomes verifiability + anti-shortcut, with only obvious brittleness flagged.

---

## Principle 3 — Anti-Shortcut Structure (requires data interaction)

**PASS** if:
- answer cannot be produced from prior knowledge alone,
- task text/options do not leak the answer,
- solving requires inspecting and computing on the provided dataset.

**FAIL** if:
- it's a textbook/definition question (dataset irrelevant),
- answer is leaked in task text, answer options, or easily accessible precomputed dataset fields,
- the agent can read the answer directly from stored results (precomputed embeddings/labels/summaries) rather than recomputing.

**Important:** Don't fail simply because the output is "just a number." Dataset-specific numbers are not guessable if shortcuts are blocked.

**Anti-shortcut hardening (preferred):**
- Remove precomputed embeddings and downstream artifacts (e.g., `adata.obsm["X_pca"]`, `adata.obsm["X_umap"]`, precomputed cluster labels, cached marker summaries) when they would enable reading answers without computation.
- For multiple-choice, ensure distractors are biologically plausible and not label-leaking ("Cluster 3 (bone formation)" is a leak; "Cluster 3" is fine).

---

## Grader ↔ Output Shape Compatibility (minimum expectations)

A grader must align with what the task asks the agent to return.

### numeric_tolerance
- Expected output: JSON numeric fields, e.g. `{"spots_after_filtering": 2947}`
- Used when the answer is a number with acceptable variation.

### multiple_choice
- Expected output: `{"answer": "A"}` (case-insensitive)
- Used when a single option letter is the answer; distractors must be plausible.

### distribution_comparison
- Expected output: `{"cell_type_distribution": {"TypeA": <pct>, ...}}`
- Used when multiple class proportions are graded together.

### marker_gene_precision_recall
- Expected output: `{"top_marker_genes": ["GENE1", ...]}`
- Used for ranked marker/DE lists; scoring typically via recall/precision thresholds.

### label_set_jaccard
- Expected output: `{"labels": ["A", "C"]}` (order-independent)
- Used for multi-select sets; penalizes missing + extra items.

---

## Using Trajectories (if provided) as evidence

Trajectories can indicate:
- **Shortcutting:** very short "answer immediately" trajectories + high pass rate.
- **Tolerance too tight:** correct-looking work + low pass rate.
- **Durability issues:** large variance across runs/methods.

Only cite trajectory evidence that is explicitly visible in the trajectory.

---

## Suggestion priority (actionable fixes)

When recommending fixes, prioritize in this order:
1. If shortcuts exist via precomputed dataset fields: **remove those fields first**.
2. Fix verifiability: specify JSON output fields; align grader with output.
3. Fix durability: widen/adjust tolerances; anchor inputs; avoid unstable artifacts; adjust eval_type if mismatched.
4. Fix anti-shortcut: remove answer leaks; make distractors plausible; avoid "answer in question".
