# Output Structure Refactor Plan

> Branch scope: `main` (pure EEG MATLAB/EEGLAB pipeline)
>
> Goal: redesign output organization so that subject-level results, batch-level merged data, analysis task outputs, and paper-delivery artifacts are clearly separated, easier to navigate, and easier to maintain.

---

## 1. Why refactor the current output structure?

The current output layout is functional, but it is starting to carry too many responsibilities at once.

### Current pain points

1. **Mixed responsibilities under `summary/`**
   - merged tables
   - QC outputs
   - analysis-2 task outputs
   - paper figures
   - reports
   - methods snapshots

2. **Single-subject and group-level outputs use different organizational logic**
   - subject-level outputs are grouped by subject
   - batch/task outputs are grouped mostly by historical script naming

3. **`raw` vs `qc` is not consistently represented as a directory-level concept**
   - sometimes embedded in filenames
   - sometimes represented as subdirectories
   - sometimes implicit

4. **Paper-ready outputs are not fully separated from analytical intermediate outputs**
   - `paper_fig/` behaves partly like a deliverable folder, partly like an analysis byproduct

5. **Future growth risk**
   - adding Task8/Task9 or additional figure families will further overload `summary/`
   - downstream users must already know repository internals to find the right file quickly

---

## 2. Refactor goals

The new output structure should satisfy the following goals:

### G1. Clear semantic layers
Separate outputs into four levels:
- subject-level
- batch-level
- analysis-level
- paper-delivery-level

### G2. Consistent location rules
A user should be able to infer where a file lives from its role, without memorizing script history.

### G3. Stronger `raw` / `qc` distinction
For analysis results, `raw` and `qc` should be first-class directory layers, not just filename suffixes.

### G4. Backward-compatible migration path
Avoid breaking all downstream scripts at once.

### G5. Support future task expansion
Task directories should scale naturally:
- task1_block2_restart
- task2_scene_block_diff
- task3_trialindex_lmm
- task4_core_lmm_suite
- task5_peakindex_invertedu
- task6_obeta_special
- task7_individual_checks
- future task8/task9/...

---

## 3. Proposed target structure

Recommended top-level layout:

```text
bandpower_outputs/
  runs/
    run_YYYYMMDD_HHMMSS/
      subjects/
        <subject_id>/
          tables/
          figures/
          qc/
          report/
      batch/
        merged/
        qc/
        reports/
        audit/
      analysis/
        task1_block2_restart/
          raw/
            tables/
            figures/
            reports/
          qc/
            tables/
            figures/
            reports/
        task2_scene_block_diff/
          raw/
            tables/
            figures/
            reports/
          qc/
            tables/
            figures/
            reports/
        task3_trialindex_lmm/
          raw/
            tables/
            figures/
            reports/
          qc/
            tables/
            figures/
            reports/
        task4_core_lmm_suite/
          raw/
            tables/
            figures/
            reports/
          qc/
            tables/
            figures/
            reports/
        task5_peakindex_invertedu/
          raw/
            tables/
            figures/
            reports/
          qc/
            tables/
            figures/
            reports/
        task6_obeta_special/
          raw/
            tables/
            figures/
            reports/
          qc/
            tables/
            figures/
            reports/
        task7_individual_checks/
          raw/
            tables/
            figures/
            reports/
          qc/
            tables/
            figures/
            reports/
      paper/
        raw/
          tables/
          figures/
        qc/
          tables/
          figures/
        report/
  latest_run.txt
```

---

## 4. Recommended semantic layers

## 4.1 `subjects/`

Purpose:
- single-subject outputs only
- one folder per `.set`

Recommended layout:

```text
subjects/
  sub01/
    tables/
      sub01_bandpower_roi.csv
      sub01_bandpower_summary.csv
      sub01_bandpower_tests.csv
      sub01_scene_level.csv
    figures/
      sub01_roi_occ.png
      sub01_scene_sequence.png
      ...
    qc/
      sub01_qc.csv
      sub01_pairs_check.csv
      sub01_marker_report.csv
    report/
      sub01_report.md
      input_set_path.txt
      config_used.json
```

Design rule:
- anything that belongs to exactly one subject stays here
- no merged tables here
- no task-level group analyses here

---

## 4.2 `batch/`

Purpose:
- cross-subject aggregation and QC filtering
- no inferential task packages here

Recommended layout:

```text
batch/
  merged/
    all_subjects_scene_level.csv
    all_subjects_pairs_check.csv
    per_subject_recovery_metrics.csv
    global_bandpower_summary.csv
  qc/
    all_subjects_scene_level_qc.csv
    all_subjects_pairs_check_qc.csv
    qc_exclusion_subjects.csv
    qc_scene_valid_counts.csv
  reports/
    summary_report.md
    qc_filter_report.md
    methods_snapshot.md
  audit/
    qc_report.json
```

Design rule:
- batch = merged data + QC + reproducibility docs
- batch should not contain task-specific inferential outputs

---

## 4.3 `analysis/`

Purpose:
- all statistical/interpretive task outputs
- task-oriented organization

Recommended layout pattern:

```text
analysis/
  taskX_name/
    raw/
      tables/
      figures/
      reports/
    qc/
      tables/
      figures/
      reports/
```

This is more readable than the current pattern:

```text
analysis-2/taskX/.../tables/<tag>/...
```

because it makes `raw` and `qc` visible earlier in the path.

### Example

```text
analysis/
  task4_core_lmm_suite/
    raw/
      tables/
      figures/
      reports/
    qc/
      tables/
      figures/
      reports/
```

Design rule:
- any output whose meaning depends on a specific analysis task belongs here
- all task-level CSV/PNG/MD outputs live under the corresponding task

---

## 4.4 `paper/`

Purpose:
- manuscript-facing or presentation-facing deliverables
- should contain only polished outputs, not exploratory intermediates

Recommended layout:

```text
paper/
  raw/
    tables/
    figures/
  qc/
    tables/
    figures/
  report/
    MASTER_REPORT.md
    figure_manifest.md
```

Design rule:
- paper = curated output layer
- should not be mixed with analytical staging tables

---

## 5. Recommended naming rules

A directory refactor alone will help, but naming consistency will make the structure much easier to maintain.

## 5.1 Keep stable names for core merged tables
Keep these well-known filenames for compatibility:

- `all_subjects_scene_level.csv`
- `all_subjects_pairs_check.csv`
- `all_subjects_scene_level_qc.csv`
- `all_subjects_pairs_check_qc.csv`
- `per_subject_recovery_metrics.csv`
- `qc_exclusion_subjects.csv`
- `qc_scene_valid_counts.csv`

These names are already used throughout the repository and downstream scripts.

## 5.2 Standardize figure names using semantic tokens
Recommended figure naming pattern:

```text
<scope>__<topic>__<metric>__<tag>.png
```

Examples:
- `group__complexity_effect__O_alpha__raw.png`
- `task3__trialindex_trend__O_theta__qc.png`
- `task7__audit_scorecard__summary__raw.png`

Advantages:
- easier sorting
- easier grep/filtering
- easier downstream manifests

## 5.3 Keep task-specific names where they already encode scientific meaning
For example, names like:
- `peakindex_summary_raw.csv`
- `block2_restart_stats_qc.csv`

are already meaningful and do not need aggressive renaming if directory structure becomes cleaner.

---

## 6. Recommended migration strategy

This refactor should **not** be implemented as a one-shot destructive change.

### Phase 1 — Path abstraction
Introduce a more complete output routing layer first.

Add helpers such as:
- `pipeline.get_run_dir(...)`
- `pipeline.get_subject_dir(...)`
- `pipeline.get_subject_table_dir(...)`
- `pipeline.get_subject_figure_dir(...)`
- `pipeline.get_subject_qc_dir(...)`
- `pipeline.get_subject_report_dir(...)`
- `pipeline.get_batch_dir(...)`
- `pipeline.get_batch_merged_dir(...)`
- `pipeline.get_batch_qc_dir(...)`
- `pipeline.get_batch_report_dir(...)`
- `pipeline.get_batch_audit_dir(...)`
- `pipeline.get_analysis_task_dir(...)`
- `pipeline.get_analysis_task_tag_dir(...)`
- `pipeline.get_paper_dir(...)`

At this stage:
- existing filenames can remain mostly unchanged
- existing logic can gradually switch to helper-based paths

### Phase 2 — Directory relocation with compatibility shims
Once all writers go through helper functions:
- redirect writes to the new layout
- optionally mirror or copy selected files to old locations during transition

### Phase 3 — Deprecate legacy layout
After downstream scripts are updated:
- remove old mirrored paths
- keep only the new structure

---

## 7. Minimal viable refactor scope

To reduce risk, the first refactor pass should focus on **directory structure**, not scientific logic.

### In-scope for first pass
1. Introduce output routing helpers
2. Move subject outputs into `subjects/<sid>/...`
3. Move merged/QC/methods reports into `batch/...`
4. Move task outputs into `analysis/taskX/{raw|qc}/...`
5. Move paper figures into `paper/{raw|qc}/...`

### Out-of-scope for first pass
1. Statistical model changes
2. Metric definition changes
3. Figure content redesign
4. Table schema redesign
5. Renaming every legacy file

This keeps the first pass operationally safe.

---

## 8. Files/functions that will need updates

### Highest priority
- `+pipeline/prepare_output.m`
- `+pipeline/get_output_root.m`
- `+pipeline/get_summary_dir.m`
- `+pipeline/get_fig_dir.m`
- `+pipeline/get_table_dir.m`
- `run_eeg_bandpower_pipeline.m`
- `summarize_bandpower_outputs.m`

### Batch/report helpers likely affected
- `+pipeline/qc_filter_and_report.m`
- `+pipeline/write_global_report_md.m`
- `+pipeline/write_methods_snapshot.m`
- `+pipeline/write_subject_report_md.m`

### Analysis task writers likely affected
- `+pipeline/analyze_block2_restart.m`
- `+pipeline/analyze_scene_block_diff.m`
- `+pipeline/analyze_trialindex_lmm.m`
- `+pipeline/analyze_core_lmm_suite.m`
- `+pipeline/analyze_peakindex_invertedu.m`
- `+pipeline/analyze_obeta_special.m`
- `+pipeline/analyze_individual_checks.m`

### Figure writers likely affected
- `+pipeline/plot_group_summaries.m`
- `+pipeline/plot_group_recovery_summaries.m`
- `+pipeline/plot_group_scene_sequences.m`
- `+pipeline/plot_group_scene_by_factors.m`
- `+pipeline/plot_group_topoplots.m`
- `+pipeline/plot_paper_figures.m`

---

## 9. Backward compatibility recommendations

To reduce breakage:

1. Keep core merged CSV filenames unchanged.
2. Keep core subject CSV filenames unchanged.
3. Change **directory layout first**, not table schema.
4. Add a short migration note in README and DEVELOPMENT docs.
5. Consider a temporary compatibility export layer for:
   - `paper_fig/`
   - `summary/` root-level reports
   - `analysis-2/...` legacy task paths

---

## 10. Proposed implementation order

### Step 1
Create new path helper functions.

### Step 2
Refactor subject output writing.

### Step 3
Refactor batch merged/QC/report output writing.

### Step 4
Refactor task1–task7 output roots.

### Step 5
Refactor paper figure output roots.

### Step 6
Add/update smoke tests for path layout.

### Step 7
Update README / DEVELOPMENT docs.

---

## 11. Suggested acceptance criteria

A successful output refactor should satisfy:

1. A new user can locate all single-subject outputs without reading code.
2. A new user can locate all merged/QC outputs without reading code.
3. Every analysis task has a self-contained output directory.
4. `raw` and `qc` are directory-level concepts for analysis outputs.
5. Paper outputs are isolated from intermediate analysis artifacts.
6. Existing core CSV schemas remain unchanged in the first pass.
7. Smoke tests still pass after refactor.

---

## 12. Recommendation

This refactor is worth doing **before** adding more analysis tasks or more figure families.

The repository is at the right stage for output-structure cleanup:
- complex enough that structure matters
- still small enough that migration is manageable
- already using helper functions, which makes path abstraction feasible

---

## 13. Short summary

Recommended direction:
- **subjects/** for per-subject outputs
- **batch/** for merged/QC/methods/audit
- **analysis/** for task-based results with `raw/` and `qc/`
- **paper/** for manuscript-facing outputs

Recommended strategy:
- **refactor path routing first**
- **preserve filenames and schemas in pass 1**
- **migrate gradually with compatibility support**
