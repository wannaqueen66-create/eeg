# DEVELOPMENT

## Architecture (current)

- `run_eeg_bandpower_pipeline.m`
  - orchestration entry (batch loop + stage-aware error context)
- `+pipeline/parse_input.m`
  - input resolution (file/folder/GUI)
- `+pipeline/load_config.m`
  - config load + defaults
- `+pipeline/prepare_output.m`
  - subject output tree setup + config snapshot
- `+pipeline/resolve_output_dir.m`
  - path resolver
- `+pipeline/write_config_snapshot.m`
  - persistence helper

## Output structure refactor (staged)

The repository is migrating toward clearer output semantics:

- `subjects/` for per-subject outputs
- `batch/` for merged/QC/report/audit outputs
- `analysis/` for task-oriented inferential outputs
- `paper/` for manuscript-facing deliverables

Current implementation uses a staged migration strategy:
- add path helpers first
- migrate writers gradually
- treat `batch/` as the canonical batch-level surface in `main`
- keep `summary/` only as a compatibility layer where practical

Key helpers introduced for this migration:
- `pipeline.get_run_dir`
- `pipeline.get_subject_dir`
- `pipeline.get_subject_table_dir`
- `pipeline.get_subject_figure_dir`
- `pipeline.get_subject_qc_dir`
- `pipeline.get_subject_report_dir`
- `pipeline.get_batch_dir`
- `pipeline.get_batch_merged_dir`
- `pipeline.get_batch_qc_dir`
- `pipeline.get_batch_report_dir`
- `pipeline.get_batch_audit_dir`
- `pipeline.get_summary_dir` (compatibility alias only; prefer `get_batch_dir` in new code)
- `pipeline.get_analysis_task_subdirs`
- `pipeline.get_paper_subdirs`
- `pipeline.list_subject_ids`

## Contribution guidelines

1. Keep entry function signature stable:
   - `run_eeg_bandpower_pipeline(input_path, config_path)`
2. Prefer adding small focused helpers under `+pipeline/`.
3. Preserve output compatibility where possible:
   - existing CSV names
   - key column names used downstream
4. Add/update smoke tests in `tests/` when changing input/config/output behavior.

## Quick check before push

```matlab
addpath(genpath(pwd));
run('tests/run_smoke_tests.m');
```
