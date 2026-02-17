# Changelog

## 2026-02-17
- fix: remove duplicated in-script config-default block; defaults now come from `pipeline.load_config` only.
- fix: correct view-gray pairing in `export_summary_tables` by using table-row indices (via shared helper), avoiding condition-index mismatch.
- fix: wire `summary_files` collection during batch run so `global_summary=true` can generate `global_bandpower_summary.csv`.
- feat: enable optional zip packaging when `zip_output=true` by calling `zip_output_files` at end of each dataset run.
- refactor: extract shared `find_view_gray_pairs` helper and reuse in pair-id assignment / paired plots / pairs-check / summary tests.
- refactor: remove unused local helper functions (`get_events_sorted`, `load_cfg`, local `resolve_output_dir`, `validate_cycle_counts`, local `write_config_snapshot`).
- improve: unify topoplot Welch parameters with pipeline-level `wlen/nover/nfft` for consistency.

## 2026-02-08
- fix: define `base` before first use in `run_eeg_bandpower_pipeline.m`.
- refactor: remove redundant late `fileparts(fn)` reassignment for `base`.
- note: no intended behavior/output changes.
- refactor: introduce `+pipeline` modules (`parse_input`, `load_config`, `prepare_output`, `resolve_output_dir`, `write_config_snapshot`) and route entry script to use them.
- note: intended to preserve behavior/output while reducing entry-script responsibilities.
- refactor: add per-file `try/catch` orchestration in entry pipeline with stage labels (`load_eeg`, `setup_bands_roi`, `segment_state_machine`, `compute_bandpower_qc`, `build_tables_export_csv`, `plot_figures`).
- improve: error logs now include dataset path + failing stage context to speed debugging in batch runs.
- test: add lightweight smoke tests in `tests/run_smoke_tests.m` for config/input/output helpers.
- docs: add `DEVELOPMENT.md`, README architecture section, and test instructions.
