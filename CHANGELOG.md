# Changelog

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
