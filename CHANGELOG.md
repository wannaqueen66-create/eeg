# Changelog

## 2026-02-08
- fix: define `base` before first use in `run_eeg_bandpower_pipeline.m`.
- refactor: remove redundant late `fileparts(fn)` reassignment for `base`.
- note: no intended behavior/output changes.
- refactor: introduce `+pipeline` modules (`parse_input`, `load_config`, `prepare_output`, `resolve_output_dir`, `write_config_snapshot`) and route entry script to use them.
- note: intended to preserve behavior/output while reducing entry-script responsibilities.
