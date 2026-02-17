# Changelog

## 2026-02-17
- feat: add preprocessing entry script `run_eeg_preprocess_pipeline.m` (bandpass, notch, reref, ICA, optional ICLabel reject).
- feat: add default preprocessing config `preprocess_config.json`.
- docs: update README/README.zh-CN with preprocessing usage.
- fix: correct low-gamma relative-power export to use high-band denominator path (`out_rel40`) instead of 1-30Hz path.
- fix: auto-expand high-band denominator upper bound when `bands.low_gamma` exceeds configured `totalBand40`.
- harden: add denominator-coverage guards in `band_power_v2` (return NaN when denominator does not cover target band).
- test: add smoke checks for low-gamma denominator mapping.
- test: add smoke check for scene_level schema fields and low_gamma support markers.
- docs: add statistical interpretation notes (gamma vs EMG artifact risk) to README and README.zh-CN.
- release: prepare reproducibility tag `v0.3.0`.
- feat: add low-gamma band support (`bands.low_gamma`, default `[30,45]`) and export ROI low-gamma columns (`F/P/O_low_gamma`).
- feat: enrich `scene_level.csv` with `block_id`, `cycle_in_block`, `pair_id`, low/high beta + low-gamma columns, and TAR/TBR/BA ratio metrics.
- docs: rewrite `README.md` as bilingual (EN+ZH) with TOC; add dedicated Chinese doc `README.zh-CN.md`.
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
