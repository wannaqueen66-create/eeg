# DEVELOPMENT

## Architecture (current)

- `run_eeg_bandpower_pipeline.m`
  - orchestration entry (batch loop + stage-aware error context)
- `+pipeline/parse_input.m`
  - input resolution (file/folder/GUI)
- `+pipeline/load_config.m`
  - config load + defaults
- `+pipeline/prepare_output.m`
  - output tree setup + config snapshot
- `+pipeline/resolve_output_dir.m`
  - path resolver
- `+pipeline/write_config_snapshot.m`
  - persistence helper

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
