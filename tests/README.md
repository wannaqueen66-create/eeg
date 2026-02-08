# Tests (Lightweight)

This folder contains lightweight smoke tests for the modularized entry pipeline.

## Run

In MATLAB (repo root):

```matlab
addpath(genpath(pwd));
run('tests/run_smoke_tests.m');
```

## Coverage (current)

- `pipeline.load_config` default behavior
- `pipeline.parse_input` extension validation
- `pipeline.prepare_output` directory creation + `config_used.json`

## Notes

- These are smoke tests (fast checks), not full numerical validation tests.
- Full regression tests with real `.set` datasets should be added later.
