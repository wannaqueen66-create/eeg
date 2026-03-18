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
- `pipeline.prepare_output` directory creation + config snapshot
- staged output path helpers for:
  - subject dirs
  - batch dirs
- static schema checks for key scene-level exports
- low-gamma denominator mapping checks
- curated main node report generation (`report/README_NODE.md`) when curated output helpers are present
- curated main builder hook checks for all 4 main nodes (structure-level check only; curated outputs are disabled by default in current `main` unless explicitly enabled)

## Notes

- These are smoke tests (fast checks), not full numerical validation tests.
- Full regression tests with real `.set` datasets should be added later.
