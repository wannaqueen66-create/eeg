# EEG × Eye-tracking Integration Plan

## Goal
Integrate eye-tracking features into the existing EEG batch-summary workflow **without destabilizing** the current EEG `.set` processing pipeline.

## Current confirmed experiment assumptions
- Eye-tracking CSV files are **already manually cut and confirmed as scene `view` segments**.
- Eye-tracking side does **not** need marker-based segmentation in v1.
- EEG side still uses marker/state-machine segmentation from `.set` files.
- EEG subject identifiers are **Chinese participant names** in `.set` filenames.
- Eye CSV file names also contain Chinese participant names.
- Scene identity is stored primarily in the **parent folder name**, such as:
  - `1-1-1_组1-C1W45`
  - interpreted as `order_id=1`, `round_id=1`, `scene_index_in_round=1`

## Recommended architecture

### Keep EEG segmentation unchanged
- `run_eeg_bandpower_pipeline.m` continues to handle EEG marker-state-machine segmentation and per-subject exports.
- Do **not** parse raw eye CSV inside the EEG segmentation entry script in v1.

### Add an eye scene-level branch
- Build a separate eye-tracking summary table from raw eye CSV files.
- One CSV is assumed to represent **one subject × one scene × one view segment**.
- Scene identity is parsed from the parent folder.
- Subject identity is parsed from the CSV filename.
- Output an analysis-friendly table keyed by:
  - `subject_id`
  - `scene_id`

### Merge at batch summary stage
- Use `summarize_bandpower_outputs.m` to merge eye features into merged EEG scene-level tables.
- New merged output:
  - `all_subjects_scene_level_with_eye.csv`

## New files added in v1
- `scripts/build_eye_scene_level.py`
  - builds eye scene-level table from raw eye CSV files
- `+pipeline/merge_eye_scene_features.m`
  - merges eye scene-level table into EEG scene-level table

## Expected eye input columns (v1 priority)

### Quality
- `Tracking Ratio[%]`
- `Validity Left`
- `Validity Right`

### Core features
- `Pupil Diameter Left[mm]`
- `Pupil Diameter Right[mm]`
- `Gaze Velocity[px/ms]`
- `Fixation Index`
- `Fixation Duration[ms]`
- `Saccade Index`
- `Saccade Duration[ms]`
- `Saccade Amplitude[px]`
- `Saccade Velocity Average[px/ms]`
- `Saccade Velocity Peak[px/ms]`
- `Blink Index`
- `Blink Duration[ms]`

### Optional metadata retained
- `User`
- `Record Name`
- `性别`
- `右眼镜片度数`
- `左眼镜片度数`

## Output schema (eye scene-level, v1)
Typical columns:
- `subject_id`
- `subject_name`
- `scene_id`
- `order_id`
- `round_id`
- `scene_index_in_round`
- `scene_folder`
- `scene_label`
- `WWR`
- `Complexity`
- `eye_source_file`
- `eye_n_rows`
- `eye_view_start_ms`
- `eye_view_end_ms`
- `eye_view_duration_ms`
- `eye_tracking_ratio`
- `eye_valid_left_ratio`
- `eye_valid_right_ratio`
- `eye_mean_pupil_left_mm`
- `eye_mean_pupil_right_mm`
- `eye_mean_pupil_mm`
- `eye_view_mean_gaze_velocity_px_ms`
- `eye_view_fix_count`
- `eye_view_mean_fix_dur_ms`
- `eye_view_sacc_count`
- `eye_view_mean_sacc_dur_ms`
- `eye_view_mean_sacc_amp_px`
- `eye_view_mean_sacc_vel_px_ms`
- `eye_view_peak_sacc_vel_px_ms`
- `eye_view_blink_count`
- `eye_view_mean_blink_dur_ms`

## Merge logic
Default merge keys:
- `subject_id`
- `scene_id`

Behavior:
- left join on EEG table
- preserve all EEG rows
- auto-prefix eye columns with `eye_` when needed to avoid name collisions

## Config options
Add to `config.json`:
```json
{
  "eye_merge_enabled": false,
  "eye_summary_path": "",
  "eye_merge_keys": ["subject_id", "scene_id"]
}
```

## Example workflow

### 1) Build eye scene-level table
```bash
python scripts/build_eye_scene_level.py \
  --input /path/to/eye_root \
  --output /path/to/eye_outputs/summary/all_subjects_eye_scene_level.csv
```

### 2) Enable merge in `config.json`
```json
{
  "eye_merge_enabled": true,
  "eye_summary_path": "/path/to/eye_outputs/summary/all_subjects_eye_scene_level.csv"
}
```

### 3) Run EEG batch summary
```matlab
summarize_bandpower_outputs('path/to/eeg_folder', 'config.json')
```

### 4) Read merged output
- `all_subjects_scene_level.csv`
- `all_subjects_scene_level_with_eye.csv`

## What v1 does NOT do yet
- event-level / trial-level eye segmentation from markers
- EEG-eye millisecond-level synchronization
- pair-level eye recovery table (`all_subjects_pairs_check_with_eye.csv`)

## Recommended next steps
1. Validate 3-5 real eye folders + filenames end-to-end.
2. Confirm that folder-based `scene_id = (round_id-1)*6 + scene_index_in_round` matches EEG design mapping in all runs.
3. Add optional QC thresholds for eye-tracking quality (e.g. minimum tracking ratio).
4. Add pair-level merged table in v2 if needed.
