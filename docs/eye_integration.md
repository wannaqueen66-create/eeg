# EEG × Eye-tracking Integration Plan

## Goal
Integrate eye-tracking features into the existing EEG batch-summary workflow **without destabilizing** the current EEG `.set` processing pipeline.

## Recommended architecture

### Keep EEG segmentation unchanged
- `run_eeg_bandpower_pipeline.m` continues to handle EEG marker-state-machine segmentation and per-subject exports.
- Do **not** parse raw eye CSV inside the EEG segmentation entry script in v1.

### Add an eye scene-level branch
- Build a separate eye-tracking summary table from raw eye CSV files.
- One CSV is assumed to represent **one subject × one scene**.
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

### Synchronization / metadata
- `Recording Time Stamp[ms]`
- `Triggle Send` / `Triggle Receive`
- `Event Label`
- `Annotation`

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

## Output schema (eye scene-level, v1)
Typical columns:
- `subject_id`
- `scene_id`
- `eye_source_file`
- `eye_n_rows`
- `eye_event_start_ms`
- `eye_event_end_ms`
- `eye_duration_ms`
- `eye_tracking_ratio`
- `eye_valid_left_ratio`
- `eye_valid_right_ratio`
- `eye_mean_pupil_left_mm`
- `eye_mean_pupil_right_mm`
- `eye_mean_pupil_mm`
- `eye_mean_gaze_velocity_px_ms`
- `eye_fix_count`
- `eye_mean_fix_dur_ms`
- `eye_sacc_count`
- `eye_mean_sacc_dur_ms`
- `eye_mean_sacc_amp_px`
- `eye_mean_sacc_vel_px_ms`
- `eye_peak_sacc_vel_px_ms`
- `eye_blink_count`
- `eye_mean_blink_dur_ms`

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
  --input /path/to/eye_csvs \
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
1. Standardize eye CSV naming so `subject_id` and `scene_id` can be parsed reliably.
2. Confirm how scene IDs map across EEG / eye-tracking / experiment design table.
3. Add optional marker-window extraction for eye `view` vs `gray` segments.
4. Add pair-level merged table in v2.
