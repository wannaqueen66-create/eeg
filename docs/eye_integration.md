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
- Design labels (`scene_name`, `WWR`, `Cond`, `Complexity`) should still follow the same long design table used by EEG.

## Recommended architecture

### Keep EEG segmentation unchanged
- `run_eeg_bandpower_pipeline.m` continues to handle EEG marker-state-machine segmentation and per-subject exports.
- Do **not** parse raw eye CSV inside the EEG segmentation entry script in v1.

### Add an eye scene-level branch
- Build a separate eye-tracking summary table from raw eye CSV files.
- One CSV is assumed to represent **one subject × one scene × one view segment**.
- Subject identity is parsed from the CSV filename.
- Base scene index is parsed from the parent folder.
- Design labels should later come from the same design long table used by EEG.

### Merge at batch summary stage
- Use `summarize_bandpower_outputs.m` to merge eye features into merged EEG scene-level tables.
- New merged output:
  - `all_subjects_scene_level_with_eye.csv`

## Priority metrics for assisting EEG

### Tier 1: artifact / QC assistance (highest priority)
These are the most useful eye-derived helpers for EEG cleaning interpretation and QC:
- `eye_tracking_ratio`
- `eye_valid_left_ratio`
- `eye_valid_right_ratio`
- `eye_view_blink_count`
- `eye_view_blink_rate_per_min`
- `eye_view_mean_blink_dur_ms`
- `eye_view_total_blink_dur_ms`
- `eye_view_blink_burden_pct`
- `eye_view_sacc_count`
- `eye_view_sacc_rate_per_min`
- `eye_view_mean_sacc_amp_px`
- `eye_view_mean_sacc_vel_px_ms`
- `eye_view_peak_sacc_vel_px_ms`
- `eye_mean_openness_pct`
- `eye_min_openness_left_pct`
- `eye_min_openness_right_pct`
- `eye_mean_eyelid_dist_mm`

### Tier 2: cognitive/state interpretation (secondary)
Useful for joint analysis and interpretation, but less direct for artifact judgment:
- `eye_mean_pupil_mm`
- `eye_view_fix_count`
- `eye_view_mean_fix_dur_ms`
- `eye_view_mean_gaze_velocity_px_ms`

## Output schema (eye scene-level, v1)
Typical columns:
- `subject_id`
- `scene_id`
- `block_id`
- `cycle_in_block`
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
- `eye_view_blink_count`
- `eye_view_blink_rate_per_min`
- `eye_view_mean_blink_dur_ms`
- `eye_view_total_blink_dur_ms`
- `eye_view_blink_burden_pct`
- `eye_view_sacc_count`
- `eye_view_sacc_rate_per_min`
- `eye_view_mean_sacc_dur_ms`
- `eye_view_mean_sacc_amp_px`
- `eye_view_mean_sacc_vel_px_ms`
- `eye_view_peak_sacc_vel_px_ms`
- `eye_mean_openness_pct`
- `eye_mean_eyelid_dist_mm`
- `eye_mean_pupil_mm`
- `eye_view_fix_count`
- `eye_view_fix_rate_per_min`
- `eye_view_mean_fix_dur_ms`

## Merge logic
Default merge keys:
- `subject_id`
- `scene_id`

Behavior:
- left join on EEG table
- preserve all EEG rows
- do not invent a second naming system for scenes; let design labels still come from the same design table path used by EEG

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

## Suggested uses in EEG analysis
- Use blink/saccade/openness/tracking metrics as **QC flags**.
- Use blink/saccade metrics as **artifact-interpretation covariates** when frontal EEG patterns look suspicious.
- Use pupil/fixation metrics as **secondary state/attention indicators** in multimodal analysis.
- Do **not** treat eye-tracking as direct EMG measurement; it is best used as indirect evidence for blink/eye-movement-related contamination.

## What v1 does NOT do yet
- event-level / trial-level eye segmentation from markers
- EEG-eye millisecond-level synchronization
- pair-level eye recovery table (`all_subjects_pairs_check_with_eye.csv`)

## Recommended next steps
1. Validate 3-5 real eye folders + filenames end-to-end.
2. Confirm that folder-based `scene_id = (block_id-1)*6 + cycle_in_block` matches EEG design mapping in all runs.
3. Add optional eye QC thresholds (e.g. minimum tracking ratio, maximum blink burden).
4. Add pair-level merged table in v2 if needed.
