# Eye-tracking QC Recommendations for Assisting EEG

This note is for the current EEG × eye-tracking workflow, where:
- EEG is segmented by marker/state machine from `.set`
- eye-tracking CSV files already represent scene-level `view` segments
- eye features are merged at the batch-summary stage

## Purpose
Eye-tracking is **not** used here as direct EMG measurement.
Instead, it is used to support EEG quality control and artifact interpretation, especially for:
- blink-related contamination
- eye-movement-related contamination
- low-confidence scene segments due to poor eye-tracking quality

## Recommended priority metrics

### Tier 1 — best for EEG QC / artifact support
- `eye_tracking_ratio`
- `eye_valid_left_ratio`
- `eye_valid_right_ratio`
- `eye_view_blink_count`
- `eye_view_blink_rate_per_min`
- `eye_view_blink_burden_pct`
- `eye_view_mean_blink_dur_ms`
- `eye_view_sacc_count`
- `eye_view_sacc_rate_per_min`
- `eye_view_mean_sacc_amp_px`
- `eye_view_mean_sacc_vel_px_ms`
- `eye_view_peak_sacc_vel_px_ms`
- `eye_mean_openness_pct`
- `eye_min_openness_left_pct`
- `eye_min_openness_right_pct`
- `eye_mean_eyelid_dist_mm`

### Tier 2 — useful for state/attention interpretation
- `eye_mean_pupil_mm`
- `eye_view_fix_count`
- `eye_view_fix_rate_per_min`
- `eye_view_mean_fix_dur_ms`
- `eye_view_mean_gaze_velocity_px_ms`

## Practical interpretation

### 1) Blink-heavy scenes
Potentially more likely to contain frontal ocular contamination when:
- `eye_view_blink_count` is high
- `eye_view_blink_rate_per_min` is high
- `eye_view_blink_burden_pct` is high
- `eye_min_openness_*` is low

Use case:
- explain suspicious frontal EEG fluctuations
- flag scenes for review before over-interpreting frontal power changes

### 2) Eye-movement-heavy scenes
Potentially more likely to contain eye-movement-related contamination when:
- `eye_view_sacc_count` is high
- `eye_view_sacc_rate_per_min` is high
- `eye_view_mean_sacc_amp_px` is high
- `eye_view_peak_sacc_vel_px_ms` is high

Use case:
- interpret scenes with strong visual search / scanning behavior
- avoid over-attributing frontal changes to cognitive factors alone

### 3) Low-confidence eye segments
Low-confidence scenes when:
- `eye_tracking_ratio` is low
- `eye_valid_left_ratio` / `eye_valid_right_ratio` are low
- openness / eyelid values look implausible or missing

Use case:
- flag eye-derived features as low confidence
- avoid strong conclusions from pupil/fixation metrics in these scenes

## Suggested exploratory flags (do NOT treat as hard clinical cutoffs)
These are only starting points for internal QC review and should be tuned on your own data distribution.

### Eye quality flags
- `low_tracking_flag`: `eye_tracking_ratio < 60`
- `low_validity_flag`: left or right validity ratio < `0.6`

### Blink burden flags
- `high_blink_rate_flag`: `eye_view_blink_rate_per_min` above cohort 75th/90th percentile
- `high_blink_burden_flag`: `eye_view_blink_burden_pct` above cohort 75th/90th percentile

### Saccade activity flags
- `high_sacc_rate_flag`: `eye_view_sacc_rate_per_min` above cohort 75th/90th percentile
- `high_sacc_velocity_flag`: `eye_view_peak_sacc_vel_px_ms` above cohort 90th percentile

## Recommended first-pass workflow
1. Generate `all_subjects_eye_scene_level.csv`
2. Merge into `all_subjects_scene_level_with_eye.csv`
3. Plot distributions of:
   - blink rate
   - blink burden
   - tracking ratio
   - saccade rate
4. Flag top/bottom outliers
5. Re-check suspicious EEG scenes with these flags in mind
6. In models, add only a **small number** of eye covariates first

## Suggested first covariates for EEG models
If you want a minimal covariate set, start with:
- `eye_tracking_ratio`
- `eye_view_blink_rate_per_min`
- `eye_view_sacc_rate_per_min`

If you want a slightly richer set:
- `eye_tracking_ratio`
- `eye_view_blink_burden_pct`
- `eye_view_peak_sacc_vel_px_ms`
- `eye_mean_pupil_mm`

## What not to do first
- do not put every eye variable into one model
- do not interpret pupil as artifact by default
- do not treat eye-tracking as direct EMG measurement
- do not use hard thresholds before checking your own cohort distributions

## Recommended wording in reports
- "Eye-tracking-derived blink and saccade metrics were used as auxiliary QC / artifact-interpretation indicators for EEG, not as direct EMG measures."
- "Eye-tracking quality indicators (tracking ratio / validity) were used to assess confidence in multimodal scene-level interpretation."
