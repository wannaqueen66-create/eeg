# EEG × Eye Quickstart

This is the fastest practical path after enabling eye merge.

## Inputs
You should already have:
- `all_subjects_eye_scene_level.csv`
- `all_subjects_scene_level_with_eye.csv`
- `design_path` configured for EEG scene labels

## What the new exploratory analysis does
It screens basic scene-level links between EEG and eye metrics using:
- correlation table
- minimal LMMs (when `fitlme` is available)
- a short markdown report

## Default EEG metrics
- `F_theta`
- `O_alpha`
- `O_beta`

## Default eye metrics
- `eye_tracking_ratio`
- `eye_view_blink_rate_per_min`
- `eye_view_sacc_rate_per_min`
- `eye_mean_pupil_mm`

## Suggested first reading order
1. `eye_qc_report.md`
2. `eye_qc_scene_summary.csv`
3. `analysis-eye/raw/reports/eye_eeg_report.md`
4. `analysis-eye/raw/tables/eye_eeg_correlation_screen.csv`

## What to look for first
- Are frontal EEG metrics associated with blink/saccade variables?
- Does tracking ratio look like a major confound?
- Do main scene effects still make sense after checking eye QC flags?
- Is pupil linked to EEG in a plausible state-related way rather than as artifact?

## Good first sanity checks
- inspect rows where `eye_qc_severe_flag == 1`
- compare model results with and without eye covariates
- avoid over-interpreting isolated correlations without checking QC context
