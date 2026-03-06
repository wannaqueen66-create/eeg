# EEG × Eye Minimal Modeling Template

This document gives a practical first-pass modeling strategy after generating:
- `all_subjects_scene_level.csv`
- `all_subjects_scene_level_with_eye.csv`

## Goal
Use a **small** set of eye-tracking variables to:
1. support EEG QC / artifact interpretation
2. test whether main EEG effects remain stable after accounting for eye behavior

Do not start with a giant all-variable model.

## Recommended role split

### A. QC / review variables
Best used for flagging, diagnostics, and cautious interpretation:
- `eye_qc_needs_review`
- `eye_qc_severe_flag`
- `eye_qc_flag_count`
- `eye_tracking_ratio`
- `eye_view_blink_burden_pct`

### B. Minimal artifact-support covariates
Best first covariates to add to EEG models:
- `eye_tracking_ratio`
- `eye_view_blink_rate_per_min`
- `eye_view_sacc_rate_per_min`

### C. Secondary state/attention covariates
Add only if needed after the minimal model works:
- `eye_mean_pupil_mm`
- `eye_view_fix_rate_per_min`
- `eye_view_mean_fix_dur_ms`

## Suggested modeling ladder

### Model 0 — EEG only
Example:
- EEG metric ~ WWR * Complexity + (1 | subject_id)

Use this as your baseline reference model.

### Model 1 — add minimal eye QC covariates
Example:
- EEG metric ~ WWR * Complexity + eye_tracking_ratio + eye_view_blink_rate_per_min + eye_view_sacc_rate_per_min + (1 | subject_id)

Purpose:
- check whether the main EEG effects remain after controlling for basic eye-behavior variation

### Model 2 — add one state variable (optional)
Example:
- EEG metric ~ WWR * Complexity + eye_tracking_ratio + eye_view_blink_rate_per_min + eye_view_sacc_rate_per_min + eye_mean_pupil_mm + (1 | subject_id)

Purpose:
- test whether pupil-linked arousal/load explains additional variance

## What to compare across models
- coefficient direction and magnitude for WWR / Complexity
- p-values / confidence intervals of main effects
- whether eye covariates absorb suspicious frontal effects
- whether model fit improves meaningfully

## When to use QC flags directly
Use QC flags mainly for:
- sensitivity analysis
- subgroup review
- excluding clearly problematic scenes in a supplementary robustness check

Example sensitivity subset:
- keep only rows where `eye_qc_severe_flag == 0`

Do not make this your only analysis.

## Good first EEG outcomes to pair with eye covariates
Most natural candidates:
- frontal theta
- occipital alpha
- occipital beta
- recovery delta metrics if later merged at pair-level

## Cautions
- blink and saccade covariates can be correlated
- fixation and saccade variables can be redundant
- pupil may reflect cognition/state, not artifact
- tracking ratio is partly a data-quality variable, not a pure behavior variable

## Recommended first pass for the current project
If you want the simplest workable version, start with:
- EEG metric ~ WWR * Complexity + eye_tracking_ratio + eye_view_blink_rate_per_min + eye_view_sacc_rate_per_min + (1 | subject_id)

Then compare against:
- EEG metric ~ WWR * Complexity + (1 | subject_id)

## Recommended wording
- "We used eye-tracking-derived blink and saccade measures as auxiliary covariates to test whether EEG effects were robust to eye-behavior variation."
- "Eye-tracking QC flags were used for sensitivity analysis and cautious interpretation, not as sole hard exclusion rules."
