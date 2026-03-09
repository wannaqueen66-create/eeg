# New Main Redesign Plan

> Goal: freeze the current repository state into a `raw` branch, then rebuild a cleaner, more focused `main` branch without radically changing the core MATLAB EEG processing logic.

---

# 1. Why redesign the current `main`?

The current repository already contains substantial functionality, but the structure has become too broad and too mixed:

- core EEG pipeline logic
- output migration work-in-progress
- multiple analysis tasks
- descriptive figures
- inferential models
- QC audit logic
- paper-facing exports
- historical/compatibility layers

This makes the repository powerful, but also:
- harder to understand quickly
- harder to maintain
- harder to present as a clean research pipeline
- harder to align with a clear paper story

So the redesign should not start from zero, but should **refocus** the project.

---

# 2. High-level strategy

## Step A — Freeze current state into `raw`
Create a dedicated `raw` branch that preserves the current repository state.

### Purpose of `raw`
- archival / historical reference
- preserve all current outputs, task structure, and compatibility layers
- allow later lookup of older logic without cluttering the new `main`

`raw` should be treated as:
- a stable archive branch
- not the primary branch for ongoing cleanup

---

## Step B — Rebuild a cleaner `main`
The new `main` should:
- keep the core MATLAB EEG logic broadly intact
- reduce output and analysis sprawl
- prioritize clarity over completeness
- present results through two main categories:
  - **descriptive analysis**
  - **inferential analysis**

---

# 3. Core design principle for the new `main`

The new `main` should answer:

1. What are the descriptive EEG patterns?
2. What are the statistically supported effects?

This suggests a top-level split between:
- `descriptive`
- `inferential`

rather than exposing the repository mainly through `task1` / `task2` / `task3` / ...

The task-oriented structure can still exist internally in code if useful, but should no longer dominate the external output organization.

---

# 4. What should the new `main` keep?

## Keep the core EEG processing logic
The following should remain conceptually intact:

- marker-state-machine segmentation
- ROI-based EEG bandpower extraction
- QC metrics
- scene-level export
- pair / recovery export
- batch-level merged tables
- design table attachment

In other words, do **not** rebuild the science pipeline from scratch.

## Keep the core four metrics as the default focus
Retain the current core metric set:
- `O_alpha`
- `O_theta`
- `O_beta`
- `F_theta`

These are already the most coherent repository-wide core set.

---

# 5. What should the new `main` simplify?

## 5.1 Simplify output categories
Instead of exposing many task folders first, outputs should be framed as:

- **single-subject outputs**
- **descriptive analysis outputs**
- **inferential analysis outputs**
- **QC / logs / audit**

## 5.2 Simplify default grouping strategy
For the first redesigned `main`, only keep:

- `overall`
- `experience`

as primary output categories.

This means:
- no need to remove `SportFreq` logic entirely
- but remove it from the default main output surface
- keep it optional, secondary, or archived in `raw`

### Current implementation direction
The redesigned main branch should prefer:
- direct summary tables for `overall`
- direct grouped summary tables for `experience`
- direct inferential summary tables and PNGs for `overall`
- direct inferential summary tables and PNGs for `experience`
- mirrored detailed task files only as a secondary support layer

while keeping `SportFreq` accessible only through detailed task outputs unless explicitly re-promoted later.

## 5.3 Simplify the visible research narrative
The new `main` should focus on:

### Descriptive
- overall EEG patterns
- overall recovery patterns
- experience-group descriptive comparisons

### Inferential
- condition effects
- experience-group modulation effects
- a smaller, cleaner set of confirmatory tests

---

# 6. Proposed new output structure for redesigned `main`

```text
bandpower_outputs/
  subjects/
    <subject_id>/
      tables/
      figures/
      qc/
      report/

  descriptive/
    overall/
      tables/
      figures/
      report/
    experience/
      tables/
      figures/
      report/

  inferential/
    overall/
      tables/
      figures/
      report/
    experience/
      tables/
      figures/
      report/

  qc/
  logs/
  audit/
```

This is intentionally more compact than the current staged task-heavy layout.

---

# 7. Recommended analysis split

## 7.1 Descriptive analysis

### A. Overall
Purpose:
- provide a clean descriptive summary of the full sample

Recommended content:
- overall condition means (WWR / Complexity)
- overall scene-level summaries
- overall recovery summaries
- overall sequence plots (if still retained)
- overall QC summaries

### B. Experience group
Purpose:
- compare `Experience High` vs `Experience Low`

Recommended content:
- group mean plots by Complexity
- group mean plots by WWR
- recovery descriptive plots by group
- selected scene sequence plots by group

### Explicitly not required in first-pass descriptive layer
- `SportFreq` default descriptive outputs
- a large number of exploratory task-specific plots

---

## 7.2 Inferential analysis

### A. Overall inferential
Purpose:
- answer the main condition-effect questions in the full sample

Recommended content:
- WWR main effect
- Complexity main effect
- WWR × Complexity interaction
- optional WWR trend / quadratic check
- optional PeakIndex summary if retained

### B. Experience inferential
Purpose:
- test whether experience modifies the main effects

Recommended content:
- Group main effect
- WWR × Group
- Complexity × Group
- WWR × Complexity × Group
- selected robustness checks

### Optional / secondary inferential items
These can be kept in code, but not emphasized in the new `main` output by default:
- SportFreq group inferential outputs
- highly exploratory scene-specific tests
- broader task museum-style exports

---

# 8. Which current task families should be retained, merged, or downgraded?

## Retain in concept
### Current Task4
Keep as the core inferential backbone.

### Current Task3
Retain, but frame as a control / sequence-check component, not as a top-level public output category.

### Current recovery outputs
Retain, especially the `view -> gray` descriptive and inferential summaries.

## Keep but downgrade in visibility
### Current Task1
Block2 restart can remain as a control analysis, but should not dominate the visible output tree.

### Current Task2
Specific scene block-difference analysis can remain optional or supplemental.

### Current Task5
PeakIndex / inverted-U can remain, but should be framed as a focused secondary analysis rather than always occupying a top-level output slot.

### Current Task6
Core-metric special robustness can remain as a secondary robustness layer.

### Current Task7
Individual outlier / influence audit should remain, but be framed as audit/support rather than as a central analysis product.

---

# 9. Recommended output priorities in the new `main`

## Highest priority outputs
These should be easiest to find and best maintained:

### Subject level
- scene-level table
- pair/recovery table
- qc table
- subject report

### Descriptive
- overall descriptive tables/plots
- experience descriptive tables/plots

### Inferential
- overall inferential summary
- experience inferential summary

## Lower priority / secondary outputs
- detailed task-specific sub-analyses
- sportfreq-first outputs
- highly exploratory paper exports
- compatibility output mirrors

---

# 10. Suggested repository structure changes

## Code organization suggestion
Keep core processing entry points but reorganize analysis helpers around the new split:

### Processing layer
- input / config / loading
- segmentation
- power extraction
- QC
- single-subject export

### Descriptive layer
- overall summaries
- experience summaries
- recovery descriptive summaries

### Inferential layer
- overall inferential models
- experience inferential models
- secondary robustness models

### Audit layer
- QC reports
- outlier/influence checks

This does not require deleting current task code immediately.
It may start as a re-orchestration problem rather than a full rewrite.

---

# 11. What should be removed from the first redesigned `main` surface?

Not necessarily deleted from the repository, but removed from the default visible top-level output flow:

- `SportFreq` as a default primary result branch
- task-by-task public output sprawl
- too many parallel figure families
- over-detailed paper-output layers before descriptive/inferential layers are stabilized

---

# 12. Suggested implementation order

## Phase 1 — Branch management
1. Create `raw` from current `main`
2. Confirm `raw` is pushed and frozen as archival branch

## Phase 2 — New `main` output architecture
1. Simplify output root structure
2. Replace task-dominant external layout with:
   - `subjects`
   - `descriptive`
   - `inferential`
   - `qc/logs/audit`

## Phase 3 — Re-orchestrate analysis writing
1. Route current overall summaries into `descriptive/overall`
2. Route current experience summaries into `descriptive/experience`
3. Route selected models into `inferential/overall`
4. Route selected experience models into `inferential/experience`

### Status note
A first implementation pass may begin as a curated mirror layer built on top of the existing task outputs,
rather than a full immediate rewrite of all internal task logic.
This is acceptable as a transition step.

## Phase 4 — Reduce visible surface area
1. De-emphasize SportFreq in default outputs
2. Keep advanced task outputs optional or internal
3. Update README and output guide

## Phase 5 — Validation
1. Re-run with design table
2. Verify all major descriptive outputs exist
3. Verify inferential outputs exist and match intended scope
4. Confirm no broken references remain

---

# 13. Acceptance criteria for the redesigned `main`

The redesign is successful if:

1. A new user can understand the output structure without reading internal analysis code.
2. The top-level outputs clearly separate descriptive vs inferential results.
3. `overall` and `experience` become the default visible result categories.
4. The repository no longer feels like a collection of many partially overlapping task export systems.
5. The core MATLAB EEG logic is preserved.
6. The `raw` branch safely preserves the current full repository state.

---

# 14. Suggested final positioning

## `raw` branch
- full historical repository
- complex but complete
- archival / compatibility / reference

## redesigned `main`
- cleaner research-facing branch
- focused on the primary story
- simpler outputs
- easier to maintain and explain
- prefer curated top-level reading flow:
  - `descriptive/overall`
  - `descriptive/experience`
  - `inferential/overall`
  - `inferential/experience`

---

# 15. Short summary

Recommended next move:

- freeze the current repository as `raw`
- rebuild `main` around two analysis categories:
  - **descriptive**
  - **inferential**
- keep only two default visible result branches:
  - **overall**
  - **experience**
- preserve the core EEG MATLAB logic, but simplify the repository surface
- enforce a consistent visible folder rule:
  - `tables/` for CSV only
  - `figures/` for PNG only
  - `report/` for README/markdown only

This approach avoids throwing away the current work while making the future `main` far more coherent.
