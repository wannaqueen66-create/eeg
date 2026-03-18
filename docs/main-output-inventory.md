# Main Output Inventory / main 输出盘点

## Table of Contents / 目录
- [1. Purpose / 用途](#1-purpose--用途)
- [2. Canonical Batch Root / 规范 batch 根目录](#2-canonical-batch-root--规范-batch-根目录)
- [3. Output Role Classes / 输出角色分类](#3-output-role-classes--输出角色分类)
- [4. Current Directory Inventory / 当前目录盘点](#4-current-directory-inventory--当前目录盘点)
- [5. Duplication Risk Map / 重复风险图谱](#5-duplication-risk-map--重复风险图谱)
- [6. Safe Preservation Rules / 安全保留规则](#6-safe-preservation-rules--安全保留规则)
- [7. Next Cleanup Targets / 下一步清理目标](#7-next-cleanup-targets--下一步清理目标)

## 1. Purpose / 用途

This document inventories the current `main`-branch output structure and assigns each directory a role:
- `source`
- `curated-summary`
- `mirror-copy`
- `navigation-only`

The goal is to enable de-duplication **without losing any PNG figure or related data**.

本文件对当前 `main` 分支输出结构进行盘点，并给每个目录标记角色：
- `source`（源结果层）
- `curated-summary`（摘要层）
- `mirror-copy`（镜像副本层）
- `navigation-only`（导航层）

目的是为后续去重提供依据，同时确保 **PNG 图和相关数据完全保留**。

---

## 2. Canonical Batch Root / 规范 batch 根目录

Current main-branch canonical batch root:

```text
<input_folder>/bandpower_outputs/runs/current/batch/
```

This is created/used by:
- `pipeline.get_batch_dir`
- `summarize_bandpower_outputs`

Single-subject canonical root:

```text
<input_folder>/bandpower_outputs/runs/current/subjects/<subject_id>/
```

---

## 3. Output Role Classes / 输出角色分类

### `source`
Original result-bearing outputs. These are canonical and must be preserved.

### `curated-summary`
Summary-level outputs that intentionally reorganize source results for easier reading.

### `mirror-copy`
Convenience copies of existing outputs. Useful for navigation, but high-risk for duplication and version drift.

### `navigation-only`
README / index / file map / task count style guidance files. These provide orientation but are not scientific source outputs.

---

## 4. Current Directory Inventory / 当前目录盘点

## 4.1 Subject-level / 被试级

### `subjects/<subject_id>/tables/`
- Role: `source`
- Contains: subject-level source CSV outputs
- Keep: yes

### `subjects/<subject_id>/figures/`
- Role: `source`
- Contains: subject-level source PNG outputs
- Keep: yes

### `subjects/<subject_id>/qc/`
- Role: `source`
- Contains: subject QC/source support files
- Keep: yes

### `subjects/<subject_id>/report/`
- Role: `source`
- Contains: subject reports
- Keep: yes

---

## 4.2 Batch source layer / batch 源结果层

### `batch/merged/`
- Role: `source`
- Produced by: `summarize_bandpower_outputs`
- Typical contents:
  - `all_subjects_scene_level.csv`
  - `all_subjects_pairs_check.csv`
  - `per_subject_recovery_metrics.csv`
- Keep: mandatory

### `batch/qc/`
- Role: `source`
- Produced by: `qc_filter_and_report`
- Typical contents:
  - `qc_exclusion_subjects.csv`
  - `qc_report.json`
  - QC-filtered outputs
- Keep: mandatory

### `batch/analysis/`
- Role: `source`
- Produced by: task-level analysis writers
- Contains:
  - `task*/raw/...`
  - `task*/qc/...`
- Includes original task PNG/tables/reports
- Keep: mandatory

### `batch/paper/`
- Role: `source`
- Produced by: paper-facing plotting/export code
- Contains paper-oriented output sets
- Keep: mandatory

### `batch/audit/`
- Role: `source`
- Produced by: batch audit/report helpers
- Contains method snapshots, path/input/config audit material
- Keep: mandatory

### `batch/reports/`
- Role: mixed (`source` + `navigation-only`)
- Contains:
  - canonical batch report files (`batch_report.md`, `qc_filter_report.md`, etc.)
  - curated navigation notes (`README_CURATED_MAIN.md`, `README_MAIN_ENTRY.md`)
- Keep: yes, but later split mentally into source reports vs navigation notes

---

## 4.3 Curated summary layer / 摘要层

### `batch/descriptive/overall/`
- Role: `curated-summary`
- Produced by: `build_curated_main_outputs`
- Contains selected summary tables/figures/reports based on `AllScene`, `AllPairs`, QC variants
- Keep: yes
- Note: should remain summary-only

### `batch/descriptive/experience/`
- Role: `curated-summary`
- Produced by: `build_curated_main_outputs`
- Contains selected Experience-group summary tables/figures/reports
- Keep: yes
- Note: may currently include copied figure material from older grouped figure roots

### `batch/inferential/overall/`
- Role: `curated-summary`
- Produced by: `build_curated_main_outputs`
- Contains overall inferential summary tables/figures/reports
- Keep: yes

### `batch/inferential/experience/`
- Role: mixed (`curated-summary` + partial `mirror-copy`)
- Produced by: `build_curated_main_outputs`
- Contains:
  - summary-level inferential outputs
  - **mirrored selected task-based experience outputs** via `mirror_task_branch(...)`
- Keep: yes for now
- De-duplication note: this is a high-priority review area because it mixes summary outputs and mirrored task outputs

---

## 4.4 Unified mirror layer / 统一镜像层

### `batch/tables/`
- Role: `mirror-copy`
- Produced by: `build_batch_unified_surface`
- Built by copying from:
  - `merged/`
  - `qc/`
  - `descriptive/.../tables`
  - `inferential/.../tables`
- Keep content for now: yes
- Future target: convert to table index layer instead of full duplicate copies

### `batch/figures/`
- Role: `mirror-copy`
- Produced by: `build_batch_unified_surface`
- Built by copying from:
  - `descriptive/.../figures`
  - `inferential/.../figures`
  - `figures/topo_scene`
  - `figures/group`
  - `figures/recovery`
  - `figures/scene_sequence`
- Keep content for now: yes
- Future target: convert to figure index layer instead of full duplicate PNG copies

### `batch/analysis/task_outputs/`
- Role: `mirror-copy`
- Produced by: `build_batch_unified_surface`
- Built by copying from `batch/analysis/`
- Duplication note: this is a direct mirror of canonical task results
- Future target: likely removable once a stable navigation index exists

---

## 4.5 Navigation layer / 导航层

### `batch/README.md`
- Role: `navigation-only`
- Produced by: `write_batch_navigation_readme`
- Keep: yes

### `batch/reports/README_CURATED_MAIN.md`
- Role: `navigation-only`
- Keep: yes

### `batch/reports/README_MAIN_ENTRY.md`
- Role: `navigation-only`
- Keep: yes

### `descriptive/.../report/README.md`
### `inferential/.../report/README*.md`
- Role: `navigation-only` or lightweight summary explanation
- Keep: yes, but may later be reduced if redundant

---

## 5. Duplication Risk Map / 重复风险图谱

## High risk / 高重复风险
1. `batch/tables/`
2. `batch/figures/`
3. `batch/analysis/task_outputs/`
4. `batch/inferential/experience/` mirrored task content

## Medium risk / 中等重复风险
1. multiple README/report wrappers in `reports/`, `descriptive/*/report/`, `inferential/*/report/`
2. grouped figures copied from older figure roots into curated directories

## Low risk / 低重复风险
1. `merged/`
2. `qc/`
3. `analysis/`
4. `paper/`
5. subject-level outputs

These are source-bearing and should be treated as canonical.

---

## 6. Safe Preservation Rules / 安全保留规则

### Never delete during early de-dup / 初期去重绝不删除
- any PNG under `analysis/`, `paper/`, or subject-level source dirs
- any task table/stat file under `analysis/`
- any merged/QC source tables
- any source report that contains unique interpretation

### Eligible for future replacement / 后续可替换为索引
- repeated copies under `batch/tables/`
- repeated copies under `batch/figures/`
- repeated `analysis/task_outputs/` mirrors

### Replacement model / 替代模式
Instead of copying whole files, future unified layers should generate:
- `tables_index.csv`
- `figures_index.csv`
- `task_output_index.csv`
- lightweight markdown navigation maps

---

## 7. Next Cleanup Targets / 下一步清理目标

### Target A
Document exact file-generation paths for `batch/tables/` and `batch/figures/` mirror copies.

### Target B
Mark `inferential/experience/` subtrees that are true summaries vs mirrored task payloads.

### Target C
Propose code changes so future runs:
- still generate source outputs
- still generate summaries
- stop generating full duplicate mirror layers by default

### Target D
Only after review, decide whether old mirror-copy directories remain as optional export surfaces.

---

## Practical summary / 实际总结

Current `main` output can be understood as:

```text
source layer:
  merged/ qc/ analysis/ paper/ audit/ subjects/

summary layer:
  descriptive/ inferential/

mirror layer:
  tables/ figures/ analysis/task_outputs/

navigation layer:
  batch README + curated readmes + report entry notes
```

The de-duplication direction should therefore be:

> Preserve source + preserve summaries + shrink mirror copies into index layers.

中文对应：

> 保留源层，保留摘要层，把镜像层逐步收缩成索引层。
