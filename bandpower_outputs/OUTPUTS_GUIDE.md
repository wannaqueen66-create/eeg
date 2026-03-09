# OUTPUTS GUIDE

> This file explains the intended output directory structure for the current `main` branch.
>
> 当前 `main` 分支的输出目录正在从旧结构逐步迁移到新的 staged layout。
> 因此你在实际运行结果里，可能会同时看到：
> - 新结构（推荐优先看）
> - 旧结构兼容路径（过渡期保留）

---

# 1. Recommended reading order / 推荐查看顺序

If you only want the most useful places to inspect results, use this order:

如果你只想快速找到最重要的结果，建议按这个顺序看：

1. `runs/current/subjects/`
   - single-subject outputs / 单被试结果
2. `runs/current/batch/merged/`
   - merged cross-subject tables / 跨被试合并表
3. `runs/current/batch/qc/`
   - QC-filtered tables / QC 过滤后结果
4. `runs/current/batch/reports/`
   - summary + methods + QC reports / 汇总与方法说明
5. `runs/current/batch/analysis/`
   - task-oriented inferential outputs / 按 task 分类的正式分析结果
6. `runs/current/batch/paper/`
   - paper-ready figures/tables / 面向论文与汇报的整理输出

---

# 2. High-level structure / 总体结构

Recommended staged layout:

```text
bandpower_outputs/
  runs/
    current/
      subjects/
      batch/
        descriptive/
        inferential/
```

Meaning:
- `subjects/` = per-subject outputs / 单被试输出
- `batch/` = cross-subject outputs / 跨被试输出
- `descriptive/` = cleaner descriptive result surface / 更简洁的描述性结果层
- `inferential/` = cleaner inferential result surface / 更简洁的显著性分析层

---

# 3. Subject-level outputs / 单被试输出

Typical structure:

```text
runs/current/subjects/<subject_id>/
  tables/
  figures/
  qc/
  report/
```

## 3.1 `tables/`
Main single-subject tables.
主要单被试表格。

Common files:
- `<base>_bandpower_roi.csv`
- `<base>_bandpower_summary.csv`
- `<base>_bandpower_tests.csv`
- `<base>_scene_level.csv`

What they mean:
- `bandpower_roi.csv`
  - segment-level main feature table
  - 分段级主结果表
- `bandpower_summary.csv`
  - condition summary statistics
  - 条件汇总统计
- `bandpower_tests.csv`
  - compact key statistical outputs
  - 关键统计摘要
- `scene_level.csv`
  - view-only scene-level table used downstream
  - 仅保留 view 段的场景级表，是后续 batch 分析主输入

## 3.2 `figures/`
Main single-subject figures.
主要单被试图。

Common figure families:
- ROI condition bars
- paired alpha recovery
- occ PSD
- paired scatter
- block comparison
- QC distributions
- beta split
- three-stage chain
- scene sequence
- topoplots

## 3.3 `qc/`
Single-subject QC and audit outputs.
单被试 QC 与审计输出。

Common files:
- `<base>_qc.csv`
- `<base>_pairs_check.csv`
- `<base>_marker_report.csv`

What they mean:
- `qc.csv`
  - per-segment QC metrics
  - 每段 QC 指标
- `pairs_check.csv`
  - view-gray pairing audit table
  - view-gray 配对审计表
- `marker_report.csv`
  - invalid marker transitions (if any)
  - 非法 marker 转移记录（若存在）

## 3.4 `report/`
Human-readable single-subject report and traceability files.
人类可读的单被试摘要与追溯文件。

Common files:
- `<base>_report.md`
- `input_set_path.txt`
- `config_used.json`

---

# 4. Batch-level outputs / 批量输出

Typical structure:

```text
runs/current/batch/
  merged/
  qc/
  reports/
  audit/
  analysis/
  paper/
  descriptive/
  inferential/
```

## 4.1 `merged/`
Cross-subject merged tables.
跨被试合并表。

Common files:
- `all_subjects_scene_level.csv`
- `all_subjects_pairs_check.csv`
- `per_subject_recovery_metrics.csv`
- `global_bandpower_summary.csv` (optional)

Recommended interpretation:
- `all_subjects_scene_level.csv`
  - main modeling table for scene-level EEG analysis
  - 主建模表（scene-level）
- `all_subjects_pairs_check.csv`
  - main recovery / paired analysis table
  - 恢复/配对分析主表
- `per_subject_recovery_metrics.csv`
  - subject-level recovery summary
  - 被试级恢复摘要

## 4.2 `qc/`
QC-filtered batch outputs.
QC 过滤后的批量结果。

Common files:
- `all_subjects_scene_level_qc.csv`
- `all_subjects_pairs_check_qc.csv`
- `qc_exclusion_subjects.csv`
- `qc_scene_valid_counts.csv`

Recommended interpretation:
- `_qc.csv`
  - these are the filtered versions to prioritize for stricter reporting
  - 这些是 QC 后的严格版结果，正式汇报时建议优先看
- `qc_exclusion_subjects.csv`
  - which subjects were excluded and why
  - 哪些被试被排除、原因是什么
- `qc_scene_valid_counts.csv`
  - valid N per scene after QC
  - 每个场景在 QC 后还剩多少有效样本

## 4.3 `reports/`
Batch-level markdown reports.
批量级 markdown 报告。

Common files:
- `summary_report.md`
- `qc_filter_report.md`
- `methods_snapshot.md`

Recommended interpretation:
- `summary_report.md`
  - entry-point index for the whole batch
  - 整个批量输出的索引入口
- `qc_filter_report.md`
  - explains QC rules and filtering results
  - 解释 QC 规则和过滤结果
- `methods_snapshot.md`
  - concise methods/config snapshot for writing and reproducibility
  - 方法与参数快照，适合写作与复现

## 4.4 `audit/`
Machine-readable audit files.
机器可读审计文件。

Common files:
- `qc_report.json`

Recommended interpretation:
- mostly for reproducibility and programmatic inspection
- 主要用于复现与程序化审计，不是优先给人工阅读的入口

---

# 5. Curated redesigned-main outputs / 面向新 main 的简化输出

Typical structure:

```text
runs/current/batch/descriptive/
  overall/
  experience/

runs/current/batch/inferential/
  overall/
  experience/
```

## 5.1 `descriptive/overall/`
Main full-sample descriptive tables and reports.
全样本描述性结果入口。

Expected content:
- scene-level descriptive tables
- pair/recovery descriptive tables
- basic report files

## 5.2 `descriptive/experience/`
Experience-group descriptive outputs.
Experience 分组描述性结果入口。

Expected content:
- experience-group scene tables
- experience-group recovery tables
- simplified group summary tables
- selected experience-oriented figures (when available)

## 5.3 `inferential/overall/`
Curated full-sample inferential entry point.
全样本显著性分析的整理入口。

## 5.4 `inferential/experience/`
Curated experience-group inferential entry point.
Experience 分组显著性分析的整理入口。

This layer mirrors selected outputs from the detailed task-based analysis tree,
so users can read the redesigned main outputs without navigating every task folder first.

这一层会镜像 task 分析中的重点输出，
让使用者不必先进入每个 task 目录也能找到新的 main 主输出。

# 6. Analysis outputs / 正式分析输出（详细 task 层）

Typical structure:

```text
runs/current/batch/analysis/
  task1_block2_restart/
  task2_C1W45_block_diff/
  task3_trialindex_lmm/
  task4_core_lmm_suite/
  task5_peakindex_invertedu/
  task6_coremetric_special/
  task7_individual_checks/
```

Each task is expected to contain:

```text
taskX_name/
  raw/
    tables/
    figures/
    reports/
  qc/
    tables/
    figures/
    reports/
```

Meaning:
- `raw/` = pre-QC-filter or less strict data path
- `qc/` = QC-filtered path

## 5.1 Task overview / 各 task 含义

### `task1_block2_restart/`
- Compare Block2 first scene vs mean of Block2 scenes 2–6
- 比较 Block2 第一个场景与 Block2 后 5 个场景均值

### `task2_C1W45_block_diff/`
- Compare one target scene between Block1 and Block2
- 比较某个目标场景在 Block1 和 Block2 的差异

### `task3_trialindex_lmm/`
- TrialIndex / adaptation control model
- 试次顺序 / 适应控制模型

### `task4_core_lmm_suite/`
- Main effect / interaction / trend analysis core
- 主效应 / 交互 / 趋势分析核心层

### `task5_peakindex_invertedu/`
- PeakIndex / inverted-U test
- PeakIndex / 倒 U 检验

### `task6_coremetric_special/`
- Core-metric robustness-focused special analysis
- 核心指标专项稳健性分析（默认覆盖 O_alpha / O_theta / O_beta / F_theta）

### `task7_individual_checks/`
- Outlier / influence / robustness audit
- 个体异常值 / 影响度 / 稳健性审计

---

# 6. Paper outputs / 论文与汇报输出

Typical structure:

```text
runs/current/batch/paper/
  raw/
    tables/
    figures/
  qc/
    tables/
    figures/
  report/
```

Meaning:
- `paper/raw/` = paper-facing outputs from raw path
- `paper/qc/` = paper-facing outputs from QC-filtered path
- `paper/report/` = curated summary docs for writing/presentation

Recommended use:
- use this layer when preparing slides, manuscript figures, and final reporting
- 做论文图、汇报图、最终交付时优先看这一层

---

# 7. Group classification results / 不同人群分类结果看哪里

Current group classifications mainly include:
- `Experience` / `ExperienceGroup`
- `SportFreq` / `SportFreqGroup`

Recommended places to inspect group results:

## 7.1 Group-level figures
- `batch/analysis/.../figures/`
- `batch/figures/` (legacy-compatible may still appear)

## 7.2 Group-level model tables
- `task3_trialindex_lmm/...`
- `task4_core_lmm_suite/...`
- `task5_peakindex_invertedu/...`
- `task6_coremetric_special/...`（兼容期内也可能看到旧路径 `task6_obeta_special/...`）
- `task7_individual_checks/...`

## 7.3 Group labels in merged data
- `batch/merged/all_subjects_scene_level.csv`
- `batch/merged/all_subjects_pairs_check.csv`
- `batch/merged/per_subject_recovery_metrics.csv`

---

# 8. Legacy-compatible paths / 旧结构兼容路径

During migration, you may still encounter older paths such as:

```text
summary/
summary/fig/
summary/paper_fig/
summary/paper_fig_qc/
analysis-2/
```

Interpretation:
- these are compatibility paths kept during staged migration
- they are not the preferred long-term target structure

在迁移期间保留这些旧路径是为了兼容旧代码，
但它们不是长期推荐的主结构。

---

# 9. Practical advice / 实用建议

## If you only want the final answer quickly
### 如果你只想快速找最终结果
Look here first:
- `runs/current/batch/reports/`
- `runs/current/batch/qc/`
- `runs/current/batch/paper/`

## If you want single-subject troubleshooting
### 如果你要排查单个被试
Look here first:
- `runs/current/subjects/<subject_id>/qc/`
- `runs/current/subjects/<subject_id>/report/`

## If you want statistical model outputs
### 如果你要看正式统计结果
Look here first:
- `runs/current/batch/analysis/task*/raw/`
- `runs/current/batch/analysis/task*/qc/`

---

# 10. Short summary / 一句话总结

The current `main` branch is moving toward a staged output design:
- `subjects/` for per-subject outputs
- `batch/` for merged/QC/report/audit outputs
- `analysis/` for task-based inferential outputs
- `paper/` for manuscript-facing deliverables

当前 `main` 分支正在迁移到如下分阶段输出结构：
- `subjects/`：单被试输出
- `batch/`：合并/QC/报告/审计输出
- `analysis/`：按 task 组织的正式统计输出
- `paper/`：面向论文与汇报的整理输出
