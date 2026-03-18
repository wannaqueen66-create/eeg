# Main Output De-duplication Plan / main 输出去重方案

## Table of Contents / 目录
- [1. Goal / 目标](#1-goal--目标)
- [2. Non-Negotiables / 不可触碰原则](#2-non-negotiables--不可触碰原则)
- [3. Current Output Layers / 当前输出层](#3-current-output-layers--当前输出层)
- [4. Classification Rules / 分层规则](#4-classification-rules--分层规则)
- [5. What to Keep / 保留内容](#5-what-to-keep--保留内容)
- [6. What to De-duplicate / 去重对象](#6-what-to-de-duplicate--去重对象)
- [7. Safe Migration Strategy / 安全迁移策略](#7-safe-migration-strategy--安全迁移策略)
- [8. Immediate Next Steps / 下一步执行项](#8-immediate-next-steps--下一步执行项)

## 1. Goal / 目标

**English**

Reduce structural duplication in the `main` branch output tree **without losing any scientific content**, especially:
- PNG figures
n- CSV/tables/stats
- task-level traceability
- QC and audit evidence

The target is to make `batch/` easier to navigate while preserving all original result-bearing files.

**中文**

在 **不丢失任何科研结果内容** 的前提下，收缩 `main` 分支输出结构中的重复层。重点保护：
- PNG 图像
- CSV / tables / stats 数据
- task 级可追溯性
- QC 与 audit 证据

目标不是删结果，而是让 `batch/` 更容易导航、更不容易混版本。

---

## 2. Non-Negotiables / 不可触碰原则

### Must not delete / 绝对不能删
- `batch/analysis/` 下的 task 原始结果（图、表、report）
- `batch/merged/` 主合并表
- `batch/qc/` 主 QC 结果
- `batch/paper/` 面向论文/展示的源图
- 每个 subject 目录下的原始被试结果
- 唯一 report / stats 文件

### Allowed to simplify later / 后续可收缩对象
- 纯镜像复制出来的 `tables/`、`figures/`、`reports/` 统一入口副本
- 与源结果内容重复的导航性 README / index 文件
- curated 层中与 task 层完全重复的整份拷贝

### Rule / 总规则
**Keep source files, reduce duplicate copies.**  
保留源文件，只减少重复副本。

---

## 3. Current Output Layers / 当前输出层

Current `main` branch batch output contains multiple overlapping layers:

```text
batch/
├─ merged/          # source tables
├─ qc/              # source QC results
├─ reports/         # batch-level reports / navigation
├─ audit/           # audit snapshots
├─ analysis/        # source task outputs
├─ paper/           # source paper-facing outputs
├─ descriptive/     # curated descriptive layer
├─ inferential/     # curated inferential layer
├─ tables/          # aggregated mirror layer
├─ figures/         # aggregated mirror layer
└─ reports/         # partly report, partly mirror/navigation
```

### Problem / 问题
The duplication risk is mainly caused by:
1. task outputs existing in `analysis/`
2. summarized views re-expressed in `descriptive/` and `inferential/`
3. additional copies collected again into `batch/tables/` and `batch/figures/`

---

## 4. Classification Rules / 分层规则

### A. Source Layer / 源结果层
Directories that contain original result-bearing files:
- `batch/merged/`
- `batch/qc/`
- `batch/analysis/`
- `batch/paper/`
- `batch/audit/`
- `subjects/<subject_id>/...`

These are canonical and must be preserved.

### B. Curated Summary Layer / 摘要层
Directories intended for easier reading:
- `batch/descriptive/`
- `batch/inferential/`

These may keep summary tables/figures, but should not grow into a second full copy of the entire task tree.

### C. Mirror / Index Layer / 镜像或索引层
Directories whose main value is convenience rather than uniqueness:
- `batch/tables/`
- `batch/figures/`
- some duplicated `reports/README/index` surfaces

These are the primary targets for de-duplication.

---

## 5. What to Keep / 保留内容

### Keep as canonical source / 作为规范源保留
- task PNGs under `batch/analysis/task*/.../figures/`
- task CSV/stat tables under `batch/analysis/task*/.../tables/`
- batch merged master tables under `batch/merged/`
- batch QC tables and JSON under `batch/qc/`
- paper-facing figures under `batch/paper/`
- per-subject source outputs under `subjects/<subject_id>/...`

### Keep as curated summaries / 作为摘要保留
- selected descriptive summary tables/figures
- selected inferential summary tables/figures
- high-level navigation reports

### Not required to keep as duplicated copies / 不要求保留为重复副本
- the same PNG copied again into a separate aggregate figure mirror directory
- the same CSV copied again purely for convenience

---

## 6. What to De-duplicate / 去重对象

### Priority 1 / 优先级最高
Replace repeated copied artifacts in:
- `batch/tables/`
- `batch/figures/`

with index files such as:
- `tables_index.csv`
- `figures_index.csv`
- `result_map.md`

### Priority 2 / 第二优先级
Reduce repeated report wrappers that restate the same information without adding new interpretation.

### Priority 3 / 第三优先级
Constrain curated layers so they contain:
- summary-level tables
- summary-level figures
- summary-level readmes

but not full duplicated task result trees.

---

## 7. Safe Migration Strategy / 安全迁移策略

### Phase 0: No deletion / 第 0 阶段：不删文件
First produce an inventory only.

### Phase 1: Inventory / 第 1 阶段：盘点
Generate a machine-readable inventory that marks each output as one of:
- `source`
- `curated-summary`
- `mirror-copy`
- `navigation-only`

### Phase 2: Canonical mapping / 第 2 阶段：主入口映射
Publish a single map that tells users:
- where to read original source results
- where to read summaries
- which directories are only convenience mirrors

### Phase 3: Stop future duplication / 第 3 阶段：停止继续制造重复
Modify code paths so future runs create:
- source files in canonical directories
- indexes instead of duplicate copies in mirror directories

### Phase 4: Optional cleanup / 第 4 阶段：可选收缩
Only after inventory validation, consider removing mirror copies that are provably redundant.

---

## 8. Immediate Next Steps / 下一步执行项

1. Add this plan document to the repository.
2. Build a concrete inventory of current `batch/` output roles.
3. Mark directories/files as source vs summary vs mirror.
4. Propose code changes that stop copying duplicated PNG/CSV artifacts.
5. Do **not** delete any result-bearing files until inventory review is accepted.

---

## Practical Decision / 当前实际决策

For the next iteration, `main` should follow this rule:

> **Source results stay in place; summary results stay concise; mirror layers should become index layers instead of full duplicate layers.**

中文对应：

> **源结果层原地保留，摘要层保持精简，镜像层逐步改成索引层，而不是整份复制层。**
