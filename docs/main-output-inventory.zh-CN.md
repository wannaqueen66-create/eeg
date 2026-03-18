# main 输出盘点（中文）

## 目录
- [1. 这份盘点的作用](#1-这份盘点的作用)
- [2. 规范 batch 根目录](#2-规范-batch-根目录)
- [3. 输出角色分类](#3-输出角色分类)
- [4. 当前目录逐项盘点](#4-当前目录逐项盘点)
- [5. 哪些地方重复风险最高](#5-哪些地方重复风险最高)
- [6. 安全保留规则](#6-安全保留规则)
- [7. 下一步清理方向](#7-下一步清理方向)

## 1. 这份盘点的作用

这份文档的目的，是把当前 `main` 分支输出结构里的每一层说清楚，并标记它属于哪一种：

- `source`：源结果层
- `curated-summary`：摘要层
- `mirror-copy`：镜像副本层
- `navigation-only`：导航层

这样后面做去重时，才不会误删：
- PNG 图
- 相关 CSV / tables / stats
- task 原始结果
- QC / audit 证据

---

## 2. 规范 batch 根目录

当前 `main` 的 batch 规范主目录是：

```text
<input_folder>/bandpower_outputs/runs/current/batch/
```

单被试规范目录是：

```text
<input_folder>/bandpower_outputs/runs/current/subjects/<subject_id>/
```

后续讨论输出结构时，都以这两个根为准。

---

## 3. 输出角色分类

### source（源结果层）
真正承载科研结果的原始输出。必须保留。

### curated-summary（摘要层）
为了让人更容易读而整理出来的摘要结果。可以保留，但不应该变成第二套完整主树。

### mirror-copy（镜像副本层）
把别处已有结果再复制一份，方便导航，但容易膨胀和混版本。

### navigation-only（导航层）
README / 目录索引 / file index / task count 一类，不直接承载科研结果。

---

## 4. 当前目录逐项盘点

## 4.1 单被试目录

### `subjects/<subject_id>/tables/`
- 角色：`source`
- 含义：单被试原始表格输出
- 处理：必须保留

### `subjects/<subject_id>/figures/`
- 角色：`source`
- 含义：单被试原始图像输出
- 处理：必须保留

### `subjects/<subject_id>/qc/`
- 角色：`source`
- 含义：单被试 QC 与辅助文件
- 处理：必须保留

### `subjects/<subject_id>/report/`
- 角色：`source`
- 含义：单被试报告
- 处理：必须保留

---

## 4.2 batch 源结果层

### `batch/merged/`
- 角色：`source`
- 来源：`summarize_bandpower_outputs`
- 典型文件：
  - `all_subjects_scene_level.csv`
  - `all_subjects_pairs_check.csv`
  - `per_subject_recovery_metrics.csv`
- 处理：必须保留

### `batch/qc/`
- 角色：`source`
- 来源：`qc_filter_and_report`
- 典型文件：
  - `qc_exclusion_subjects.csv`
  - `qc_report.json`
  - QC 过滤后的结果
- 处理：必须保留

### `batch/analysis/`
- 角色：`source`
- 含义：task 原始结果树
- 内容：
  - `task*/raw/...`
  - `task*/qc/...`
  - 其中含图、表、报告
- 处理：必须保留

### `batch/paper/`
- 角色：`source`
- 含义：论文/展示导向源结果
- 处理：必须保留

### `batch/audit/`
- 角色：`source`
- 含义：审计、快照、方法记录
- 处理：必须保留

### `batch/reports/`
- 角色：混合层（`source` + `navigation-only`）
- 说明：这里既有真正的 batch 报告，也有导航性 README
- 处理：保留，但后续可进一步区分“结果报告”和“导航说明”

---

## 4.3 摘要层

### `batch/descriptive/overall/`
- 角色：`curated-summary`
- 含义：overall 描述性摘要
- 内容：表、图、说明
- 处理：保留

### `batch/descriptive/experience/`
- 角色：`curated-summary`
- 含义：experience 分组描述性摘要
- 处理：保留
- 注意：这里有一部分内容是从旧 grouped figures 根里 copy 过来的，需要后续核查

### `batch/inferential/overall/`
- 角色：`curated-summary`
- 含义：overall 推断统计摘要
- 处理：保留

### `batch/inferential/experience/`
- 角色：混合层（`curated-summary` + 部分 `mirror-copy`）
- 含义：experience 推断统计摘要
- 注意：这里还会通过 `mirror_task_branch(...)` 把部分 task 输出镜像进来
- 处理：先保留，但这是后续重点去重区域

---

## 4.4 统一镜像层

### `batch/tables/`
- 角色：`mirror-copy`
- 来源：`build_batch_unified_surface`
- 本质：把 `merged/`、`qc/`、`descriptive/`、`inferential/` 的表再拷一份
- 处理：现阶段保留，后续优先改成索引层

### `batch/figures/`
- 角色：`mirror-copy`
- 来源：`build_batch_unified_surface`
- 本质：把 `descriptive/`、`inferential/`、`topo`、`group`、`recovery`、`scene_sequence` 的图再拷一份
- 处理：现阶段保留，后续优先改成索引层

### `batch/analysis/task_outputs/`
- 角色：`mirror-copy`
- 来源：`build_batch_unified_surface`
- 本质：把 `batch/analysis/` 再镜像一份
- 处理：后续高度怀疑属于冗余镜像层

---

## 4.5 导航层

### `batch/README.md`
- 角色：`navigation-only`
- 处理：保留

### `batch/reports/README_CURATED_MAIN.md`
- 角色：`navigation-only`
- 处理：保留

### `batch/reports/README_MAIN_ENTRY.md`
- 角色：`navigation-only`
- 处理：保留

### `descriptive/*/report/README.md`
### `inferential/*/report/README*.md`
- 角色：`navigation-only` 或轻量摘要说明
- 处理：保留，但后续可收缩重复说明

---

## 5. 哪些地方重复风险最高

### 高风险
1. `batch/tables/`
2. `batch/figures/`
3. `batch/analysis/task_outputs/`
4. `batch/inferential/experience/` 中镜像进来的 task 内容

### 中风险
1. 多层 README / report
2. 某些 curated 目录里 copy 过来的 grouped figures

### 低风险
1. `merged/`
2. `qc/`
3. `analysis/`
4. `paper/`
5. 单被试目录

这些低风险目录应该视为规范源结果层。

---

## 6. 安全保留规则

### 初期去重绝不删
- `analysis/` 下的 PNG
- `paper/` 下的 PNG
- subject 级 PNG
- task 表格 / stats 文件
- merged / qc 主表
- 唯一报告文件

### 后续可改索引化
- `batch/tables/`
- `batch/figures/`
- `batch/analysis/task_outputs/`

### 未来更合理的替代物
以后不再整份复制，而是生成：
- `tables_index.csv`
- `figures_index.csv`
- `task_output_index.csv`
- 统一导航 markdown

---

## 7. 下一步清理方向

### 方向 A
继续核查 `batch/tables/` 和 `batch/figures/` 的复制路径，确认它们到底复制了哪些内容。

### 方向 B
把 `inferential/experience/` 中“真正摘要”和“镜像 task 内容”拆开看。

### 方向 C
后续代码改造方向：
- 继续保留源结果
- 继续保留摘要层
- 停止默认生成整份镜像层
- 改为生成索引层

---

## 当前总结

当前 `main` 的输出结构，可以简单理解为：

```text
源结果层（默认主入口）：
  merged / qc / reports / analysis / paper / audit / subjects

摘要层（仅在显式开启 curated 时出现）：
  descriptive / inferential

镜像层（仅在显式开启 unified mirror surface 时出现）：
  tables / figures / analysis/task_outputs

导航层：
  README / curated readme / entry notes
```

所以后续 main 去重的正确方向就是：

> 保留源层，保留可选摘要层，把镜像层默认关闭，并逐步改成索引层。
