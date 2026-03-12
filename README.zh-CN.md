# EEG 频段功率分析流水线（中文说明）

> 该文件为纯中文文档；英文+中文双语请见 `README.md`。
>
> 当前 `main` 分支仅保留**纯 EEG 的 MATLAB/EEGLAB 主流程**。
> EEG × 眼动融合相关内容已单独保存在分支 `feature/eeg-eye-integration`。

## 目录
- [项目简介](#项目简介)
- [功能概览](#功能概览)
- [实验 Marker 流程](#实验-marker-流程)
- [运行环境](#运行环境)
- [输入数据约定](#输入数据约定)
- [快速开始](#快速开始)
- [完整教程](#完整教程)
- [配置文件说明](#配置文件说明)
- [输出结果](#输出结果)
- [测试](#测试)
- [统计解释建议](#统计解释建议)
- [常见问题](#常见问题)

## 项目简介
本项目用于 VR 场景观看实验中的 EEG 频段功率分析，基于 MATLAB + EEGLAB。

## 功能概览
- Marker 状态机分段（adapt / intro / eyes_closed / eyes_open / view / questionnaire / gray / rest）
- ROI 频段功率分析（theta / alpha / beta / low_beta / high_beta / low_gamma）
- 相对功率说明：`low_gamma` 使用高频分母路径导出（默认会在需要时自动扩展 high-band 分母上界以覆盖 gamma）
- view-gray 配对分析 + QC 质量控制
- 自动导出 CSV、图表、汇总表

## 实验 Marker 流程
```text
1→2  VR适应
2→3  实验介绍
3→4  闭眼静息
4→9  睁眼基线
(7→8→9) × 6  场景观看→小问卷→灰屏
8→5  组末大问卷→休息
(7→8→9) × 6  第二组
8→6  终止大问卷
```

## 运行环境
- MATLAB（建议 R2018a 及以上）
- EEGLAB 已安装并加入 MATLAB 路径


## 输入数据约定

- 本流水线默认输入为 **已完成预处理** 的 EEGLAB `.set` 文件。
- 请直接提供已完成滤波 / 重参考 / ICA（如使用）等步骤的数据。

## 快速开始
```matlab
run_eeg_bandpower_pipeline('path/to/data.set');
```

### Marker 间隔统计
```matlab
marker_interval_stats('path/to/data.set');          % 单文件
marker_interval_stats('path/to/folder');            % 文件夹批量
```

## 完整教程
1. 准备 `.set` 数据并确认 marker 结构正常（1~9）
2. 添加 EEGLAB 路径：
```matlab
addpath('/path/to/eeglab');
```
3. 运行单文件：
```matlab
run_eeg_bandpower_pipeline('path/to/data.set');
```
4. GUI 模式：
```matlab
run_eeg_bandpower_pipeline();
```
5. 文件夹批量：
```matlab
run_eeg_bandpower_pipeline('path/to/folder');
```
6. 使用配置文件：
```matlab
run_eeg_bandpower_pipeline('path/to/data.set', 'config.json');
```

## 配置文件说明
`config.json` 关键项：
- `strict_structure`（默认 true）：若开启，则只要全局分段计数与理想结构不一致，就会**直接报错并跳过该被试**（批量模式下继续处理其他人）。
- `design_path`（默认空）：可选，传入一个 CSV（或包含 CSV 的文件夹），描述每个被试的 trial 顺序/场景信息/以及每个场景的主观评分。支持“宽表 trial01_*”和“长表（含 SubjectID/Block/Position/SceneID 等）”两种格式。
- `pairing_mode`：`strict` 或 `lenient`
- `zip_output`：是否打包输出 zip
- `global_summary`：是否导出批量汇总
- `roi`：ROI 通道定义
- `bands`：频段定义（已支持 `low_gamma`）

## 输出结果
### 当前分阶段输出布局

仓库目前正在迁移到更清晰的分阶段输出结构。
在这一阶段，新结构与部分旧式兼容路径可能会并存。

典型布局：
```text
<.set所在文件夹>/bandpower_outputs/
  └─ runs/
     └─ current/
        ├─ subjects/
        │  └─ <subject_id>/
        │     ├─ tables/
        │     ├─ figures/
        │     ├─ qc/
        │     └─ report/
        └─ batch/
           ├─ merged/
           ├─ qc/
           ├─ reports/
           ├─ audit/
           ├─ analysis/
           └─ paper/
```

说明：
- 每个被试的结果文件夹名 = 输入 `.set` 的**文件名（去掉扩展名）**，便于对照核对。
- 如果在 `config.json` 里设置了 `output_dir`，输出会写到 `<output_dir>/...`。
- 在迁移期间，一些旧式输出仍可能出现在 `summary/`、`fig/`、`paper_fig/` 等旧路径下。对 `main` 分支来说，**`batch/` 才是规范的批量输出主入口**；`summary/` 现在仅作为向后兼容概念，保留在部分 helper 和旧路径中。
- 单被试输出正在迁移到：
  - `subjects/<subject_id>/tables/`
  - `subjects/<subject_id>/figures/`
  - `subjects/<subject_id>/qc/`
  - `subjects/<subject_id>/report/`
- 批量汇总与 QC 输出正在迁移到：
  - `batch/merged/`
  - `batch/qc/`
  - `batch/reports/`
  - `batch/audit/`
- analysis task 输出正在迁移到：
  - `batch/analysis/task*/{raw|qc}/...`
- 面向论文/汇报的输出正在迁移到：
  - `batch/paper/{raw|qc}/...`

常见表格包括：
- `*_bandpower_roi.csv`
- `*_bandpower_summary.csv`
- `*_bandpower_tests.csv`
- `*_scene_level.csv`
- `*_pairs_check.csv`
- `*_qc.csv`
- `*_marker_report.csv`
- `*_scene_topo_long.csv`（单被试场景级 topo 长表，供组水平聚合）
- `*_scene_topogrid_layout_block*.csv`（2x3 topo 网格布局清单，用于核对 scene / WWR / Complexity 排位）
- `all_subjects_scene_level.csv`
- `all_subjects_pairs_check.csv`
- `per_subject_recovery_metrics.csv`
- `global_bandpower_summary.csv`（当 `global_summary=true`）

常见图像包括：
- `*_scene_topogrid_<band>_block1.png`
- `*_scene_topogrid_<band>_block2.png`
- `group_scene_topogrid_<band>_block1.png`
- `group_scene_topogrid_<band>_block2.png`
- `experience_low_scene_topogrid_<band>_block1.png`
- `experience_high_scene_topogrid_<band>_block1.png`
- `experience_highminuslow_scene_topogrid_<band>_block1.png`
- `overall_trialindex_response_<tag>.png`
- `overall_trialindex_response_by_order_<tag>.png`
- `experience_trialindex_response_<tag>.png`

TrialIndex 说明：
- TrialIndex 表示呈现顺序（1..12）。
- 如果通过 design mapping 附加了 `Order`，task3 的 TrialIndex 模型还会额外检验 `Order` 和 `Order × TrialIndex`，用来区分顺序效应与 counterbalance 顺序效应。

常见报告包括：
- `<base>_report.md`
- `batch_report.md`（规范主文件名）
- `summary_report.md`（旧兼容别名）
- `qc_filter_report.md`
- `methods_snapshot.md`

另外还有：
- `config_used.json`
- `input_set_path.txt`
- `*_outputs.zip`（当 `zip_output=true`）

## 测试
```matlab
addpath(genpath(pwd));
run('tests/run_smoke_tests.m');
```


## 统计解释建议
- `low_gamma (30–45 Hz)` 可用于沉浸感/唤醒相关探索，但对肌电伪迹更敏感（咬肌、额肌、颈部）。
- 建议将 gamma 指标与 QC 一起解释，特别是高频占比、边界段标记、异常配对比例。
- 对论文主结论，优先使用跨被试、跨场景重复稳定的配对效应，不建议单段峰值直接下结论。

## 常见问题
1. **缺少 O1/OZ/O2 通道？**
   - 检查通道命名，或在配置中修改 ROI。
2. **没有导出图表/CSV？**
   - 检查 marker 结构是否正确、片段是否足够长。
3. **如何关闭 GUI？**
   - 直接传入文件或文件夹路径。
