# EEG Bandpower Pipeline / EEG 频段功率分析流水线

MATLAB/EEGLAB pipeline for **VR scene-viewing EEG bandpower analysis**.  
用于 **VR 场景观看实验 EEG 频段功率分析** 的 MATLAB/EEGLAB 流水线。

---

## Table of Contents / 目录
- [1. Overview / 功能概览](#1-overview--功能概览)
- [2. Marker Protocol / Marker 实验流程](#2-marker-protocol--marker-实验流程)
- [3. Architecture / 架构概览](#3-architecture--架构概览)
- [4. Requirements / 运行环境](#4-requirements--运行环境)
- [5. Input Data Assumptions / 输入数据约定](#5-input-data-assumptions--输入数据约定)
- [6. Quick Start / 快速开始](#6-quick-start--快速开始)
- [7. Full Tutorial / 完整使用教程](#7-full-tutorial--完整使用教程)
- [8. Configuration / 配置文件](#8-configuration--配置文件)
- [9. Outputs / 输出结果](#9-outputs--输出结果)
- [10. Tests / 测试](#10-tests--测试)
- [11. Statistical Notes / 统计解释建议](#11-statistical-notes--统计解释建议)
- [12. FAQ / 常见问题](#12-faq--常见问题)

---

## 1. Overview / 功能概览

**English**
- Marker-state-machine segmentation (`adapt / intro / eyes_closed / eyes_open / view / questionnaire / gray / rest`)
- ROI bandpower extraction (`theta / alpha / beta / low_beta / high_beta / low_gamma`)
- Relative-power note: `low_gamma` is exported from the high-band denominator path (default high-band adjusted to include gamma upper bound when needed)
- View-gray pairing analysis + quality control (QC)
- Automatic CSV/table/figure export

**中文**
- 基于 Marker 状态机分段（`adapt / intro / eyes_closed / eyes_open / view / questionnaire / gray / rest`）
- ROI 频段功率计算（`theta / alpha / beta / low_beta / high_beta / low_gamma`）
- view-gray 配对分析与 QC 质量检查
- 自动导出 CSV、汇总表和图像

---

## 2. Marker Protocol / Marker 实验流程

```text
1→2  VR adaptation / VR适应
2→3  Intro / 实验介绍
3→4  Eyes closed baseline / 闭眼静息
4→9  Eyes open baseline / 睁眼基线
(7→8→9) × 6  Scene viewing → small questionnaire → gray
8→5  Block-end questionnaire → rest
(7→8→9) × 6  Second block
8→6  Final questionnaire (end)
```

---

## 3. Architecture / 架构概览

- Entry orchestration / 入口编排：`run_eeg_bandpower_pipeline.m`
- Input parsing / 输入解析：`+pipeline/parse_input.m`
- Config loading / 配置加载：`+pipeline/load_config.m`
- Output preparation / 输出准备：`+pipeline/prepare_output.m`
- Output path resolver / 路径解析：`+pipeline/resolve_output_dir.m`
- Config snapshot / 配置快照：`+pipeline/write_config_snapshot.m`

---

## 4. Requirements / 运行环境

- MATLAB (recommended R2018a+)
- EEGLAB installed and added to MATLAB path

---


## 5. Input Data Assumptions / 输入数据约定

**English**
- The pipeline expects **already-preprocessed** EEGLAB `.set` files.
- Please provide datasets with filtering / re-reference / ICA (if used) already completed.

**中文**
- 本流水线默认输入为 **已完成预处理** 的 EEGLAB `.set` 文件。
- 请直接提供已完成滤波 / 重参考 / ICA（如使用）等步骤的数据。

---

## 6. Quick Start / 快速开始

```matlab
run_eeg_bandpower_pipeline('path/to/data.set');
```

### Marker interval statistics / Marker 间隔统计（用于人工核查）

```matlab
% single file / 单文件
marker_interval_stats('path/to/data.set');

% folder batch / 文件夹批量
marker_interval_stats('path/to/folder');
```

Outputs:
- `<base>_marker_intervals.csv`
- `<base>_marker_transition_summary.csv`
- `ALL_marker_transition_summary.csv` (folder mode)

---

## 7. Full Tutorial / 完整使用教程

1) Prepare `.set` data with valid markers 1~9  
2) Add EEGLAB path
```matlab
addpath('/path/to/eeglab');
```
3) Single-file run
```matlab
run_eeg_bandpower_pipeline('path/to/data.set');
```
4) GUI mode
```matlab
run_eeg_bandpower_pipeline();
```
5) Folder batch mode
```matlab
run_eeg_bandpower_pipeline('path/to/folder');
```
6) Use custom config
```matlab
run_eeg_bandpower_pipeline('path/to/data.set','config.json');
```

---

## 8. Configuration / 配置文件

Use `config.json` to override defaults.

Additional (new) option:
- `strict_structure` (default `true`): if enabled, the pipeline will **error & skip** any dataset whose global segment counts do not match the expected structure (view=12, q_small=12, q_big=2, gray=12, rest=1).

```json
{
  "gray_dur_min": 3,
  "gray_dur_max": 15,
  "quest_dur_min": 5,
  "quest_dur_max": 120,
  "pairing_mode": "strict",
  "verbose": true,
  "log_file": "",
  "output_dir": "",
  "zip_output": false,
  "global_summary": false,
  "global_summary_path": "",
  "roi": {
    "front": ["F3","F4"],
    "par": ["P3","PZ","P4"],
    "occ": ["O1","OZ","O2"]
  },
  "bands": {
    "theta": [4,7],
    "alpha": [8,12],
    "beta": [13,30],
    "low_beta": [13,20],
    "high_beta": [20,30],
    "low_gamma": [30,45],
    "totalBand40": [1,40],
    "totalBand30": [1,30]
  }
}
```

---

## 9. Outputs / 输出结果

Typical structure (default):

```text
<input_set_folder>/bandpower_outputs/
  ├─ <subject_id>/
  │   ├─ csv/
  │   ├─ fig/
  │   └─ qc/
  └─ summary/
```

Notes:
- Subject folder name equals the input `.set` **base filename** (without extension) for easy verification.
- If `output_dir` is set in `config.json`, outputs go to `<output_dir>/...` instead.
- Each subject output folder contains:
  - `input_set_path.txt` (exact input file path)
  - `<base>_report.md` (human-readable summary for sharing)
- The batch-level markdown index is always written to:
  - `summary/summary_report.md`
- `summary/global_bandpower_summary.csv` is created in folder-batch mode when `global_summary=true`.

CSV includes:
- `*_bandpower_roi.csv`
- `*_bandpower_summary.csv`
- `*_bandpower_tests.csv`
- `*_scene_level.csv`
- `*_pairs_check.csv`
- `*_qc.csv`
- `*_marker_report.csv`

Also:
- `config_used.json`
- `summary/global_bandpower_summary.csv` (if `global_summary=true`)
- `*_outputs.zip` (if `zip_output=true`)

---

## 10. Tests / 测试

```matlab
addpath(genpath(pwd));
run('tests/run_smoke_tests.m');
```

---


## 11. Statistical Notes / 统计解释建议

**English**
- `low_gamma (30–45 Hz)` is informative for immersion/arousal hypotheses, but also sensitive to muscle artifacts (jaw/forehead/neck EMG).
- Prefer interpreting gamma together with QC indicators (especially high-frequency ratio and segment quality checks).
- For confirmatory conclusions, prioritize robust paired effects across repeated scenes and subjects rather than single-segment peaks.

**中文**
- `low_gamma (30–45 Hz)` 可用于沉浸/唤醒相关探索，但对肌电伪迹（咬肌、额肌、颈部）更敏感。
- 建议结合 QC 指标共同解读（尤其是高频占比与分段质量检查）。
- 在正式结论中，优先看跨场景、跨被试的稳定配对效应，不建议仅依据单段峰值判断。

---

## 12. FAQ / 常见问题

**Q1: Missing O1/OZ/O2 channels?**  
Check channel names or edit ROI labels in config.

**Q2: No CSV/figures generated?**  
Check marker quality and whether segments are long enough.

**Q3: How to avoid GUI popups?**  
Pass file/folder path directly as the first argument.

