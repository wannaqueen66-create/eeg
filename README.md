# EEG Bandpower Pipeline / EEG 频段功率分析流水线

MATLAB/EEGLAB pipeline for **VR scene-viewing EEG bandpower analysis**.  
用于 **VR 场景观看实验 EEG 频段功率分析** 的 MATLAB/EEGLAB 流水线。

---

## Table of Contents / 目录
- [1. Overview / 功能概览](#1-overview--功能概览)
- [2. Marker Protocol / Marker 实验流程](#2-marker-protocol--marker-实验流程)
- [3. Architecture / 架构概览](#3-architecture--架构概览)
- [4. Requirements / 运行环境](#4-requirements--运行环境)
- [5. Quick Start / 快速开始](#5-quick-start--快速开始)
- [6. Full Tutorial / 完整使用教程](#6-full-tutorial--完整使用教程)
- [7. Configuration / 配置文件](#7-configuration--配置文件)
- [8. Outputs / 输出结果](#8-outputs--输出结果)
- [9. Tests / 测试](#9-tests--测试)
- [10. FAQ / 常见问题](#10-faq--常见问题)

---

## 1. Overview / 功能概览

**English**
- Marker-state-machine segmentation (`adapt / intro / eyes_closed / eyes_open / view / questionnaire / gray / rest`)
- ROI bandpower extraction (`theta / alpha / beta / low_beta / high_beta / low_gamma`)
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

## 5. Quick Start / 快速开始

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

## 6. Full Tutorial / 完整使用教程

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

## 7. Configuration / 配置文件

Use `config.json` to override defaults.

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

## 8. Outputs / 输出结果

Typical structure:

```text
output_dir/subject_id/
  ├─ csv/
  ├─ fig/
  └─ qc/
```

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
- `global_bandpower_summary.csv` (if `global_summary=true`)
- `*_outputs.zip` (if `zip_output=true`)

---

## 9. Tests / 测试

```matlab
addpath(genpath(pwd));
run('tests/run_smoke_tests.m');
```

---

## 10. FAQ / 常见问题

**Q1: Missing O1/OZ/O2 channels?**  
Check channel names or edit ROI labels in config.

**Q2: No CSV/figures generated?**  
Check marker quality and whether segments are long enough.

**Q3: How to avoid GUI popups?**  
Pass file/folder path directly as the first argument.

