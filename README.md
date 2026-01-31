# EEG Bandpower Pipeline

MATLAB/EEGLAB pipeline for **VR 场景观看实验 EEG 频段功率分析**。

---

## 目录
- [功能概览](#功能概览)
- [运行环境](#运行环境)
- [快速开始](#快速开始)
- [完整使用教程](#完整使用教程)
- [配置文件](#配置文件)
- [输出结果](#输出结果)
- [常见问题](#常见问题)

---

## 功能概览
- Marker 状态机分段（adapt / transition / eyes_closed / eyes_open / view / questionnaire / gray / rest）
- ROI 频段功率（theta / alpha / beta / low/high beta）
- view-gray 配对分析 + QC 质量检查
- 自动输出 CSV 与图表


## Marker 流程（实验）
```
1→2  VR适应
2→3  实验介绍
3→4  闭眼静息
4→9  睁眼基线
(7→8→9) × 6  场景观看→小问卷→灰屏
8→5  大问卷（组末）→休息
(7→8→9) × 6  第二组
8→6  大问卷（终止）
```

## 运行环境
- MATLAB (R2018a 及以上建议)
- EEGLAB 已加入 MATLAB 路径

---

## 快速开始
```matlab
run_eeg_bandpower_pipeline('path/to/data.set');
```

---

## 完整使用教程

### 1) 准备数据
- EEG 数据需为 **EEGLAB .set** 格式
- Marker 需符合实验流程（1~9）

### 2) 加载 EEGLAB
在 MATLAB 中确保 EEGLAB 已加入路径，例如：
```matlab
addpath('/path/to/eeglab');
```

### 3) 单文件分析（无 GUI）
```matlab
run_eeg_bandpower_pipeline('path/to/data.set');
```

### 4) GUI 模式（弹窗选择文件/文件夹）
直接运行：
```matlab
run_eeg_bandpower_pipeline();
```
弹窗可选择 **单个 .set 文件** 或 **包含 .set 的文件夹**。

### 5) 批量分析（文件夹）
```matlab
run_eeg_bandpower_pipeline('path/to/folder');
```
脚本会自动遍历目录下所有 `.set` 文件并依次处理。

### 6) 使用配置文件
```matlab
run_eeg_bandpower_pipeline('path/to/data.set','config.json');
```

---

## 配置文件
配置文件为 `config.json`，用于覆盖默认参数：
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
    "totalBand40": [1,40],
    "totalBand30": [1,30]
  }
}
```

参数说明：
- `gray_dur_min/max`：灰屏时间合理范围（秒）
- `quest_dur_min/max`：问卷时间合理范围（秒）
- `pairing_mode`：`strict` 或 `lenient`
- `verbose`：是否输出详细日志
- `log_file`：日志文件路径（空则不写日志）
- `output_dir`：输出目录（相对路径或绝对路径）
- `zip_output`：是否打包输出为 zip
- `global_summary`：是否生成批量总表
- `global_summary_path`：总表保存路径（可选）
- `roi`：ROI 通道列表
- `bands`：频段配置

---

## 输出结果
脚本会在数据文件同目录输出（或 output_dir 指定目录）：

### CSV
- `*_bandpower_roi.csv`
- `*_bandpower_summary.csv`
- `*_bandpower_tests.csv`
- `*_scene_level.csv`
- `*_pairs_check.csv`
- `*_qc.csv`
- `*_marker_report.csv`

### 其他
- `global_bandpower_summary.csv`（当 global_summary=true）
- `*_outputs.zip`（当 zip_output=true）

### 图表
- ROI 条件柱状图
- view-gray 配对散点图
- PSD 曲线
- Topoplot（theta/alpha/beta）
- Block 对比图
- QC 分布图等

---

## 常见问题

### Q1: 提示找不到通道 O1/OZ/O2
请确认通道命名是否一致，或修改 ROI 定义。

### Q2: 没有生成图表/CSV
检查 `.set` 文件是否包含正确 marker，以及是否有足够长度的片段。

### Q3: 如何关闭 GUI
直接传入文件路径即可。

---

如需更多定制（ROI/频段可配置、总表汇总、自动打包输出），欢迎提需求。
