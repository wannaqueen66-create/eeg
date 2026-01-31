# EEG Bandpower Pipeline

A MATLAB/EEGLAB pipeline for **VR 场景观看实验 EEG 频段功率分析**。

## 功能概览
- Marker 状态机分段（adapt / transition / eyes_closed / eyes_open / view / questionnaire / gray / rest）
- ROI 频段功率（theta / alpha / beta / low/high beta）
- 视图/灰屏配对分析 + QC 质量检查
- 自动输出 CSV 与图表

## 依赖
- MATLAB
- EEGLAB（已加入路径）

## 使用方法

### 1) 单文件
```matlab
run_eeg_bandpower_pipeline('path/to/data.set');
```

### 2) 批量处理文件夹
```matlab
run_eeg_bandpower_pipeline('path/to/folder');
```

### 3) 使用配置文件
```matlab
run_eeg_bandpower_pipeline('path/to/data.set', 'config.json');
```

> 若不传参数，将弹窗选择 `.set` 文件。

## 配置文件（config.json）
```json
{
  "gray_dur_min": 3,
  "gray_dur_max": 15,
  "quest_dur_min": 5,
  "quest_dur_max": 120,
  "pairing_mode": "strict",
  "verbose": true,
  "log_file": ""
}
```

## 输出
- CSV：bandpower_roi / bandpower_summary / bandpower_tests / scene_level / pairs_check / qc
- 图表：ROI 条件柱状、配对散点、topoplot、PSD、block 对比、QC 分布等

## 提示
- `pairing_mode` 支持 `strict` 或 `lenient`
- `log_file` 可设置日志文件路径，留空则不写日志
