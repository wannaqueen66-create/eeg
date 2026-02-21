# EEG 频段功率分析流水线（中文说明）

> 该文件为纯中文文档；英文+中文双语请见 `README.md`。

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
- `strict_structure`（默认 true）：若开启，则只要全局分段计数与理想结构不一致（view=12, q_small=12, q_big=2, gray=12, rest=1），就会**直接报错并跳过该被试**（批量模式下继续处理其他人）。
- `pairing_mode`：`strict` 或 `lenient`
- `zip_output`：是否打包输出 zip
- `global_summary`：是否导出批量汇总
- `roi`：ROI 通道定义
- `bands`：频段定义（已支持 `low_gamma`）

## 输出结果
目录结构（默认）：
```text
<.set所在文件夹>/bandpower_outputs/
  ├─ <subject_id>/
  │   ├─ csv/
  │   ├─ fig/
  │   └─ qc/
  └─ summary/
```

说明：
- 每个被试的结果文件夹名 = 输入 `.set` 的**文件名（去掉扩展名）**，便于对照核对。
- 如果在 `config.json` 里设置了 `output_dir`，输出会写到 `<output_dir>/...`。
- 每个被试结果文件夹内会写入：
  - `input_set_path.txt`（记录本次处理的输入文件路径，方便追溯）
  - `<base>_report.md`（可直接发给别人看的 Markdown 摘要）
- 汇总文件夹会始终生成：
  - `summary/summary_report.md`（批量运行的索引/摘要）
- 批量分析用合并表（默认开启 `batch_summaries=true`）会写入 `summary/`：
  - `all_subjects_scene_level.csv`（如果提供 `design_path`，会自动附加 scene_name/WWR/Cond/Complexity 等）
  - `all_subjects_pairs_check.csv`（如果提供 `design_path`，会自动附加 scene_name/WWR/Cond/Complexity 等）
  - `per_subject_recovery_metrics.csv`（如果提供 `design_path`，会带 SportFreq/Experience 分组标签）
- 文件夹批量模式下当 `global_summary=true` 时，会额外生成：
  - `summary/global_bandpower_summary.csv`

主要文件：
- `*_bandpower_roi.csv`
- `*_bandpower_summary.csv`
- `*_bandpower_tests.csv`
- `*_scene_level.csv`
- `*_pairs_check.csv`
- `*_qc.csv`
- `*_marker_report.csv`
- `config_used.json`
- `summary/global_bandpower_summary.csv`（当 `global_summary=true`）
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
