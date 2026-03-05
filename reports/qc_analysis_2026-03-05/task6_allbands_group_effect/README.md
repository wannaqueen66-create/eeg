# Task6（四波段扩展）- Figures（QC）

本目录用于满足“Task6 也要有对应 PNG”的要求。

## 为什么这些 PNG 合法/一致？
- 对 4 个波段做“控制 WWR + Complexity 后的 Group 效应”的 LMM，Task4 的 Model1 已经完成：
  - `EEG ~ WWR + Complexity + Group + (1|Subject)`
- 因此 Task4 的 QC 图（按 experience / sportfreq 输出的每波段模型图）可直接作为 Task6（四波段 Group 效应）对应 PNG。

## PNG 文件说明
- `figures/qc/task6_overview_qc_experience.png`
- `figures/qc/task6_overview_qc_sportfreq.png`
  - 来自 Task4 的 overview，用于展示 4 个指标的整体效应概览。

- `figures/qc/experience/task6_group_effect_qc_experience_<metric>.png`
- `figures/qc/sportfreq/task6_group_effect_qc_sportfreq_<metric>.png`
  - 对应 4 个波段（O_theta、F_theta、O_alpha、O_beta）的模型结果图（QC）。

- `figures/qc/task6_obeta_model_compare_qc_<analysis>.png`
  - 保留原 Task6（O_beta special）“ModelA vs ModelB”对比图：它回答的是更具体的问题（仅组效应 vs 控制后组效应）。

## 统计汇总表
- 四波段 Group 效应汇总（含均值与 F/p/p_holm）：
  - `/root/analysis-2-EEG/task6_allbands_group_effect_qc.md`
