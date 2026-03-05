# Task8 Group Regulation Modes（QC）

目标：新增一个“专门判断不同人群是否存在不同神经调节方式”的任务。

核心判据：
- 速度差异：Task3 的 Group×TrialIndex（斜率差异）
- 条件调节差异：Task4 的 Group 相关交互（Complexity×Group、WWR×Group、三阶交互）
- 曲线形态差异：Task5 的 PeakIndex Group×Complexity（仅 O_theta/F_theta）
- 稳健性校正：Task7 的 robustness_score 与 direction_flip

输出表：
- `tables/qc/group_regulation_evidence_qc.csv`
- `tables/qc/group_regulation_scorecard_qc.csv`

## experience

- **O_theta**：Level=Weak（score=1.00，robustness=C，flip=1）；关键信号：task4_model2:Complexity:Group(p=0.0302); task5:Group×Complexity (PeakIndex)(p=0.0434)
- **F_theta**：Level=Weak（score=1.00，robustness=C，flip=0）；关键信号：task5:Group×Complexity (PeakIndex)(p=0.0444)
- **O_alpha**：Level=Weak（score=0.50，robustness=A，flip=0）
- **O_beta**：Level=Weak（score=-0.50，robustness=C，flip=0）

## sportfreq

- **O_theta**：Level=Weak（score=-0.50，robustness=C，flip=1）
- **F_theta**：Level=Weak（score=-1.00，robustness=C，flip=1）
- **O_alpha**：Level=Moderate（score=2.00，robustness=C，flip=0）；关键信号：task4_model2:Complexity:Group(p=0.00411)
- **O_beta**：Level=Weak（score=0.50，robustness=C，flip=0）；关键信号：task3:Group×TrialIndex(p=0.0397)

## 结论（可直接汇报）

- 最明确的“人群神经调节方式差异”证据：**sportfreq 的 O_alpha**（Task4: Complexity×Group，Holm 后仍显著）。
- “调节速度差异”证据：**sportfreq 的 O_beta**（Task3: Group×TrialIndex 显著）。
- Experience 下 O_theta/F_theta 的组差异更多体现在 PeakIndex 交互（Task5），但需结合稳健性审计谨慎表述。
