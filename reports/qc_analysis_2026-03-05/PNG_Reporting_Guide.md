# PNG 输出如何汇报（QC）

> 放在 PPT/论文的目标：让图不仅“展示”，而是**承载一个统计结论**。

---

## 1. 通用规则（所有 PNG 通用）

### 1.1 推荐口播结构（20 秒讲完一张图）
1) **图回答的问题**：自变量（WWR/Complexity/Group/TrialIndex/BlockDiff）是什么？因变量（O_theta 等）是什么？
2) **结论一句话**：显著/趋势 + 方向。
3) **证据一句话**：F/df/p 或 t/p（必要时补 estimate/CI）。
4) **稳健性一句话**：是否受离群影响（Task7）。

### 1.2 图注模板（复制即用）
- 图 X. QC 数据下 <指标> 在 <条件/分组> 的变化。统计使用 <LMM ANOVA/配对 t 检验>；报告 F(df1,df2)=…, p=…（Holm 校正 p=…）。误差线含义以图中定义为准。

### 1.3 PPT 排版建议
- “主结果页”每页 1 张图 + 3 行 bullet：
  - 结论（1 行）
  - 证据（1 行：F/p 或 t/p）
  - 解释/备注（1 行：方向、稳健性）

---

## 2. 按 Task 的“该放哪张图 + 怎么讲”

### Task2（scene block diff）
**建议图**：`scene_blockdiff_overview_qc_*.png`
- 讲法：
  - 先说定义：diff=Block2−Block1（目标 scene WWR45_C1）
  - 结论：SportFreq Low 的 O_theta 有下降趋势（未达 0.05）
  - 证据：t 检验 p=0.089（signrank p=0.054）

### Task3（trialindex LMM）
**建议图**：`trialindex_lmm/figures/qc/<grouping>/` 下的 overview
- 讲法：
  - 结论 1：TrialIndex 主效应稳定（O_theta、F_theta 负向；O_alpha 正向）
  - 结论 2：SportFreq 的 O_beta 出现 Group×TrialIndex（两组斜率不同）

### Task4（core LMM factor WWR）
**建议图**：
- overview：`task4_factor_overview_qc_<analysis>.png`
- 关键效应图：
  - O_theta（Complexity 与 WWR 主效应）
  - SportFreq-O_alpha（Complexity×Group 交互）

**讲法（建议主结果从 Task4 开始）**
- O_theta：
  - Complexity 主效应显著（Holm 后仍显著），HighComplexity 更高
  - WWR 主效应显著（Holm 后仍显著），WWR45 相对更高
- SportFreq-O_alpha：
  - Complexity×Group 交互显著且 Holm 后仍显著（p_holm=0.0247）
  - 解释：复杂度效应在两组方向相反（交互“交叉”形态）

### Task5（PeakIndex 倒 U）
**建议图**：`peakindex_overview_qc_*.png` + 对应指标图
- Experience：强调 Complexity×Group（O_theta、F_theta）
- SportFreq：强调 Complexity 主效应为主、交互不显著

### Task6（O_beta special）
**建议图**：`obeta_group_model_compare_qc_*.png`
- 讲法：
  - 先说均值差异（High>Low）
  - 再说结论：LMM 中 Group 主效应不显著，控制变量前后一致

### Task7（individual checks）
**建议图**：`task7_audit_scorecard.png` + 对应指标的 individual_check 图
- 讲法：
  - 说明做了 MAD>3 与 Top-N influence
  - 对 robustness_score=C 的结论，在主结论中使用“趋势/提示”，或放附录

---

## 3. 最常见错误（避免）
- 只念 p 值不讲方向；或只讲方向不报统计证据。
- 把 raw 图当结论来源（本项目要求只看 qc）。
- 交互效应不画/不讲 simple effects（至少要说明“方向相反/差异来自哪个条件”）。
- 不做稳健性提示：当 Task7 给出 C 时，不能写成“确定性结论”。

