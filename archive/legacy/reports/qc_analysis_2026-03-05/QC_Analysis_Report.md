# EEG 数据分析报告（QC 结果｜数据挖掘版）

> 路径：`/root/analysis-2-EEG/`  
> 说明：本报告**仅汇总 QC（qc）结果**；`raw` 目录/图表一律不讨论。  
> 指标：O_theta、F_theta、O_alpha、O_beta（均来自 QC 管线输出）。

---

## 0. 你说的“全面”我怎么理解（本报告的写法原则）

你要的“全面”= **把 QC 数据里“挖掘出来的可解释信息”都写出来**，而不是把所有不显著项逐条罗列。

因此本报告遵循：
- **写**：显著主效应、显著交互、稳定趋势（例如 TrialIndex）、以及“差一点显著但形态清晰、且可在讨论里作为趋势”的点。
- **不写/少写**：明显不显著、且方向/形态也不提供额外信息的结果（避免噪声淹没重点）。
- **每个 Task 都补一段“可能反映的心理/神经过程”**（你说的“结论”：这些数据结果可以对应反映什么）。
- **稳健性**：若 Task7 指示某些结论对少数被试敏感（robustness_score=C / direction_flip=1），会在解释段中降级措辞（“提示/趋势/需谨慎”）。

> 重要提醒：EEG 频段的功能解释属于“可解释框架”，不是直接因果。下面用“可能反映/提示”表述，并优先基于你这套任务因子（WWR、Complexity、Group、TrialIndex/Block）来解释。

---

## 1. 数据与方法概览（QC）

### 1.1 数据范围
- 分组维度：Experience（High/Low）、SportFreq（High/Low）
- 任务因子：WWR（15/45/75）、Complexity（Low/High）、Block（Block1/Block2）、TrialIndex（1–12）

### 1.2 统计框架（来自本次输出）
- Task1/Task2：差值（diff）+ 配对检验（见各自 qc 报告 md）
- Task3：TrialIndex LMM：`EEG ~ Group*WWR*Complexity + TrialIndex + Group:TrialIndex + (1|Subject)`
- Task4：核心 LMM（主效应/二阶/三阶），Satterthwaite ANOVA（输出 CSV 含 pValue/p_holm）
- Task5：PeakIndex（倒 U 指标）LMM：`PeakIndex ~ Complexity*Group + (1|Subject)`
- Task6：O_beta special（Group 效应稳健性）：
  - A：`EEG ~ Group + (1|Subject)`
  - B：`EEG ~ WWR + Complexity + Group + (1|Subject)`
- Task7：个体/离群影响审计（MAD>3，Top-N influence，robustness_score）

### 1.3 多重校正说明
- Task4/Task6/部分表同时提供 `pValue` 与 `p_holm`（Holm 校正）。
- 本报告：
  - **主结论优先使用 Holm 后仍显著的结果**；
  - 若仅原始 p<0.05、但 Holm 后不显著，会标注为“未校正显著/校正后不显著”，可在讨论里当“探索性发现”。

---

## 2. Task1：Block2 首次 vs 其余（restart/适应）

### 2.1 分析定义
- 每名被试、每个指标：`diff = block2_first − block2_rest_mean`
- QC 输出：
  - 报告：`task1_block2_restart/reports/qc/block2_restart_report_qc.md`
  - 表：`task1_block2_restart/tables/qc/block2_restart_subjectlevel_qc.csv`
  - 图：`task1_block2_restart/figures/qc/`

### 2.2 数据挖掘发现（QC，方向与量级）
- Experience 分组（均值 diff）：
  - Low：O_theta +0.001185；F_theta +0.000465；O_alpha −0.009363；O_beta −0.013162
  - High：O_theta −0.002692；F_theta −0.000496；O_alpha −0.011003；O_beta −0.001832
- SportFreq 分组（均值 diff）：
  - Low：O_theta +0.002908；F_theta +0.000477；O_alpha −0.018900；O_beta −0.019700
  - High：O_theta −0.005045；F_theta −0.000728；O_alpha −0.003311；O_beta +0.006315

> 是否显著请以 task1 的 qc 报告 md 里的检验结果为准；此处重点把“首段 vs 后段”的变化方向与量级挖出来。

### 2.3 结论（Task1）
- Task1 在 QC 上呈现出“Block2 内部首段与后段存在差值”的现象，但不同分组/指标方向不完全一致。

### 2.4 可能反映的内容（解释口径）
- 这一 Task 更像是在检测 **restart/重新进入任务后的瞬时状态**（首段）与 **短时稳定状态**（后段）的差异。
- 若 O_alpha diff 为负（首段 < 后段），常见解释是：**首段注意/投入更强（alpha 更抑制），随后 alpha 回升**（更放松/更自动化）。
- 由于组间方向不稳定、且需要结合显著性与 Task7 稳健性，本 Task 的解释建议放在“现象描述/补充”层级，而不要单独作为主结果。

---

## 3. Task2：特定场景（WWR45_C1）Block2−Block1 差值

### 3.1 分析定义
- 目标场景 WWR45_C1；每名被试：`diff = Block2 − Block1`
- QC 输出：
  - 报告（all/experience/sportfreq）：`task2_C1W45_block_diff/reports/qc/**/scene_blockdiff_report_qc.md`
  - 表：`task2_C1W45_block_diff/tables/qc/**/scene_blockdiff_subjectlevel_qc.csv`
  - 图：`task2_C1W45_block_diff/figures/qc/**/`

### 3.2 数据挖掘发现（QC）
（来自 qc 报告 md 的配对检验）
- SportFreq Low：
  - O_theta：mean(diff) = −0.01365，t(10)=−1.883，p=0.0891（signrank p=0.0537）→ **下降趋势**
- 其余指标/分组：未见稳定显著。

### 3.3 结论（Task2）
- 在 WWR45_C1 场景下，SportFreq Low 组的 O_theta 存在从 Block1 到 Block2 的**下降趋势**（探索性）。

### 3.4 可能反映的内容（解释口径）
- Task2 的 block 差值通常对应：
  - **熟悉化/学习（novelty 下降）**
  - **任务策略改变（更自动化）**
  - 或 **疲劳/资源下降**
- O_theta 的下降趋势可以被解释为：在第二个 block 中，该组在同一场景下的 **认知控制/新奇反应需求降低**（更熟练或更省力）。

---

## 4. Task3：TrialIndex 适应趋势（LMM）

### 4.1 模型
`EEG ~ Group*WWR*Complexity + TrialIndex + Group:TrialIndex + (1|Subject)`

### 4.2 数据挖掘发现（QC，显著趋势）

**Experience 分组**
- O_theta：TrialIndex 负向显著（Estimate = −0.001558，p=0.000178）
- F_theta：TrialIndex 负向显著（Estimate = −0.001073，p=0.000265）
- O_alpha：TrialIndex 正向显著（Estimate = +0.003569，p=0.00286）

**SportFreq 分组**
- O_theta：TrialIndex 负向显著（Estimate = −0.001546，p=1.87e−05）
- F_theta：TrialIndex 负向显著（Estimate = −0.001152，p=6.43e−06）
- O_alpha：TrialIndex 正向显著（Estimate = +0.004132，p=4.56e−05）
- O_beta：
  - TrialIndex 正向显著（Estimate = +0.002475，p=0.0239）
  - Group×TrialIndex 负向显著（Estimate = −0.003053，p=0.0397）
    - 含义：Low 组斜率约 +0.002475；High 组斜率约 +0.002475−0.003053≈ −0.00058（两组随试次变化方向/强度不同）

### 4.3 结论（Task3）
- **最稳定的“数据挖掘”发现**：随 TrialIndex 推进，
  - O_theta、F_theta **系统性下降**；
  - O_alpha **系统性上升**。
- SportFreq 的 O_beta 还提示：两组的“随试次变化斜率”不同。

### 4.4 可能反映的内容（解释口径）
- 这一组结果最像经典的 **短时适应/熟悉化/策略自动化**：
  - theta（尤其前额/中线 theta 常被用作控制/努力指标）随 trial 下降 → **控制需求下降/更自动化** 或 **资源下降/疲劳**（需结合行为/主观指标判别）。
  - alpha 上升 → **注意投入降低** 或 **抑制增强/放松增加**；在“熟悉化”解释里通常对应“更省力”。
- SportFreq 的 O_beta 组间斜率差异提示：
  - 两组在训练/运动背景相关的适应轨迹可能不同（例如某组更快稳定、另一组逐渐增强/减弱）。这类结论建议结合 Task7 稳健性与 Task4/5 的交互一起讲，避免过度单点解释。

---

## 5. Task4：核心 LMM（WWR/Complexity/Group）——主结果面板

### 5.1 Experience 分组（QC）

**O_theta（主效应稳健）**
- Complexity：F=21.78，p=4.89e−06（Holm 后仍显著 p_holm=1.47e−05）
- WWR：F=5.65，p=0.00397（Holm 后仍显著 p_holm=0.00795）
- 方向（direction_means）：
  - ComplexityHigh > ComplexityLow（0.1107 vs 0.1030）
  - WWR45 相对更高（0.1110），15/75 较低（0.1041/0.1058）

**F_theta（WWR 主效应；交互为探索性）**
- WWR：F=4.56，p=0.0113（Holm 后仍显著 p_holm=0.0340）
- WWR×Complexity：p=0.00885，但 Holm 后边缘 p_holm=0.0531（探索性）
- 形态：在 WWR45 条件下，高复杂度的 F_theta 更高（0.1080 vs 0.1012）

### 5.2 SportFreq 分组（QC）

**O_theta（主效应稳健）**
- Complexity：p=4.88e−06（p_holm=1.47e−05）
- WWR：p=0.003995（p_holm=0.00799）

**F_theta（WWR 主效应稳健；交互为探索性）**
- WWR：p=0.01138（p_holm=0.0341）
- WWR×Complexity：p=0.00888，但 Holm 后边缘 p_holm=0.0533

**O_alpha（关键交互：稳健）**
- Complexity×Group：F=8.38，p=0.00411，且 Holm 后仍显著 p_holm=0.0247
- 方向（折叠 WWR 的 Complexity×Group 均值）：
  - Group=High：HighComplexity (0.11095) > LowComplexity (0.10404)
  - Group=Low：HighComplexity (0.10055) < LowComplexity (0.12040)
  - 典型“交叉交互”：复杂度效应方向在两组相反。

### 5.3 结论（Task4）
- **O_theta 对 Complexity 与 WWR 的主效应是本次 QC 最稳健的主结果之一**（两种分组下均 Holm 后显著）。
- SportFreq 分组下，O_alpha 出现 **Complexity×Group 的稳健交互**（Holm 后仍显著），是最值得作为“组差异机制”去解释的发现。

### 5.4 可能反映的内容（解释口径）
- **Complexity → O_theta 增强**：常被解释为任务复杂度/负荷提升导致的 **控制/工作记忆/整合需求增加**。
- **WWR → O_theta/F_theta 的非单调形态（45 更高）**：更像“中等水平刺激/节律/负荷下反应最大”的 **倒 U/最优唤醒** 形态（这点可与 Task5 PeakIndex 呼应）。
- **SportFreq 下 O_alpha 的交叉交互**：可以讲成“不同运动背景组在高复杂度时采用了不同的注意/抑制策略”：
  - 一组在高复杂度时 alpha 上升（可能更强的抑制/选择性过滤/更省力策略），
  - 另一组 alpha 下降（可能更强的外显投入/注意增强）。
  具体“上升/下降代表什么”建议与你任务范式（Complexity 的操控含义）一起在讨论里定调。

---

## 6. Task5：PeakIndex（倒 U）与 Complexity×Group

### 6.1 数据挖掘发现（QC）

**Experience 分组（交互显著）**
- O_theta：
  - Complexity：F=8.61，p=0.00725
  - Group×Complexity：F=4.33，p=0.0484（fixed effects 中交互 p=0.0434）
  - 均值（Complexity×Group，折叠 round）：
    - HighComplexity：Low 组 0.01418 > High 组 0.00592
    - LowComplexity：Low 组约 0（−0.00023）≈ High 组 0.00315
- F_theta：
  - Complexity：F=12.09，p=0.00109
  - Group×Complexity：F=4.28，p=0.0439（fixed effects 中交互 p=0.0444）

**SportFreq 分组（以 Complexity 主效应为主）**
- O_theta：Complexity 主效应 p=0.0226；交互不显著
- F_theta：Complexity 主效应 p=0.00292；交互不显著

### 6.2 结论（Task5）
- Experience 分组下，PeakIndex 显示 **Complexity×Group 的交互**（O_theta 与 F_theta），提示“倒 U/最优点”的形态在不同组/复杂度下不同。
- SportFreq 分组下更偏向 **Complexity 主效应**。

### 6.3 可能反映的内容（解释口径）
- PeakIndex 的直观讲法可以是：
  - “在 WWR(15/45/75) 三水平下，EEG 指标是否呈现倒 U（中间更强）以及倒 U 的强弱/位置”。
- Experience 下的交互可解释为：
  - **在高复杂度情境中，不同经验组对 WWR 的‘最优反应区间’不同**（某组在中等/特定 WWR 下更容易出现峰值反应）。
- 这一 Task 最适合用来支撑 Task4 里“WWR45 更高/倒 U 形态”的叙事：Task4 给主效应，Task5 给“曲线形态证据”。

---

## 7. Task6：O_beta special（Group 差异稳健性）

### 7.1 数据挖掘发现（QC）
- Experience：原始均值 High 0.2214 > Low 0.1861，但 Group 不显著：
  - Model A p=0.366；Model B p=0.366
- SportFreq：原始均值 High 0.2329 > Low 0.1823，但 Group 不显著：
  - Model A p=0.158；Model B p=0.158

### 7.2 结论（Task6）
- O_beta 的组间均值差异目前属于**描述性差异**，在 LMM 下不构成可下“显著组差异”结论。

### 7.3 可能反映的内容（解释口径）
- 汇报时可以这么说：
  - “虽然均值上 High>Low，但在控制 WWR/Complexity 后组差异不显著，说明 O_beta 的组差异证据不足/不稳健。”
  - 这类结论非常适合配合 Task7 的稳健性审计一起讲，避免观众抓着均值差异过度解读。

---

## 8. Task7：个体影响与稳健性审计（QC）

### 8.1 数据挖掘发现（QC）
- Experience：MAD>3 离群点 8/96；robustness_score：
  - O_alpha = A（相对稳健）
  - O_theta/F_theta/O_beta = C（需谨慎）
- SportFreq：MAD>3 离群点 6/96；多为 C

### 8.2 关键提示（可写进 discussion/limitations）
- O_theta（Experience、SportFreq）出现离群敏感：主要 MAD>3 被试包括：宋、谢庭勇、廖思源、张景鸿。
- 若你要把“组差异”当作主结论，建议在正文里加一句：
  - “我们进行了 MAD/影响点审计；部分指标（robustness_score=C）在剔除极端个体后方向/幅度可能改变。”

### 8.3 结论（Task7）
- **稳健性信息本身就是数据挖掘的重要发现**：它告诉我们哪些结论可以强说、哪些只能弱说。

---

## 9. PNG 应该怎么汇报（把“图”变成“结论”）

详细见：`/root/analysis-2-EEG/PNG_Reporting_Guide.md`

一句话规则：
- **每张图必须对应一句结论**（主效应/交互/趋势）
- 并给 **一条统计证据**（F/df/p 或 t/p）
- 再加 **一句解释**（可能反映什么 + 是否稳健）

---

## 10. 总结（直接可用的“结果—含义”写法）

1) **适应/熟悉化轨迹（最稳健）**：TrialIndex 上 O_theta、F_theta 系统性下降，O_alpha 系统性上升，提示随试次推进出现显著适应（更自动化/更省力或疲劳下降）。
2) **任务负荷（Complexity）效应（稳健）**：O_theta 在高复杂度显著更高（两种分组下 Holm 后仍显著），支持复杂度操控确实改变了神经负荷相关指标。
3) **WWR 相关的最优区间/非单调形态（证据链：Task4+Task5）**：O_theta/F_theta 对 WWR 的反应在 45 条件更强，并在 PeakIndex 中体现出倒 U/峰值形态差异，提示存在“中等水平 WWR 反应更强”的最优区间。
4) **组差异机制（最值得讲的交互）**：SportFreq 分组下 O_alpha 存在稳健的 Complexity×Group 交叉交互，提示两组在高复杂度下采用不同的注意/抑制策略。
5) **该弱说的地方**：O_beta 的组间差异在 LMM 中不显著；同时 Task7 提示若干“组差异”结论对少数个体敏感，需要谨慎措辞。

---

## 附：关键文件索引（QC）
- 总览：`/root/analysis-2-EEG/MASTER_REPORT.md`
- 本报告（数据挖掘版）：`/root/analysis-2-EEG/QC_Analysis_Report.md`
- PNG 汇报指南：`/root/analysis-2-EEG/PNG_Reporting_Guide.md`

- Task6 扩展（四个频段 Group 效应汇总）：`/root/analysis-2-EEG/task6_allbands_group_effect_qc.md`
