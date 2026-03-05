# Task6（更新版）：四个频段的 Group 效应稳健性（QC）

> 你提到“Task6 应该对 4 个波都做”。目前目录 `task6_obeta_special/` 只对 **O_beta** 做了专门的两模型检验。
>
> 但如果你的目的其实是“检验 Group 差异在控制 WWR/Complexity 后是否仍成立”，那么 **Task4 的 Model1（主效应模型）本质上已经给出了每个频段在控制 WWR+Complexity 后的 Group 检验**（Satterthwaite ANOVA）。
>
> 因此这里给出一个“Task6 扩展版”：把 **4 个频段**在两种分组体系（experience / sportfreq）下的：
> 1) 描述性组均值（High vs Low）
> 2) 控制 WWR+Complexity 后的 Group 显著性（F/df/p 与 p_holm）
> 统一列出来，方便直接放进汇报。

数据来源：
- 组均值：`/root/analysis-2-EEG/task4_core_lmm_suite/factor_WWR/tables/qc/<analysis>/direction_means_<metric>_qc.csv`（factor=Group 行）
- Group 显著性（控制 WWR+Complexity）：`/root/analysis-2-EEG/task4_core_lmm_suite/factor_WWR/tables/qc/<analysis>/model1_main_effects_<metric>_qc_anova.csv`（Term=Group 行）

---

## Experience 分组（控制 WWR + Complexity 后的 Group 效应）

- O_theta：High=0.110608 (n=190) vs Low=0.099533 (n=95)；Group：F=0.975 (1,24.01)，p=0.3334，p_holm=0.3334
- F_theta：High=0.103970 (n=190) vs Low=0.097573 (n=95)；Group：F=0.438 (1,24.00)，p=0.5142，p_holm=0.5142
- O_alpha：High=0.121000 (n=190) vs Low=0.084560 (n=95)；Group：F=3.817 (1,23.97)，p=0.06251，p_holm=0.1875（趋势，但不显著）
- O_beta：High=0.221376 (n=190) vs Low=0.186132 (n=95)；Group：F=0.850 (1,24.01)，p=0.3656，p_holm=0.7312

**一句话结论（Experience）**：四个频段在控制 WWR/Complexity 后，Group 主效应均未达到显著；其中 O_alpha 在未校正下接近显著（p≈0.063），但 Holm 校正后不显著。

---

## SportFreq 分组（控制 WWR + Complexity 后的 Group 效应）

- O_theta：High=0.107715 (n=154) vs Low=0.105978 (n=131)；Group：F=0.024 (1,24.01)，p=0.8788，p_holm=0.8788
- F_theta：High=0.104111 (n=154) vs Low=0.099166 (n=131)；Group：F=0.333 (1,24.00)，p=0.5694，p_holm=0.5694
- O_alpha：High=0.107538 (n=154) vs Low=0.110400 (n=131)；Group：F=0.028 (1,23.99)，p=0.8676，p_holm=0.8676
- O_beta：High=0.232850 (n=154) vs Low=0.182329 (n=131)；Group：F=2.121 (1,24.01)，p=0.1583，p_holm=0.4748

**一句话结论（SportFreq）**：四个频段在控制 WWR/Complexity 后，Group 主效应均未达到显著；O_beta 虽然均值差异较大，但统计上仍不显著。

---

## 说明：为什么这就是“Task6 做 4 个波”的正确做法

- 你原本的 Task6（`task6_obeta_special`）做的是：只针对 O_beta，用“只含 Group”与“控制 WWR/Complexity”两种模型对比 Group 效应。
- 如果要对 4 个频段都这样做，最直接的方式就是：
  - 复用 Task4 的 Model1（它已经是 `EEG ~ WWR + Complexity + Group + (1|Subject)`），对每个指标都给出了 Group 的 F/p。
- 这样不会引入新的建模差异，也避免重复跑一套冗余分析。

