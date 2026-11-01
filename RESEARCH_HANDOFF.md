# 核心物理接入：逐篇交接入口

更新：2026-11-01。当前任务只整理并接通已有五篇论文的核心计算，暂不开展新选题。一次只推进一篇；原文争议与可用计算分开记录，不为消除 `PARTIAL` 扩展成完备材料或散射理论项目。

## 优先级与停止条件

| 顺序 | 论文 | 影响核心调用的缺口 | 本篇做到哪里即停 |
|---|---|---|---|
| 1，已接通 | [NMEE](papers/jia-2026-nonlin-magnetoelectric-edelstein/README.md) | PT简并/多带完整E×B计算藏在底层，缺统一的物理输入与验证入口 | 接通明确频率协议、两种作用顺序、热占据和同规范算符；独立求解/收敛通过；DC争议留明，不强行冒充作者Fig1c |
| 2，已接通 | [Spin Transport](papers/xiang-2025-spin-quantum-geometry/README.md) | 约化Eq13与完整多带响应容易混用；原SM积分和多源项有分歧 | 确认已有完整响应内核的适用范围，补实际缺失通道；以独立密度矩阵验证；给出可改模型的入口，不将材料DFT设为模型使用门槛 |
| 3，排队 | [Magnetotransport](papers/xiang-2026-quantum-geometric-map/README.md) | spin/orbital以及bulk/device能力混杂；原SM符号、三带轨道项需分开 | 优先闭合论文核心模型的指定spin/orbital通道，标出一般多带仍缺的具体项；停止扩大到普适Kubo–NEGF或全自洽无序 |
| 4，排队 | [Curvature–Flux](papers/cui-2026-curvature-flux-zeeman-qgt/README.md) | v3的频率符号、Table I及PR107无效反例影响可靠调用 | 修正反例，按真实反幺正群作用检查目标响应；保持几何/频率入口；不新增EPW、浴或材料设施 |
| 5，排队 | [IGMC](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/README.md) | 核心模型已可用；依赖导入、绝对幅度和反幺正约定需收尾 | 修复实际运行阻碍；完成指定模型单位/频率桥；作者原始数据缺失不无限追索 |

排序依据是核心调用受阻程度，不是论文重要性。每篇交付：一个清晰入口、一个可重算例子、独立核验、明确边界。下一篇从其当前代码及正文/SM重新核对，不能只接受历史验收文字。

## NMEE 本轮结果

- `nmee_lab.evaluate_dynamic`：接通SM (6)/(7)的完整有限频率E/B响应；支持精确简并、有限温度费米面项、一般多带、非零位置连接及k依赖自旋算符。
- `cumnas_dynamic_model`：原Eq5模型及A/B子晶格可观测量；明确物理k与晶格长度。
- 两种作用顺序单独输出；复数频率之和、Pauli单位、cosine输入的1/4因子明确。有限展宽不称为微观散射时间，结果不称为已提取的本征DC。
- 独立轨道空间N²求解与k有限差分、CuMnAs双重简并投影、六带热占据、规范变换和负对照通过。38项模块测试（原30+新增8）。
- 固定AC参数下，481²至721²网格相对变化5.8311e-7；分块计算避免密网格内存膨胀。原五组冻结数据不变。
- 源文：SM p3 Eq9对角项按p2的Bbar=-g mu_B B收缩后，应得到与Eq11印刷项相反的符号。这是已独立确认的局部内部不一致，不等于已得到正确的全篇DC公式。

详见该模块的[当前验证记录](papers/jia-2026-nonlin-magnetoelectric-edelstein/docs/FINAL_VERIFICATION.md)和[原文核查](papers/jia-2026-nonlin-magnetoelectric-edelstein/docs/AUTHOR_DISCUSSION_BRIEF.md)。

## Spin 本轮结果

完整一阶/二阶电场响应接通，spin current与spin polarization使用物理算符；SM四源项、规范变换及非零四带积分检验通过。76项测试在Python3.11/3.12通过。紧凑Eq13/14的作者省略项与全有限频率响应分开；材料DFT与物理tau0仍不混称完成。下一篇Magnetotransport。

## 继续工作规则

1. 当前工作为本地 `feature/five-paper-core-access` 分支，起点 `86980c87890f0ec674db37a659a88e131762f77f`；先检查实际分支、HEAD与未提交改动。不要覆盖他人工作，也不要把本地结果声称为远端CI通过。
2. 一篇处理完后再进入下一篇，只更新这个交接入口和相应论文文档，避免新增重复的版本报告。
3. 核对H、物理算符、位置连接、频率与eta、占据、单位和简并处理；不能从HR文件自行假定物理spin或position。
4. 原文逐式与独立算法冲突，分成：确切内部代数错误、协议/约定差异、我们的实现问题、尚缺原始数据。不要用测试数量代替科学证据。
5. 沿用已有冻结值作为回归证据；新增物理例子另存清楚命名的结果。源码清单因改动而更新时，记录旧清单哈希与起点提交，不能改旧数值来凑通过。
6. 只有用户重新开启选题讨论后再做idea/studies。FDT与交替磁三阶Hall尚未接入，本轮不增加第六/七篇。

## 关于根目录的统一接口

`integration/physics_interface/`仍是四模块的元数据/设计契约，不是数值调度器。调用当前NMEE应使用其模块入口；不要把根 `contracts.evaluate` 的占位接口当作已实现。其他论文完成核心接入后，再决定是否需要一个很薄的统一入口；现在不重做平台架构。
