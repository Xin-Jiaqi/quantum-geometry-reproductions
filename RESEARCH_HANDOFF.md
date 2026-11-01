# 核心物理接入：逐篇交接入口

更新：2026-11-01。当前任务只整理并接通已有五篇论文的核心计算，暂不开展新选题。一次只推进一篇；原文争议与可用计算分开记录，不为消除 `PARTIAL` 扩展成完备材料或散射理论项目。

## 优先级与停止条件

| 顺序 | 论文 | 影响核心调用的缺口 | 本篇做到哪里即停 |
|---|---|---|---|
| 1，已接通 | [NMEE](papers/jia-2026-nonlin-magnetoelectric-edelstein/README.md) | PT简并/多带完整E×B计算藏在底层，缺统一的物理输入与验证入口 | 接通明确频率协议、两种作用顺序、热占据和同规范算符；独立求解/收敛通过；DC争议留明，不强行冒充作者Fig1c |
| 2，已接通 | [Spin Transport](papers/xiang-2025-spin-quantum-geometry/README.md) | 约化Eq13与完整多带响应容易混用；原SM积分和多源项有分歧 | 确认已有完整响应内核的适用范围，补实际缺失通道；以独立密度矩阵验证；给出可改模型的入口，不将材料DFT设为模型使用门槛 |
| 3，已接通 | [Magnetotransport](papers/xiang-2026-quantum-geometric-map/README.md) | spin/orbital以及bulk/device能力混杂；原SM符号、三带轨道项需分开 | 优先闭合论文核心模型的指定spin/orbital通道，标出一般多带仍缺的具体项；停止扩大到普适Kubo–NEGF或全自洽无序 |
| 4，已接通 | [Curvature–Flux](papers/cui-2026-curvature-flux-zeeman-qgt/README.md) | v3的频率符号、Table I及PR107无效反例影响可靠调用 | 修正反例，按真实反幺正群作用检查目标响应；保持几何/频率入口；不新增EPW、浴或材料设施 |
| 5，已接通 | [IGMC](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/CORE_ACCESS.md) | 核心模型已可用；依赖导入、绝对幅度和反幺正约定需收尾 | 修复实际运行阻碍；完成指定模型单位/频率桥；作者原始数据缺失不无限追索 |

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

完整一阶/二阶电场响应接通，spin current与spin polarization使用物理算符；SM四源项、规范变换及非零四带积分检验通过。76项测试在Python3.11/3.12通过。紧凑Eq13/14的作者省略项与全有限频率响应分开；材料DFT与物理tau0仍不混称完成。Magnetotransport随后已收拢，见下。

## Magnetotransport 本轮结果

统一示例直接调用spin PHE/MNHE、orbital PHE/MNHE与通用几何。修正旧v4入口对Eq56的误读：d_ab f包含fsecond链式项，两行并不冲突；历史数值另标。Eq49→52因子2与Eq20→31相对符号问题继续单列。一般三带轨道扩展不是原文两带近似的验收前提。SI测试按物理常数版本比较，源求解器/旧数值不改。1214项全量检查中仅一项文件指纹待更新，保留旧指纹后该项复测通过；新增入口结果已记录。

## Curvature–Flux 本轮结果

几何/通量、完整带间GHE/KME可调用；补原文σ/2的KME低频入口。PDF p4/p12/p13核对：Eq16与Eq112动态项一致；Eq19与Eq126一阶项确有内部冲突，轨道Liouville支持附录。无效m′和C3径向反例已用物理不变性测试排除。375项通过，四带BZ 64²→96²相对变化6.86e-8。Table I全表与真实金属散射不设为本篇模型门槛。

## IGMC 本轮结果

Rashba/Dirac的C、D与完整非共振响应已有可靠内核，新增一处重算入口，列明Hz/rad·s⁻¹和SI换算。费米海/费米面以及波包/格林函数/轨道Liouville一致。SM Eq34/46局部系数问题重读确认，但不推断作者Fig2(d)程序。纠正旧Round8说明：完整1 THz与低频结果的微小差别是有限频率项，不是101²网格误差。Python3.12显式载入完整本地toolkit后76项全部通过；干净Python3.11环境75通过、1可选外部toolkit跳过。Round1–8冻结核验和Round9本地完整性门槛通过；历史“LOCAL_DRAFT”是旧验收器标签，不代表本轮远端状态。

## 五篇当前使用边界

**核心物理已在指定模型/响应协议下接入，可以支持后续idea检验；没有把五篇全部实验/材料/图像标为完全复现。** 优先直接调用各模块入口，改变模型或协议后只验证相应假设，不预先扩大为通用材料平台。

| 论文 | 本轮直接重算入口 | 本轮无需追平的内容 |
|---|---|---|
| NMEE | 模块`nmee_lab.evaluate_dynamic`及README示例 | 作者Fig1c本征DC、完整散射极限 |
| Spin | `scripts/reproduce_core.py` | 材料DFT、将有限η自动认作tau0 |
| Magnetotransport | `scripts/reproduce_core.py` | 超出原文省略三带过程的完整轨道理论、统一bulk/device散射 |
| Curvature–Flux | `scripts/reproduce_core.py` | Table I全表、完整金属带内/耗散KME、真实材料SI |
| IGMC | `scripts/reproduce_core.py`及CORE_ACCESS.md | 作者Fig2(d)原始幅度来源、Eq14完整约定、一般多带共振/散射 |

本轮论文PDF正文及SM/附录均只读取于本地Zotero；没有复制进代码库。源文争议与我们误读的更正已分别记录。停止增加新物理设施，等待用户开启idea讨论。

## 继续工作规则

1. 当前工作为本地 `feature/five-paper-core-access` 分支，起点 `86980c87890f0ec674db37a659a88e131762f77f`；先检查实际分支、HEAD与未提交改动。不要覆盖他人工作，也不要把本地结果声称为远端CI通过。
2. 一篇处理完后再进入下一篇，只更新这个交接入口和相应论文文档，避免新增重复的版本报告。
3. 核对H、物理算符、位置连接、频率与eta、占据、单位和简并处理；不能从HR文件自行假定物理spin或position。
4. 原文逐式与独立算法冲突，分成：确切内部代数错误、协议/约定差异、我们的实现问题、尚缺原始数据。不要用测试数量代替科学证据。
5. 沿用已有冻结值作为回归证据；新增物理例子另存清楚命名的结果。源码清单因改动而更新时，记录旧清单哈希与起点提交，不能改旧数值来凑通过。
6. 只有用户重新开启选题讨论后再做idea/studies。FDT与交替磁三阶Hall尚未接入，本轮不增加第六/七篇。

## 关于根目录的统一接口

`integration/physics_interface/`仍是四模块的元数据/设计契约，不是数值调度器。调用当前NMEE应使用其模块入口；不要把根 `contracts.evaluate` 的占位接口当作已实现。五篇已有各自数值入口，目前不增加统一调度层，也不重做平台架构。
