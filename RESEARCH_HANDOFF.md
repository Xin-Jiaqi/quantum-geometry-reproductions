# RESEARCH HANDOFF — Quantum Geometry Reproductions

**给新的 ChatGPT / Codex / 研究者的交接入口。** 本文为 2026-10-21 的快照，不是实时的科学验收结论。每次开始新任务，务必重新读取 live GitHub 的 README.md、RELEASE_STATUS.md，以及所选论文目录中的 README 和最新 acceptance 文档；其他研究线程可能同时更新其余论文。

## 研究任务

本仓库的目标是逐篇独立复现量子几何与输运论文的**主要模型、可验证数值方法和物理能力**，不是在存在原文公式/单位分歧时人为凑出全部发表图像。现在要把可靠的各论文计算方法用于自主科学研究，尤其是群论、磁性/层交换对称性、spin quantum geometry、Edelstein 与非线性自旋响应，之后再考虑有真实微观跃迁的滑移、堆叠、插层调控。

四篇论文已上传各自有独立测试的模型/方法子集，但**四篇科学验收状态均为 PARTIAL**。软件 CI PASS 不等于所有公式、实验协议、真实材料与绝对量级已被物理复现。

## 四篇独立论文，不能互换物理协议

| Registry ID | 论文模块 | 已有可用子集 | 保持开放 |
| --- | --- | --- | --- |
| jia2026-nmee | [Jia PRL 2026 NMEE](papers/jia-2026-nonlin-magnetoelectric-edelstein/) | Kane–Mele、CuMnAs、PT 投影、printed/static 分支、Wannier Stage-1 只读预检 | 驱动有限温金属 DC 不等于 equilibrium static；HR-only 无完整物理算符 |
| xiang2025-igmc | [Xiang PRL 2025 IGMC](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/) | Rashba/Dirac Zeeman geometry、磁场驱动 AC 电流、Kubo、磁点群 | Fig.2(d) 绝对量级和时间反演/频率约定 |
| xiang2025-spin-qgt | [Xiang/Jin/Wang PRL 2025 Spin Quantum Geometry](papers/xiang-2025-spin-quantum-geometry/) | Dirac spin-QGT/NDSC/reduced Eq.13，多带有限 eta Liouville，PT 双重简并，Laurent，SM Eq.2 四源项 | Fig.1(a) 约 2 倍、SM Eq.55/56 和 Eq.69/70，真正金属 intrinsic tau^0 |
| xiang2026-qgmap | [Xiang PRB 2026 magnetotransport](papers/xiang-2026-quantum-geometric-map/) | Dirac/Rashba、E/B 双频 Liouville、部分 spin/orbital 响应、SCBA/Kubo/Ward | 物理 steady-state EB/PHE、全 conserving Kubo–Streda、实材料量级 |

IGMC 的 Zeeman geometry/磁场驱动 AC、电场驱动自旋量子几何、NMEE 和 2026 PRB 的 E/B 混合磁输运是不同问题。不得把四套 solver 的响应公式、物理算符或 eta 约定混成一个未定义的 QGT。

## 新会话按这 7 步继续

1. 从 GitHub 读取本文件、[接口规范](docs/PHYSICS_INTERFACE_SPEC.md)、[并行规则](docs/CONCURRENT_RESEARCH_WORKFLOW.md) 和 integration/physics_interface/protocol_registry.json。记录使用的当前 quantum-geometry-reproductions **完整 Git SHA**，不要只记录 main。
2. 选择一个 **paper stable ID + protocol ID**，阅读该论文最新 README、实际计算入口、验收和未解决问题。若用户提供原文 PDF + Supplemental Material，应先查证使用的公式、单位和符号；GitHub 不包括没有公开再分发许可的 PDF/SM，没有材料就不要声称读过。
3. 在论文独立 Python 环境安装并运行其自有回归脚本，保存具体参数和 stdout/结果 JSON。先做到对原计算结果的只读复核，不修改冻结基准。
4. 明确 H(k) 基底、速度/位置/Berry 连接、自旋与层算符、k 的 Cartesian/fractional、能量与自旋单位、占据、温度、驱动 E/B 频率、简并、eta 的物理含义和响应张量定义。缺失关键算符就 BLOCKED。
5. 对新物理仅在 studies/ 下创建独立目录和 Git 分支。固定纸面代码 SHA、Group Theory Toolkit SHA、物理参数/单位及独立测试数据，不在 papers/ 里偷改已有公式。
6. 新结果用同一 H(k) 的原文 baseline + 独立算法 + 对称性、简并/gauge 与收敛核对；出现不一致建立 OPEN_DISCREPANCY，保留 printed 与 derived 两条分支。
7. 对其他模块的任何代码修改和 push/merge 都要独立 PR、用户授权、作用域审查与目标模块 CI；不要覆盖他人并行修改。

## 当前交付的统一接口是什么

已完成的是**版本化元数据/物理协议发现接口**、Python 设计契约、只读 CLI 和独立 CI。这一层不直接替代四套已有计算代码，也**尚未实现统一数值 evaluate 调度器**，后者只能在明确物理协议并通过每篇论文的交叉验证后逐步提供 adapter。

从仓库根目录运行：

~~~bash
python -m integration.physics_interface.cli list
python -m integration.physics_interface.cli inspect xiang2025-spin-qgt
python -m integration.physics_interface.cli plan xiang2025-spin-qgt spinqgt.finite_eta
python -m integration.physics_interface.cli validate
python -m unittest discover -s integration/physics_interface/tests -v
~~~

其中 list/inspect/plan/validate **不会执行论文物理计算**，plan 只返回人工检查的命令和协议。真实物理复现还需选中论文目录的 verify 和 tests。每次接口迭代后先跑测试，不能以 registry 的 PASS 代替 physics 的 PASS。

## 给新 ChatGPT 的复制提示

> 请先查看本仓库根目录 RESEARCH_HANDOFF.md、docs/PHYSICS_INTERFACE_SPEC.md、integration/physics_interface/protocol_registry.json、当前 RELEASE_STATUS.md；再审阅我提供的论文 PDF/SM，核对选中论文的源公式及当前 Git SHA。请选择具体响应协议，不要把 eta→0 的有限常数直接当物理 intrinsic DC，也不要从 wannier90_hr.dat 单独推算 spin/Zeeman 响应。先复现仓库冻结计算，然后仅在 studies/ 里发展新模型，独立检查对称性、算符、单位、简并、收敛与基准。如果需要 GitHub 写操作，先明确变更文件和授权。

## 下一步科研工作

先针对一个已经严格验证的 TB 模型建立**只读 adapter**（保留原 paper solver），并做独立数值对照；再调用 [Group Theory Operations Toolkit](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit) 已有的磁群、层群、层交换、张量基计算。现有工具库 v0.18.0 的已知提交 SHA 为 551a0daf8d4ba5a5f24b2811694bfe438616de48，后续每个 study 要重新固定其实际 SHA。二阶**自旋流**张量的磁反幺正与频率变换必须另行推导，不能生搬普通一阶 Edelstein。

之后的滑移、插层模型要保证物理相对位移确实改变原子轨道间跃迁 T(k,delta)，并满足对称性协变、层交换、周期性、能带与 k 网格收敛。一个可规约的整体相因子不构成真实滑移物理。真实材料量级需要同规范 spin/position/velocity 等投影算符以及必要的 DFT/Wannier、结构能量或电荷数据。

## 对并行开发的明确约束

本轮不改任何 papers/ 模块、已有 studies/ 或它们的测试与 CI；只增加根目录本文件、两份 docs/ 架构说明、integration/physics_interface/ 和自己的 scoped CI。共同的 README.md 和 RELEASE_STATUS.md 此轮也**不动**，避免与其他并行工作竞争。后续若对现有论文能力/入口有新的修改，以 live GitHub 为准，单独 PR 更新 registry。

最终原则：**软件验收和科学验收独立；公式分歧保留；可复用物理算符与计算协议要可追溯；新物理放 studies/。**
