# BD36 | FAST IDEA GATE：先验证 LEE，再考虑 Wannier （2026-10-23）

> 目的：避免为一个尚无真实材料正信号的 idea 反复修理完整 DFT→Wannier→.spn 工具链。本文是**新的验证计划**，不表示任何新 DFT 已经执行，也不扩大用户先前授权。

## 已知事实与风险

- BiDB 2253 及平行滑移结构位点审计已通过：原构型 D3d、一般刚性滑移 Ci，全族可保持反演。**一般滑移允许**层反向电流致面外自旋密度 `S_-^z`，**不是已经算到非零**。
- BD34 的六个零点是人为 Fourier registry 场引入、不是 PtSe2 的真实材料独立预测。
- para 真实环境预检 **FAIL**：当前账户既无 Pt/Se POTCAR，又无可用 .spn 输出的 VASP 构建。**para Phase-0 当前路径停止，不再为旧 VASP 补丁做工程修复，不再尝试提交 2253 基准的旧任务。**该已批准的一次 para Phase-0 从未使用，且不能迁移成对其他平台/构型的授权。
- 主问题：真实电子结构中，`j_parallel -> (S_top^z - S_bottom^z)` 是否存在稳定、非忽略不计的响应？优先回答材料是否有候选正信号，再决定参数、零点与材料库筛选。

## 四级关卡（严格按信息增益顺序，失败就停）

| Gate | 最小任务 | GO 条件 | NO-GO / 升级 |
|---|---|---|---|
| **G0 环境，无计算** | 只读调查 IMR VASP 6.4.3+ / PT+Se 合法 POTCAR / LORBIT、`PROCAR` 与 SOC；或 WSL VASP6.6+ 合法 Pt POTCAR。对每个平台最多一轮确认；不要修 para VASP 5.4.4。| 找到现成可用的自旋轨道 DFT 执行环境 | 没有就停 DFT，转向 **公开 QE fully relativistic Pt/Se UPF** 的非生产性可行性研究或现有材料数据；不擅自装/编译/提交作业 |
| **G1 真实自旋纹理候选** | **先做一组一般滑移** BiDB 2253+`u=(0.08,0.19)`：确认间隙与固定 x,y、只弛豫 z 的几何可行性，再做非磁 SOC 有限 k 点，使用 `LORBIT=11` 输出 `PROCAR`，无需 Wannier/`.spn`。初始 9×9×1 用来查电子结构/投影及 k 纹理。| 正确识别上下层原子、SOC doublets 与非零的 Kramers-**双重态迹** `s_-^z(k)`；控制 `s_+^z(k)≈0` 和 `s_-^z(-k)≈-s_-^z(k)` | 如果数据仅有大投影误差/不可分双重态，诚实标为 **INCONCLUSIVE**，不宣称材料无效；先检查局域投影定义 |
| **G2 粗粒度电流响应** | 在候选有信号的费米能附近电子/空穴参考掺杂，增密 k 网格（例如 21×21×1 → 33×33×1），用双重态能量有限差分速度、`PROCAR` 上下层 m_z 双重态求和，近似常弛豫时间计算 `Q` 与 Drude `D`，仅导电且 `D` 非奇异才取 `F=QD^{-1}`。加入 2253 高对称 C3 零响应对照。| 物理上允许的 `F_-^z` 在两个 k 网格/展宽窗口下稳定，且高对称对照接近数值零；报告有量纲与投影缺失 | 若接近噪声，**停该候选**；不能推广否定理论或拓扑。此阶段结果为快速材料筛查，非论文精确响应 |
| **G3 精细计算** | **仅 G2 GO 后**才做 SOC Wannier+spin/layer 算符校准、非阿贝尔 PT 规范、足够密 k 网格、多 registry 零点搜索以及真实响应定量。 | 获得可复核电子结构与层分辨响应 | G3 工作另行申请具体授权；不自动开展 |

## 关键物理量与验收（防止把 spin texture 当成 Edelstein）

给任意滑移构型及每个 `k` 的完整 PT Kramers doublet `D_{νk}`，先用 VASP `PROCAR` SOC 输出的原子投影 m_z 定义*近似*层自旋纹理：

`s_{-,ν}^z(k) = Σ_{n∈D_{νk}} [Σ_{atom∈top} m^z_{atom,nk} − Σ_{atom∈bottom} m^z_{atom,nk}]`.

- **必须以完整双重态求和**，不允许用随机规范的单条简并能带的自旋作为物理结果。双重态投影迹在局域算符投影固定时是规范不变的。
- PT 预言总层自旋迹 `s_+^z(k)≈0`；T 预言 `s_-^z(k) = −s_-^z(−k)`（非磁结构、匹配双重态）。在孤立高对称点，允许 k-resolved 纹理非零，但 C3 对宏观 `j→S_-^z` 的张量积分为零，不能误把单个 k 非零当成 Edelstein 效应。
- 若采用常 `τ` Boltzmann 近似，`Q_a ∝ Σ_{νk} v_{νa}(k) s^z_{-,ν}(k) [-f'(ε_{νk}−μ)]`；`D_ab ∝ Σ_{νk} 2 v_{νa}v_{νb}[-f']`；`F=QD^{-1}`，使 `τ` 在比值中消去。需声明投影算符、单位、掺杂、温度/展宽、速度差分方法、PT 分组容差、完整 BZ 求和。
- VASP 官方文档确认 `LNONCOLLINEAR` + `LORBIT=11` 的 `PROCAR` 含各离子的三个磁矩方向投影，但**它是离子局域投影近似**，不能等同全空间层分辨自旋或 Kubo 结果。源：https://vasp.at/wiki/index.php/PROCAR 及 https://vasp.at/wiki/LORBIT。
- 如果没有合适 VASP PAW 许可路径，Quantum ESPRESSO PS Library 有 Pt、Se 的 **PBE fully relativistic UPF**（公开），但 **QE `projwfc.x` 在 SOC 下输出 `j,m_j` 投影，不可直接当 `S_z`**；必须另行验证 spinor wavefunction → layer-spin 的正确接口。源：
  https://pseudopotentials.quantumespresso.org/legacy_tables/ps-library/pt
  https://pseudopotentials.quantumespresso.org/legacy_tables/ps-library/se
  https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html

## 本阶段应保留的定理与未证之事

- 定理/选择定则：中心对称平行滑移 PtSe2 家族允许层反向 `S_-^z` 但禁止层同向 `S_+^z`；C3 固定构型 `F_-^z=0`。前述至少一组三点额外零点，仅在 `G/M/N` 是非退化指数 +1 零点、其余零点孤立、全滑移响应场光滑周期等假设成立时才适用。
- 未证：真实 PtSe2 一般滑移下 `F_-≠0`、数量级、符号、额外零点位置、材料可测性与相对已有滑移 REE 工作的新颖性。

## 执行政策（消除用户中转）

1. 本 Gate 的**研究组织与 READ-ONLY 环境探测**可由研究代理按现有只读范围继续，所有进展直接写入 **PR #99**；不用用户转发 AI 之间的审查意见。
2. **没有获得新的、有界的计算授权**前，禁止在 IMR、WSL、QE、para 提交任何生产作业。先将平台候选、具体作业与最小核时提请批准，不要沿用已失效的 para 2253 Phase-0 授权。
3. 一份 gate 表足够：`G0_PLATFORM_AVAILABLE / G1_SPIN_TEXTURE / G2_TRANSPORT / G3_QUANTITATIVE` = PASS/FAIL/INCONCLUSIVE/NOT_RUN。不要继续生成 R1.4/R1.5/R1.6 多轮全局检查文档。
4. 候选 2253+`u=(.08,.19)` 一旦 G2 NO-GO，暂停该材料方向，转向其他真实物理材料候选，而不是试图通过调节预设 Fourier 模型再次得到需要的零点。
