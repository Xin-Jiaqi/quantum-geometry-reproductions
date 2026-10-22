# AgF₂ B1 初筛后的磁群与张量选择规则：零 VASP 科学审查

**2026-10-22｜SYMMETRY_CONDITIONALLY_PASS；REAL_MATERIAL_IGMC_NOT_COMPUTED。** 本轮无新增 VASP/PBS 作业，也没有访问私有原始 OUTCAR。

## 1. B1 初筛状态

来源：PR #38 commit 247af15662c160b32e142cb2ecc0ae982aa7af3f 的 IDEA_SCREENING_VERDICT.md、RUN_RECEIPT.json、TRAJECTORY.csv 与 POSCAR.screen_last_complete。11 个离子步全部电子自洽；VASP 报告 EDIFFG=-0.08 eV/Å 达到，末态 max|F|=0.0782 eV/Å，TOTEN=-35.42794489 eV/Ag4F8 12原子，最近跨层 F–F=3.002 Å；Ag 磁矩 [-0.567,+0.567,+0.566,-0.566] μB，net M≈0；初判 gap 1.3875 eV (12×8×1 未收敛)。计算 8 ranks × 1.62 h，无重复正式 VASP，单项 SCREENING_PASS。**这不是已证明的局部极小值、磁基态、磁易轴、滑移双势阱或 IGMC。**

## 2. 用真实 screening 坐标独立核验结构操作与磁矩配对

直接读取公开的 12 原子最终 POSCAR，正交晶胞 a=3.8270249701 Å，b=5.8796657458 Å，c=26.93854144 Å，Ag 1–4，F 5–12。该 POSCAR 坐标原点下的一项结构 glide 是

\[
g_x=\{m_x|(\tfrac12,\tfrac12,0)\}:\qquad
(x,y,z)\mapsto(-x+a/2,\ y+b/2,\ z).
\]

按二维 PBC 最近像做同元素一一站点匹配（真空 z 不跨胞），**12/12 正确，最大匹配残差仅 1.28×10⁻⁵ Å**。原子对为 Ag 1↔2、3↔4；F 5↔8、6↔7、9↔12、10↔11。相符于历史来源预检 LG12 pb11 (3D 等价 Pc #7)。附加独立可复算纯 Python 脚本 [check_magnetic_symmetry_screen.py](check_magnetic_symmetry_screen.py)。

对于 Ag 共线标量磁矩 [- + + -]，每对站点磁矩反号。**如果**采用来源 PRL/此前 magnetic preflight 的物理 Néel x 方向假设，实自旋是轴矢量，单位 m_x 对 Sx 不变，而 T 反号；于是 Ag 的反号配对符合内部**反幺正 m'_x**，MPG 候选 m'_x={E,m'_x}。

**重要条件：**此次计算是 ISPIN=2、LSORBIT=.FALSE. 的纯共线无 SOC 算例。VASP 的自旋量子化坐标不能被视为已测得物理 x 易轴。无 SOC 的自旋可整体旋转；只有选定 x 方向以后才可以把上述结果用作真实 SOC 材料磁点群的预判。SOC 若产生不同易轴、canting 或其它对称性破缺，必须重新评估。

## 3. 单态 α 的全部允许分量

定义低频绝缘体静态平衡磁电极化率

\[
\alpha_{ij}^{PB}=\left.\frac{\partial P_i}{\partial B_j}\right|_{B=0}.
\]

P 是时间偶极矢量，B 是时间奇轴矢量；磁操作 (R,T^θ) 约束

\[
\alpha'=(-1)^\theta\det(R)\,R\alpha R^\mathsf T.
\]

对于内部 m'_x，R=diag(-1,+1,+1)，θ=1；因此允许张量是

\[
\boxed{\alpha=
\begin{pmatrix}
\alpha_{xx}&0&0\\
0&\alpha_{yy}&\alpha_{yz}\\
0&\alpha_{zy}&\alpha_{zz}
\end{pmatrix}}.
\]

**αxx、αyy 都是允许的**，但允许性并不证明其非零。值得强调测量方向：**αzx=∂Pz/∂Bx 禁止**；若研究面外 Pz 的线性磁电变化，只能考虑该群下的 αzy（By）或 αzz（Bz）。单态自发 Pz 在 m'_x 下允许，但**还未验证实际非零极化**。

## 4. Type-I / Type-II 跨 B1→B2 磁连接

根据已有 Zhu 2025 来源与源构造预检，Type-I 层间 Néel 反平行的域连接为幺正 m_z；Type-II 层间 Néel 平行、来源对应 4% 双轴应变的域连接为反幺正 m'_z。R_z=diag(1,1,-1)。相同外参和坐标系下：

\[
\alpha_2^{\mathrm{Type-I}}=-R_z\alpha_1R_z^\mathsf T,\quad
\alpha_2^{\mathrm{Type-II}}=+R_z\alpha_1R_z^\mathsf T.
\]

| Quantity | m'_x 单态 | Type-I m_z 域间 | Type-II m'_z 域间 |
| --- | --- | --- | --- |
| αxx=∂Px/∂Bx | 允许 | **反号** | 同号 |
| αyy=∂Py/∂By | 允许 | **反号** | 同号 |
| αzz=∂Pz/∂Bz | 允许 | **反号** | 同号 |
| αyz=∂Py/∂Bz | 允许 | 同号 | **反号** |
| αzy=∂Pz/∂By | 允许 | 同号 | **反号** |
| αzx=∂Pz/∂Bx | 禁止 | 恒零（在该群下） | 恒零（在该群下） |
| 自发 Pz | 允许 | **反号** | **反号** |
| AHE σxy (T-odd axial z) | 允许 | 同号 | **反号** |

**关键澄清：** Type-I Pz 反号**不代表** αzy 必须反号，后者实际上 Type-I 同号。不要把 FE switching、magnetic-order switching、AHE switching 与所有 α 分量的切换混为同一件事。Type-I/II 是**不同磁连接**，不是把本次 0% 无 SOC 构型自动视为两种磁基态。

把 screened B1 原子坐标经 m_z + 合适平面平移逐一映射即可免费构造**SYMMETRY_IMAGE B2**，连同自旋的轴矢量/时间反演规则一起赋值。对无外场、具备此晶体镜面对称的 Hamiltonian，这是严格对称相关结构，域间能量相同是对称性结果，**不是独立 B2 DFT 证据，更不能证明滑移铁电双势阱或势垒。** B1 此时本身也只有 0.08 eV/Å 筛选精度。

## 5. 最重要的物理风险：非 SOC 共线 spin-only 的响应可严格为零

当前 B1 粗 DFT 完全**没有 SOC**。在一个更明确的、限定条件的电子模型中：

- 固定的共线 AFM 自洽势，且无 SOC、无 spin-flip hopping；
- uniform-g、均匀 spin Zeeman \(H_B=-\mu\mathbf B\cdot\boldsymbol\sigma\)，B 足够弱，不改变能带占据；
- 只研究**clamped magnetic texture 的 single-particle spin-Zeeman** 的绝缘体 polarization response；电极化算符为 spin scalar，不包含磁矩自洽旋转、磁致伸缩、orbital-B 或空间依赖 g。

取磁序 n 为自旋量子化轴，有 \(H_0(k)=H_\uparrow(k)\oplus H_\downarrow(k)\)。

1. \(B\parallel n\) 时，Zeeman 项分别是 spin-block 内比例于单位阵的均匀能量平移，occupied states 不变，所以 \(\partial P/\partial B_\parallel=0\)。
2. \(B\perp n\) 时，Zeeman 一阶扰动只混合两自旋块，而 polar charge/position operator 自旋对角，因此其线性期望迹为零，\(\partial P/\partial B_\perp=0\)。

**结论仅适用于上述具体 rigid, collinear, SOC-free spin-only 模型：所有 α 的 spin-Zeeman 静态线性分量为零，即使 m'_x 群论允许它们。** 不能推广为“所有无 SOC magnetoelectricity 消失”，因为 magnetostructural feedback、非共线 magnetic order、自旋相关 g、轨道耦合等超出此模型。八带 toy model 有 spin-dependent interlayer hopping，也不能不检查自旋对易性便与无 SOC AgF₂ 等同。Xiang et al. PRL 134, 116301 (2025) 的机制针对 **spin-orbit coupled materials**。

**这改变下一个试算的优先级：** 再运行多个纯无 SOC 的磁畴单点，对论证 nonzero spin-Zeeman IGMC 的帮助很有限。若准备再花少量 WSL 算力，优先考虑**一个 B1、SOC-on 的轻量单点**，确认所选 x 磁矩方向下能隙与内部 MPG 是否仍然合理；然后再讨论具体响应核/ Wannier。这一步即使通过也不能证明 α 非零，只有必要性信息。

## 6. 下一阶段、明确决策和边界

- **G1 当前通过**：screened B1 结构的实际 glide m_x 与条件磁态 m'_x 对应；αxx/yy/zz 等在该候选磁群下对称允许；磁连接符号如表（**SYMMETRY_CONDITIONALLY_PASS**）。
- **G2 待检验**：含 SOC 的单态磁易轴、SOC 磁群、绝缘 gap，及 nonzero spin-Zeeman mixed k–B response（全未计算）。
- **G3 待检验**：真正的 B1/B2 极化双稳态、滑移势垒和 Pz（未计算，不能由对称像跳过）。

**零 VASP 结论：AgF₂ 值得继续，但下一次的高信息收益并非机械补 B2，而是先明确 SOC 后是否存在目标自旋混合与允许响应。任何新 VASP 仍需独立预算授权。**

### 资料与验证路径

1. 实际公开结构：PR #38 / stage2_real_dft/b1_wsl_relax/idea_screen_1/POSCAR.screen_last_complete 与 RUN_RECEIPT.json。
2. 历史磁对称性来源：同项目 stage6_harness/stage1_review/STRUCTURES_AND_METHODS.md，包含逐原子源对称检查及原论文磁连接定义。
3. 响应机制：Xiang et al., *Intrinsic Gyrotropic Magnetic Current from Zeeman Quantum Geometry*, PRL 134, 116301 (2025), DOI https://doi.org/10.1103/PhysRevLett.134.116301。
4. 本目录 check_magnetic_symmetry_screen.py 负责站点一一映射、磁轴假设验证和 α 张量符号表。
