# BD36 | PtSe₂ 双层 Layer Edelstein Effect 研究计划（PLAN.md, rev R1.1）

## 1. 目标与科学基线

检验真实 PtSe₂ 双层中面内 DC 电流 j∥ 是否诱导层分辨自旋密度 **S⁻z = S_top^z − S_bottom^z**（Layer Edelstein Effect, LEE），及其滑移(registry)依赖性。明确区分：

- **真实电子结构结果**（本计划要做的）：PBE+D3+SOC 第一性原理 + 真实层/原子投影的 Wannier 自旋算符 + 数值响应。
- **非结果**：BiDB 数据库标量（gap_pbe 2253=0.219 / 2255=0.614 eV）、任何 `hr.dat` 推断的层投影 spin、BD34 人造 Fourier 向量场。

## 2. 结构集（READ_ONLY 审计已通过，见 STRUCTURE_SYMMETRY_AUDIT.md）

| 阶段 | 结构 | 对称性 | 角色 |
|---|---|---|---|
| Round 1 | 2253 双层（未滑移）| D3d | **高对称零响应基准**（i+C3 禁止 S±z）|
| Round 1 | 2253+u=(+0.08,+0.19) 刚性滑移 | Ci | **一般滑移**（i 允许 S⁻z，C3 缺失）|
| Round 1 | 2252 单层 | D3d | 独立物理对照（无层间耦合，方法/基组标定）|
| Deferred | 2253+u=(−0.08,−0.19) 反向滑移 | Ci | 符号检验（S⁻z(u) = −S⁻z(−u) 预期）|
| Deferred | 2255（registry 双层）| D3d | 第二高对称基准（gap 0.614 eV）|
| Deferred | 2254（异取向双层）| C3v | 独立取向族，不并入滑移族分析 |

## 3. 计算协议（零场非磁 SOC；逐级有界；R1.1 修正）

每结构逐级独立、有界、幂等：

1. **registry-preserving 弛豫**（非磁、PBE+D3、IBRION=2/ISIF=2/NSW≤60、EDIFFG=−0.01）：
   - 双层：**选择性动力学 x,y 固定、z 自由（F F T）**（VASP 语义：T=可动/自由、F=固定）→ 面内分数滑移精确不变，只弛豫层距/翘曲；单层全弛豫 T T T；
   - 单层 2252：全弛豫（T T T，独立物理对照）；
   - 弛豫后**必须**跑 `scripts/verify_relaxed.py`：复核反演/C3/镜面、滑移残差（≤1e-3）、层距/晶胞变化（ISIF=2 晶胞不变），并记录弛豫前后层隙。
2. **SOC 静态**（LNONCOLLINEAR=T, LSORBIT=T, SAXIS=0 0 1, 微小初始磁矩 0.02 z, PBE+D3, ENCUT=500, k 9×9×1 基组 → 收敛表 7/9/11/13）→ 总能量、能带、净磁矩（非磁 SOC 应为 0）、费米能位置。
3. **SOC+W90 接口**（vasp_ncl + wannier90.win 文件，ISYM=−1 全 k 网格）→ `.eig/.amn/.mmn/.spn`、Wannier centers、插值能带 → 响应计算。

**VASP-W90 接口流程（R1.1 修正，4 步）**：
1. 静态 SOC 收敛（ICHARG=1 续算）后，用 `scripts/make_win_kpoints.py IBZKPT|OUTCAR` 从**实际 VASP k 列表**生成 win 的 `begin kpoints`（ISYM=−1 时 9×9×1 全 BZ = 81 点；禁止手写 k 列表）；
2. `wannier90.x -pp` 生成 `.nnkp`；
3. VASP 接口运行：`LWANNIER90=.TRUE.` + `LWRITE_SPN=.TRUE.`（.spn 由 VASP 写出；win 中**不设置** `spin=true`，INCAR 中**不嵌入** `WANNIER90_WIN`——均为已标识的不兼容/易错配置）；
4. `wannier90.x` 后处理 → `hr.dat`、centers、插值检验。

**费米能位置**：2253 电荷中性有 ~0.22 eV 隙（BiDB 数据库值）。DC 响应需导电态 → 响应阶段以电子/空穴参考掺杂明确指定费米能靠近导带/价带；σ 奇异/禁带时**不得**换算 χ^j。

## 4. Wannier 层分辨自旋算符方案（R1.1 修正）

- **投影基组统一为 Pt:d + Se:p**（无 s）：双层 num_wann = 2×(5×2 + 3×4) = **44 spinor**；单层 2×(5+6) = **22 spinor**。num_bands=64（44+20 缓冲），**最终由投影能带检验决定**。
- **spin 矩阵来源（两级区分，见 §5）**：`LWRITE_SPN=.TRUE.` 由 VASP 写出 `seedname.spn`（44×44×3 Pauli 矩阵元）；win 不设 `spin=true`。`.spn` 解析与校验（维度/Hermite/规范一致性）用 `scripts/read_spn.py`。
- **层投影**：Wannier centers（`wannier_centres.xyz` / wout z 坐标）按层间隙中点分顶/底层 → 层算符 S_layer^z = P_layer σ_z P_layer，P_layer = Σ_{n∈layer}|w_n⟩⟨w_n|；S⁻z = S_top − S_bottom，S⁺z = S_top + S_bottom。
- **必须检验**：hr.dat 平移/厄米/规范、.spn 一致性、WF 层投影（center z 分布与原子层吻合）、带插值误差（≤5 meV 目标，记录实测）。

## 5. 响应设计（总自旋 vs 层分辨自旋；PT 简并规范安全）

**两个算符必须区分，不能混用：**
- **总自旋算符** S_i = <w_n|σ_i|w_m>（`.spn`，Wannier 基 Pauli 矩阵元）——描述 WFs 的自旋结构，本身**不分辨层**；
- **层分辨自旋算符** S_layer^z = P_layer S_z P_layer（P_layer 由 Wannier centers 按层投影）——LEE 观测量 S⁻z = S_top − S_bottom、S⁺z = S_top + S_bottom 都由它构造。
- 原子投影（按原子轨道）与层投影（按 WF center z）**不同**：本任务用**层投影**（WF center 判层），原子投影仅作交叉检验。

**PT 简并（规范安全）处理：**
- SOC 下中心对称结构 (iT)²=−1 → 每个 k 点 Kramers 二重简并 → 逐带 spin 期望值无定义。
- 响应计算采用**简并子空间投影迹**（gauge-invariant，不选带规范）：
  - 在 Bloch 基同规范下 Fourier 变换 H(k)、S(k)（含层投影）；
  - 每 k 点：对 H(k) 对角化，取 PT 简并子空间（对能量简并对），
    χ(S, E) = (e/2N_k) Σ_k Tr_degen[ P_k S(k) P_k (∂H/∂k 投影顶点) ] 按 Kubo 静态极限实现；
  - **必须验证的判据**：SOC-off 对照（S±z 响应 ≡ 0）、S⁺z=0（i 禁戒）、2253 零基准（C3 禁戒）、滑移 (+u) 非零允许、反向滑移符号、层间耦合敏感性、k/展宽收敛、零点无关性。

## 6. 数值验证门槛（每结构必过）

1. SOC-off 对照：无 SOC 时 S±z 响应必须为 0。
2. S⁺z=0 检验：反演对称结构 S⁺z 响应 = 0（i 禁戒）。
3. 零基准：2253（D3d）S⁻z 响应 = 0（C3 禁戒）。
4. 滑移族：(+u) S⁻z 允许非零；反向滑移 S⁻z(u)=−S⁻z(−u)（Deferred）。
5. 层间耦合：2252 单层对照；层距敏感性（可选）。
6. 数值收敛：k 7/9/11/13、σ 展宽、零点无关性。
7. 验收：对称性“允许”≠“算得非零”；只报已执行事实（RUN_FINISHED / NUMERIC_CONVERGED / FE_VERIFIED 分列）。

## 7. 静态检查与自动化（R1.1 新增）

- `scripts/static_check_inputs.py`：每结构 17 项 + Round-1 预算检查（投影数/num_wann、k 点数、POSCAR↔win 结构一致性、INCAR 标签兼容性、资源上限、无内部绝对路径、无 POTCAR 入库）；0 FAIL 才可进入授权队列。
- `scripts/verify_relaxed.py`：弛豫后 CONTCAR 复核（反演/C3/镜面/滑移残差/层距）。
- `scripts/make_win_kpoints.py`：从实际 VASP k 列表生成 win kpoints（接口一致性）。
- `scripts/read_spn.py`：.spn 解析/维度/Hermite/规范校验（输出后运行）。

## 8. 交付物（本分支 Draft PR #99）

- `STRUCTURE_SYMMETRY_AUDIT.md`（✓ R1 完成，R1.1 复核通过）
- `INPUTS/`（✓ R1.1 重建：SD 标志/投影/接口/win 原子坐标/k 点策略全部修正，冻结哈希）
- `COST_AND_PERMISSION_GATE.md`（资源与授权清单，含 IMR 核查结果与首个可执行方案）
- `PLAN.md`（本文档）、`LOGS/`（审计/构建/静态检查日志）
- 后续（授权后）：数值、原始证据、失败诊断、`RESULTS/`

## 9. 执行顺序与授权门

1. R1.1 修正（本提交）→ 用户审阅。
2. 用户批准**具体平台与作业资源**（见 COST_AND_PERMISSION_GATE.md）后：平台就地验证作业（SOC+W90+spin 最小接口测试）→ Round-1 作业 → 监控 → 数值验收 → 结果归档。
3. Deferred 阶段（反向滑移、2255、2254）在 Round-1 验收后单独批准。
