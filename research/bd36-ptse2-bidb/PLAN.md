# BD36 | PtSe₂ 双层 Layer Edelstein Effect 研究计划（PLAN.md）

## 1. 目标与科学基线

检验真实 PtSe₂ 双层中面内 DC 电流 j∥ 是否诱导层分辨自旋密度 **S⁻z = S_top^z − S_bottom^z**（Layer Edelstein Effect, LEE），及其滑移(registry)依赖性。明确区分：

- **真实电子结构结果**（本计划要做的）：PBE+D3+SOC 第一性原理 + 真实层/原子投影的 Wannier 自旋算符 + 数值响应。
- **非结果**：BiDB 数据库标量（gap_pbe 2253=0.219 / 2255=0.614 eV）、任何 `hr.dat` 推断的层投影 spin、BD34 人造 Fourier 向量场。

## 2. 结构集（READ_ONLY 审计已通过，见 STRUCTURE_SYMMETRY_AUDIT.md）

| 阶段 | 结构 | 对称性 | 角色 |
|---|---|---|---|
| Round 1 | 2253 双层（未滑移）| D3d | **高对称零响应基准**（i+C3 禁止 S±z）|
| Round 1 | 2253+u=(+0.08,+0.19) 刚性滑移 | Ci | **一般滑移**（i 允许 S⁻z，C3 缺失）|
| Round 1 | 2252 单层 | D3d | 独立物理对照（无层间耦合，S⁻z 无定义，用于基组/方法标定）|
| Deferred | 2253+u=(−0.08,−0.19) 反向滑移 | Ci | 符号检验（S⁻z(u) = −S⁻z(−u) 预期）|
| Deferred | 2255（registry 双层）| D3d | 第二高对称基准（gap 0.614 eV）|
| Deferred | 2254（异取向双层）| C3v | 独立取向族，不并入滑移族分析 |

所有滑移为**刚性**（NOT relaxed）；计算先做“固定滑移的有限结构优化”。

## 3. 计算协议（零场非磁 SOC；逐级有界）

每结构三级（均为独立、有界、幂等作业）：

1. **registry-preserving 弛豫**（非磁、PBE+D3、IBRION=2/ISIF=2/NSW≤60、选择性动力学 x,y 冻结 z 自由）：保证面内分数滑移精确不变，只弛豫层距/翘曲。晶胞取 BiDB 源胞（真空 18.7 Å > 15 Å）。EDIFFG=−0.01 eV/Å。
2. **SOC 静态**（LNONCOLLINEAR=T, LSORBIT=T, SAXIS=0 0 1, 微小初始磁矩 0.02 z, PBE+D3, ENCUT=500, Gamma 9×9×1 基组 → 收敛表扩展 7/9/11/13）→ 总能量、能带、磁矩（非磁 SOC 应无净磁矩）、费米能位置。
3. **SOC+W90 接口**（vasp_ncl + 嵌入 wannier90.win，模式同已校验的 para HgI2 L1 benchmark）：产生 hr.dat、`.spn`、Wannier centers、插值能带 → 响应计算。

**费米能位置**：2253 在电荷中性有 ~0.22 eV 隙（BiDB 数据库值）。DC 响应需导电态 → 在响应阶段以参考掺杂（电子/空穴）明确指定费米能靠近导带/价带，**不做**“费米能落在禁带内”的假导电处理。σ 奇异/禁带时不得换算 χ^j。

## 4. 基组与 Wannier 层分辨自旋算符方案

- **投影候选**：Pt 5d + Se 3p：双层 22 空间轨道 → **44 spinor**（num_wann=44, spinors=true）；单层 11 → 22。初始 num_bands=64（44+20 缓冲），**最终由投影能带检验决定**（投影态局域度、froz 窗口内带序、插值误差）。
- **spin 矩阵**：wannier90 `spin = true` → `wannier90.spn`（<w_n|σ_i|w_m>，dimensionless Pauli），逐项校验：矩阵维度 44×44、Hermite、σ_z 与 hr 同规范；若 .spn 不可用，使用有记录的等价矩阵（UNK 基自旋期望 + 明确文档）。
- **层投影**：Wannier centers（`wannier_centres.xyz` / wout 中 z 坐标）按层间隙中点分顶/底层；层算符
  S_top^z = Σ_{n∈top}|w_n⟩⟨w_n| σ_z 的 Wannier 矩阵 → S⁻z = S_top − S_bottom，S⁺z = S_top + S_bottom。
  **必须检验**：hr.dat 平移/厄米/规范、.spn 一致性、WF 层投影（中心 z 分布与原子层吻合）、带插值误差（≤ ~1-5 meV 目标，记录实测）。
- **响应**：Bloch 基同规范 Fourier 变换 S(k)、H(k)；**PT 简并子空间投影迹**（或等价已验证的 gauge-covariant Kubo）计算 χ^E(S⁻z)、χ^E(S⁺z) 与纵向 σ；仅在非奇异 σ 导电态换算 χ^j = χ^E·σ⁻¹。

## 5. 数值验证门槛（每结构必过）

1. **SOC-off 对照**：无 SOC 时自旋响应必须为 0（含层分辨）→ 检测算符/流程错误。
2. **S⁺z=0 检验**：所有反演对称结构 S⁺z 响应 = 0（i 禁戒）→ 数值验证选择定则。
3. **零基准**：2253（D3d）S⁻z 响应 = 0（C3 禁戒）→ 高对称零响应基准。
4. **滑移族**：(+u) 滑移 S⁻z 允许非零；反向滑移 S⁻z(u)=−S⁻z(−u)（符号检验，Deferred 阶段）。
5. **层间耦合敏感性**：单层 2252 对照（无第二层，S⁻z 退化检查）；改变层距的敏感性测试（可选）。
6. **数值收敛**：k 网格 7/9/11/13、展宽 σ、偶极修正（不需要：P1 镜像族无净偶极；记录在案）、零点（算符位置无关性检验）。
7. **验收**：不得把“对称性允许”说成“算得非零”；只报告已执行事实（RUN_FINISHED ≠ NUMERIC_CONVERGED ≠ FE_VERIFIED 分列）。

## 6. 交付物（本分支 Draft PR）

- `STRUCTURE_SYMMETRY_AUDIT.md`（✓ 已完成）
- `INPUTS/`（✓ 已生成：6 结构 × relax/soc/soc_w90/KPOINTS/POSCAR_SD/w90.win + PBS/Slurm 草稿，冻结哈希）
- `COST_AND_PERMISSION_GATE.md`（资源与授权清单）
- 运行日志：`LOGS/`（审计输出、构建输出、未来计算日志）
- 后续（授权后）：数值、原始证据、失败诊断、`RESULTS/`

## 7. 执行顺序与授权门

1. 本计划 + 审计 + 输入 → 用户审阅（Draft PR）。
2. 用户批准**具体平台与作业资源**（见 COST_AND_PERMISSION_GATE.md）后：Round 1 作业（2253、+u 滑移、2252 单层）→ 监控 → 数值验收 → 结果分支归档。
3. Deferred 阶段（反向滑移、2255、2254）在 Round 1 验收后单独批准。
