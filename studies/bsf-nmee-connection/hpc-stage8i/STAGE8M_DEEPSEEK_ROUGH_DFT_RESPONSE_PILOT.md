# Stage8m | DeepSeek 自主执行任务：真实 WSe2 粗粒度 DFT 响应与 idea/novelty GO-NO-GO

**用户目标（优先级）：** 不要把项目困在无限的解析器/单带自旋检查中。用户明确认为，用已有真实 WSe2 DFT 数据尽快判断 idea 是否有物理信号、是否值得发展，是当前重要目标。优先得到一个诚实、可证伪的**粗粒度材料响应估计**，再决定是否投入更密网格、完整 PAW layer operator 和生产级 Wannier。

**研究仓库/分支：** `Xin-Jiaqi/quantum-geometry-reproductions`, `research/bsf-nmee-stage7-sw-gate-20261020`。所有新增结果只属于独立研究 `studies/bsf-nmee-connection/`。Jia et al. PRL 原论文复现必须保持冻结，Draft PR #3 不合并。

## A. 物理问题与不许越界的声明

我们当前研究的是一个明确命名的 **0 K equilibrium/static spin-only-Zeeman derivative of the mixed spin-electric Berry curvature**，不是自动等同于 Jia et al. 的 driven-DC NMEE、也不是 current-driven Edelstein accumulation。

\[
\Gamma_{iz;j}^{\rm stat}=\left.\partial_{B_z}\int_{\rm BZ}
\frac{d^2k}{(2\pi)^2}\,i\,\mathrm{Tr}\,
P[\partial_{\lambda_i}P,\partial_{k_j}P]\right|_{\lambda=B=0},
\]
其中 \(H(k,\lambda,B)=H_0(k)+\lambda_iS_i-g\mu_B B_zS_z\)，\(\lambda_i\) 是与所定义自旋算符耦合的能量源。响应的绝对单位必须在 spin normalization、electron charge、二维面积与 reciprocal basis convention 确认后才给出。

已知的真实结构事实：
- AB：真实 Stage8h CONTCAR，C3v，total in-plane static tensor 有一个独立 scalar 系数，**允许非零但未计算**。
- 半滑移 midpoint：真实 z-only relaxed 几何保留 layer-exchanging glide；total static tensor **严格为零**，layer-odd tensor 只被允许，不意味着一定非零。
- 两个原 SOC SCF 都有可用 EIGENVAL/WAVECAR；旧 PBS 退出 255、PROCAR 截断，因此需只读重后处理。禁止把 21×21 sampled positive gap 称为 rigorously converged continuous-BZ gap。

## B. 先做 Stage8l 的零新作业核验（必须真正执行，不可只更新代码）

1. Fetch 最新分支、记录 HEAD SHA；运行 `stage8l_crosscheck_band_ledger.py`。保留原 Stage8k JSON，生成版本化 corrected evidence JSON：七个非简并新记录标为 band35 的实际对应 band36；简并 glide 样本不做任意单带匹配。
2. 在 HPC 原目录只读使用 Stage8i 已存在 WAVECAR/EIGENVAL，SHA 必须等于：
   - AB WAVECAR `1d1a5dbf858503d9b9d07a85088bac8c3088af44a1598c0472724473739919c3`
   - AB EIGENVAL `2112a13fac886fb5f6bdc72ce0b0549d08c7ccf47c49b7382c5239e932cd4618`
   - glide WAVECAR `750b6dd5806c69bb17180242349bc40fc5cb3e86ef6ea24c1434d2fe47058127`
   - glide EIGENVAL `b70c4a88873833d80e82c1a8b76327716d5d2db2f7b9b6f15a55b7e1a5e89c94`
3. 使用 `stage8l_pseudo_interband_gate.py` 分别对 AB/glide 做只读 screening，传入上面 expected SHA。先核验每个目标 k 精确落在网格、WAVECAR 与 EIGENVAL 的 k 与全部带能量一致、band 33–36 (1-based) 真正占据且 37–44 真正空。任何一项失败就 stop，不改容差去“通过”。
4. 记录各选定 k 的 32 个 \(\langle v|S_{x,y}|c\rangle_{\rm pseudo}\) 与占据/空带完整 block Frobenius norm。个态 pair 在简并子空间内 gauge-dependent，优先比较完整 block norm；明确这是 Euclidean-normalized pseudo spinor gate，**不是物理 PAW 自旋、速度矩阵元或 NMEE**。
5. 在零作业结果中给出一个明确判断：是否存在非简并跨隙 pseudo spin-flip channel，AB 与 glide 是否有实质差异，哪些目标 k 是真正网格点。不能仅凭它宣布 idea 成立或新颖。

## C. 真正向材料响应迈进：现有 WAVECAR 的低成本 VASP2WANNIER90 后处理

Stage8i 原 INCAR 已使用 `SAXIS=0 0 1`、`NBANDS=96`、`LORBIT=11`、`ISYM=0`、21×21×1，原 WAVECAR/EIGENVAL 的 SHA 已冻结。VASP 6.4.3 按官方文档支持 `LWRITE_SPN`（版本 >6.4.2），其默认 SAXIS 约束与本项目原始 `SAXIS=0 0 1` 一致。官方说明：
- `ALGO=None` 可从已有 WAVECAR 做 post-processing，不做 orbital update；`NELM=1`、`LWAVE=.FALSE.`、`LCHARG=.FALSE.` 用于避免写回原文件。
- `LWANNIER90=.TRUE.` 输出 `wannier90.mmn`, `.amn`, `.eig` 等；`LWRITE_SPN=.TRUE.` 在非共线计算中输出 `wannier90.spn` spin-matrix elements；`LWRITE_MMN_AMN=.TRUE.` 显式要求 overlap/projection files。

参考官方文档：
- https://vasp.at/wiki/ALGO
- https://vasp.at/wiki/LWANNIER90
- https://vasp.at/wiki/LWRITE_MMN_AMN
- https://vasp.at/wiki/LWRITE_SPN
- https://vasp.at/wiki/LORBIT

### 作业预算与提交边界

**用户已要求继续做粗粒度 DFT idea-feasibility 验证。** 允许的最高预算是 **最多 2 个 post-processing-only PBS 作业**，分别 AB 与 glide，必须使用原 WAVECAR，不进行 SCF 迭代、离子弛豫或新 k-mesh 计算。若站点环境/输入兼容不确定、VASP2WANNIER90 编译能力未证实、输入 SHA 不符或 dry-run 失败，停止提交，保留自动恢复状态，给出唯一具体阻断项。不得转而重做 SCF、NEB 或增加 k-mesh 作业。

对每个结构：
1. 重新确认原 WAVECAR/CHGCAR/POSCAR/POTCAR/KPOINTS、VASP 二进制版本、POTCAR 许可、`NBANDS=96` 与 441 个全 BZ k 点；所有源文件只读。新建唯一 postprocess 目录，不能覆盖原目录或原 WAVECAR。
2. 在独立目录创建基于原输入的后处理 INCAR，目标设置为 `ISTART=1`, `ALGO=None`, `NELM=1`, `LWAVE=.FALSE.`, `LCHARG=.FALSE.`, `NSW=0`, `LWANNIER90=.TRUE.`, `LWRITE_SPN=.TRUE.`, `LWRITE_MMN_AMN=.TRUE.`, `LORBIT=11`。保留原 `ENCUT`, `NBANDS`, `LSORBIT`, `LNONCOLLINEAR`, `SAXIS`, `ISYM`, `KPOINTS`, `POSCAR`, 合法本地 POTCAR；按 VASP 6.4.3 实际手册和现场输出确定 `ICHARG`，不要盲目覆盖默认。设置/确认 `LWANNIER90_RUN=.FALSE.`，本阶段只输出矩阵数据，不运行 Wannier minimization。
3. **先在副本上 dry-run/小型格式预检。** 检查 `vasp_ncl` 确有 VASP2WANNIER90 + SPN 支持。让 VASP 生成/修补 `.win` 时，验证 `num_bands=96`、`num_wann` 与实际 matrix dimensions、441 k points、spinor 标志、`spn_formatted=true`、`mp_grid=21 21 1`、KPOINTS 顺序完全一致。不得把 `num_wann=96` 当作已知事实强行设置；只有确认 full-band Bloch basis 文件可写出才使用该设置。若 AMN 投影配置阻止生成 MMN/SPN，停止并诊断，不用随机 projections 硬凑。
4. 只在预检通过且原输入 SHA/许可/站点队列有效时，分别提交最多 1 个 AB 和 1 个 glide 的 postprocess-only 作业。它们不允许更新波函数/电荷：检查输出前后 WAVECAR SHA 相同；记录 PBS ID、walltime、退出码、stdout 结束标志与完整输出文件 sizes/SHA。
5. 如果矩阵文件成功产生，原始大文件留在 HPC，不推送 Git。只提交 sanitised manifest、输入/输出 SHA、解析检查、简短数值和可复现分析代码。

## D. 粗粒度真实材料响应估计（本阶段主要科学产出）

目标不是再做一份 spin texture，而是尽量用真实 DFT Bloch-space data 做一次 **粗网格、明确标记为 exploratory** 的 \(\Gamma^{\rm stat}\) 估计，决定是否值得投入生产级 Wannier。优先考虑不做 MLWF localization 的直接 Bloch-basis Wilson-plaquette 方法：

1. 用 VASP 输出的 `.eig`/EIGENVAL 建立每个 k 的对角 \(H_0(k)\)，用 `.spn` 的完整 band-off-diagonal \(S_x(k),S_z(k)\) 矩阵构造 \(H(k,\lambda,B)=H_0+\lambda S_x-g\mu_BBS_z\)。先核对 SPN 的矩阵规范与 normalization（Pauli \(\sigma\) 还是 \(\sigma/2\)），将对角 spin matrix elements 与 PROCAR/Stage8k 可靠的非简并方向做 consistency check；不要把 `SPN` 文件名字当作 PAW 完备性的证明。
2. 从 `.mmn` 读取全 96 带相邻 k 的 Bloch-overlap matrices \(M(k,k')\)，必须处理 reciprocal-boundary G vectors，不能把所有 links 当作普通相邻数组。用 perturbation eigenvectors \(U_{\rm occ}(k,\lambda,B)\) 构造 covariant links：
   \[
   L(k,k')=U_{\rm occ}(k,\lambda,B)^\dagger M(k,k')U_{\rm occ}(k',\lambda,B),
   \]
   λ links 则用相同 k 点原始 Bloch basis 中的 \(U_{\rm occ}(\lambda)^\dagger U_{\rm occ}(\lambda+\delta\lambda)\)。对占据整个 36-band manifold 使用 determinant Wilson plaquette，避免单条带规范依赖。每个 link 用 SVD/polar-unitary 修正并记录最小奇异值，不能静默忽略接近奇异的 links。
3. 先写小型合成哈密顿量 regression：同一个代码必须重现现有 Stage7i 的独立 Kubo/Wilson static \(\Gamma\) 基准（见 `docs/STAGE7I_GEOMETRY_DFT_BRIDGE.md` 和 `bsf_nmee/mixed_wilson_geometry.py`）。若没有通过这个回归，禁止运行/引用材料响应结果。
4. 对每个 B 值先算 \(\mathcal F_{\lambda_x k_x}\) 和 \(\mathcal F_{\lambda_y k_y}\) 的 plaquette phase，再用 \(B=\pm0.125,\pm0.25\) T 的中心差分估 \(\partial_B\mathcal F\)。源步长先试 \(\delta\lambda=0.05,0.10,0.20\) meV；所有 λ/B corners 都检查 36/37 gap 仍开、occupied projector dimension=36。必须报告 source-step stability 与 gap guards，不挑能产生漂亮信号的步长。
5. 用原 21×21×1 完整 BZ 先得到 **coarse DFT response estimate**，以实际 POSCAR 晶格把 fractional reciprocal derivatives 转换为 Cartesian \(k_x,k_y\)，再按二维 reciprocal-cell area 计算 BZ 权重。利用现有 96-band 矩阵做 nested remote-band-window check（occupied 36 固定，空带截断分别到 48/64/80/96）；Stage7g 已显示远带可强烈抵消，若结果随 band cutoff 明显变化则不得称稳健。报告全部 2×2 tensor entries、C3v isotropy/off-diagonal symmetry residual、有限差分符号、原始离散 BZ 求和与归一化约定。绝不能把 21×21 当作 k-mesh convergence 或发表级数值。
6. Glide midpoint total static response 是精确几何对称性 null test；数值总响应应接近零。定义 ratio = abs(Gamma_glide_total) / max(abs(Gamma_AB_total), noise_floor)，报告是否与数值噪声相容。不要从总自旋 `SPN` 推导 layer-odd response。
7. 在没有可靠 PAW operator convention、完整 BZ mesh convergence 或同一物理 Hilbert-space basis links 前，数值状态只能叫 `COARSE_BLOCH_SPACE_STATIC_RESPONSE_ESTIMATE`。它是判断 idea 是否有意思的早期证据，不是 production prediction。 该估计是在固定零场 Kohn–Sham Hamiltonian中加入 spin-only 源的 frozen-KS 响应，不包含外场诱导的自洽电荷/交换关联重整化，也不含 orbital magnetic coupling。

### 不要为了“非零”而放松科学门槛

- AB 结果必须对 λ/B 步长稳定，且 tensor 满足 AB 的 C3v 对称约束；否则判为数值不可信/需诊断。
- glide midpoint 的 total response 应是 symmetry zero；如果显著非零，首先修正 reciprocal-boundary links、k-point ordering、spin source normalization、symmetry/geometry，不可直接解释成物理效应。
- 若 AB coarse response 小于数值噪声或随步长反号/剧烈变化，停止大规模 Wannier；优先判断是否是远带抵消、矩阵规范或 mesh 粗糙造成。只有稳定且显著的结果才进入下一门槛。
- 一个有限的 \(\langle v|S_x|c\rangle\) 或有限的 spin-flip norm 不等于有限 \(\Gamma\)，更不等于新颖性。
- 当前 novelty prior art 已由 `docs/STAGE10_PEER_LANDSCAPE_20261021.md` 梳理。要讨论创新，必须比较“静态平衡响应”与 LEE、FE Rashba、spin-QGT、层奇 NLH、已有 insulating NMEE。不要把“第一个算出某个材料的数值”自动包装成新效应。
- 真正独特且尚待验证的假说是：**真实 3R-WSe2 滑移路径中，glide 对称性强制总静态响应为零，却允许可由物理 layer projector 定义的 layer-odd 响应有限；AB 端点的 total scalar static response 可作为对照。** Stage8m 先评估 AB total 数值是否值得继续；layer-odd midpoint 是随后必须单独做的物理算符门槛。

## E. 必须输出的交付物

1. 一页中文报告：band index erratum、cross-gap pseudo spin gate、是否成功输出真实 `MMN/SPN/EIG`、矩阵 normalization 与 PAW status、AB 粗响应、glide total null-test、误差/收敛、novelty 判断和明确 STOP/GO。
2. 新的 versioned sanitised JSON，包含 input/output SHA256、VASP version、PBS IDs（如果有）、完整参数、matrix dimensions、mesh、finite-source-step series、gap/occupied-subspace guards、2×2 response tensor、symmetry residual、单位约定、known limitations。
3. 新的 response script + synthetic Stage7i regression tests；测试必须验证规范变换、TR、tensor symmetry、source-step guard、glide null test和错误 MMN link 的 fail-closed 行为。
4. 更新 `RESEARCH_INDEX.md` 和导航性 paper research README。保持 PR #3 Draft。
5. 最后一句明确回答：**粗粒度真实材料数据是否显示足够稳定/显著的响应，值得再投入 layer-odd operator + dense Wannier？** 不能以“还需要更多计算”作为唯一结论。

## 明确禁止

- 不得提交/上传 POTCAR、WAVECAR、CHGCAR、PROCAR、MMN、SPN、AMN、私人 OUTCAR 或任何凭据到 GitHub。
- 本阶段最多两次 `ALGO=None` postprocessing-only VASP 作业；**零新 SCF、零 NEB、零新增 k-mesh SOC、零磁场自洽有限差分、零生产级 Wannier 任务**。
- 不改历史 Stage8i/8k JSON，不改论文原复现。
- 如矩阵规范/输出文件/站点编译能力不确定，停止并报告具体阻断，不编造 response，也不重跑 SCF。
