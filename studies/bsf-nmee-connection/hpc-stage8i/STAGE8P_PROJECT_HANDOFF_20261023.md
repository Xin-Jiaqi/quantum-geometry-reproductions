# BSF × NMEE 原创研究交接文档 | Stage8o-real → Stage8p / 下一轮真实材料响应

**日期：**2026-10-23  
**仓库：** [Xin-Jiaqi/quantum-geometry-reproductions](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions)  
**研究分支：** \`research/bsf-nmee-stage7-sw-gate-20261020\`  
**对应 Draft PR：** [#3](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/3)，**不得自行合并**  
**本 handoff 编写前最后已验证代码 HEAD：** \`60da6b6bb8ef96eee8f8fd131320160a6501da08\`；对应两条 GitHub CI **SUCCESS**；BSF 工作流 **194 tests passed**。发布本 handoff 自身会产生更新的 commit；以后请以远端 HEAD 为准。

> 本文件旨在让新 ChatGPT/DeepSeek harness **直接接手科研工作，不让用户重新解释背景或复制粘贴多轮命令**。不可将历史来源报告当作经过独立验收的真实物理响应。必要时引用完整源码及机器证据。日本 IMR 原始受许可 DFT 数据不进入公开 GitHub。

## 0. 一句话科学状态

**真实 3R-WSe₂ 的 AB 与 glide-half-slide 已有 z-only 几何、SOC SCF、以及 EIG/MMN/SPN 后处理的真实 DFT 数据；Stage8m 静态 NMEE 估算未通过 AB C₃v 与 glide total-null 验收。Stage8n 曾错误归因于 VASP ISYM=0 数据对称性，Stage8o-real 用真实晶体 Seitz 操作和 441/441 k 配对否定了该“排他性归因”的原始证据。独立 Kubo/Wilson 非零材料结果仍未交叉通过。idea 尚开放，不启动生产级计算。**

严格状态：
\`MATERIAL_RESPONSE_SYMMETRY_GATE_FAILED__VASP_ROOT_CAUSE_NOT_ESTABLISHED__NONZERO_KUBO_MATERIAL_CROSSCHECK_PENDING__IDEA_OPEN\`

### 项目三条主线

1. **Idea：** stacking/sliding symmetry + intrinsic equilibrium/static magnetic-source derivative of mixed spin–electric quantum geometry，关注总自旋响应和物理定义的 layer-staggered 响应在 glide 中点的不同选择定则。
2. **真实材料验证：**先粗 DFT 判断信号与能否独立复现；重点物理证据包括 PAW 自旋/速度、混合 Berry 响应、占据－空带、远带抵消和数值规范协变。
3. **同行工作：**原研究 [Stage10 peer landscape](../docs/STAGE10_PEER_LANDSCAPE_20261021.md) 已评估已有 LEE、LNHE、ferroelectric Rashba、spin QGT 与 insulating NMEE。**“sliding FE 会切换响应”或“有非零自旋混合”本身不足以宣称创新**。论文新意尚未确立。

## 1. 物理定义与硬性防混淆

研究中的冻结 KS 平衡/静态几何核（并非 Jia 等 PRL 的 driven-DC 协议）：
\[
\Gamma_{iz;j}^{\mathrm{stat}} =
\left.\partial_{B_z}\int_{\mathrm{BZ}}\frac{d^2 k}{(2\pi)^2}\,
i\,\mathrm{Tr}\!\left[
P[\partial_{\lambda_i}P,\partial_{k_j}P]
\right]\right|_{\lambda_i=B_z=0}.
\]

**约定需统一检查：**若 VASP \`SPN\` 输出的是 Pauli 矩阵 \(\sigma_\alpha\)，物理自旋源需明确 \(S_\alpha=\sigma_\alpha/2\)，Zeeman 也须同时使用 \(H_B=-g\mu_B B_z S_z\)。不同子程序 \(-\mu_B B\sigma_z\) 与 \(-2\mu_B B\sigma_z\)，或 \(\lambda\sigma_x\) 与 \((\lambda/2)\sigma_x\) **不可直接比较**。数值总核按目前微分定义的未加电荷、二维厚度等因子的量纲是 \(\mathrm{\AA^{-1}\,eV^{-1}\,T^{-1}}\)，旧报告所写 \(1/(\mathrm{eV\,T\,\AA^2})\) 不成立。

- **AB C₃v：**对 \(M_{ij}=\Gamma^{\mathrm{stat}}_{iz;j}\) 的 in-plane **Cartesian** tensor，应满足 \(M=\gamma I_2\)，因此 \(M_{xx}=M_{yy}\)、\(M_{xy}=M_{yx}=0\)。具体诊断前检查斜倒格矢、源指标和 Jacobian 定义。
- **真实半滑移 nonsymmorphic \(\{m_z|\mathbf t\}\)：**总（uniform）静态 \(M\) 严格为零；**layer-odd** 响应可以在对称性上允许，但从未计算其非零数值。glide 对 in-plane \(k\) 为恒等，但对 axial \(S_x,S_y\) 变号；真正电子检验需 sewing/带窗口封闭。
- **严格 \(S_z\) 守恒**意味着当前横向响应 no-go；Stage8k 的非简并**伪自旋子**非零支持 AB 中 \([H,S_z]\neq0\)（spin-diagonal PAW 广义重叠下的定性判断），**不证明**真实材料 NMEE 非零。
- 不得把普通 driven-current Edelstein、equilibrium-static \(\Gamma\)、论文源 PRL driven-DC NMEE、layer-staggered 响应混为一个物理量。
- 原论文复现 \`papers/jia-2026-nonlin-magnetoelectric-edelstein/\` **PARTIAL、冻结**。该目录的 \`research/README.md\` **仅做导航链接**，不向原复现混入独立研究代码。

## 2. 已经完成的真实 DFT 和重要证据

| 阶段 | 实际工作 | 目前可引用的结论 |
|---|---|---|
| Stage8a–f | 3R-WSe₂ 主候选、3R-MoS₂ 对照、结构/群论路径筛选 | AB C₃v、midpoint layer-exchanging glide 是具体材料可检验结构；其允许/禁止规则不证明非零效应 |
| Stage8g–h | 日本 IMR VASP 6.4.3 PBS 真实 3 点无 SOC z-only relax（s=0,0.25,0.5） | 受限能量 0/+13.633/+24.117 meV/胞；中点 glide 保留。**不是 NEB/切换势垒** |
| Stage8i | 真实两次 SOC SCF（AB 491186、glide 491187） | 21×21 采样间接 gap 1.088/1.149 eV；PBS exit 255，PROCAR 截断，EIGENVAL/WAVECAR 可读；**不是全 BZ 连续 gap 保证** |
| Stage8j | 更正错误 K′ 和自旋劈裂 | 真 K′ 是 \(−K=(-1/3,-1/3)\)，旧 “Kp” \((-1/3,+1/3)\) 在 Γ–M 内部；K band35/36 87.5/35.3 meV 是近同向自旋相邻带差，不是已证实的反向自旋 Ising splitting |
| Stage8k–l | 从真实 WAVECAR 独立读取 spinor，记录能量与 EIG 匹配；修正 35/36 机器账本错位 | 定性的非简并 transverse **pseudo-spin** 线索；7 个 Stage8k 新 JSON 字段实际记录 band36，非标签 band35；原始文件保留并加勘误 |
| Stage8m | **真实两次**原 WAVECAR 的 VASP2WANNIER90 后处理，PBS 492698/492699、Exit 0；96 bands、441 k；EIG/MMN/SPN | 冻结-KS 21×21 初步 Wilson 数值**未通过对称性**：AB C₃v 张量各向异性，glide total 远离零。产生的是 **不合格 exploratory 数值**，不是材料预测 |
| Stage8n | 首轮诊断，声称 VASP ISYM=0 数据不对称导致 | **此排他性根因主张被后续 Stage8o 方法审查否定**：用错 Cartesian mirror、只比 12/441 k 对、漏旋转 axial spin、误标 raw Wilson phase，Kubo 对照差异巨大 |
| Stage8o-real | 真 Seitz 操作、全 441 k、轴自旋 3×3 Gram，新增 Kubo 代码 | 真晶体对称下 Gram 中位 ~2×10⁻⁷，但 worst residual 可达 AB ~2×10⁻³、glide ~9×10⁻³；部分近简并非常严重；**真实 WSe₂ Wilson/Kubo/projector 三方数值未返回可复核结果**，根因仍开放 |
| **Stage8p（新）** | 独立非零两带解析模型，resolvent/projector/Wilson 三方 benchmark 与相位测试 | 新代码 \`stage8p_nonzero_benchmark.py\` 及测试 \`test_stage8p_nonzero_benchmark.py\` 已提交；**194 tests + 原论文复现 CI 两条 SUCCESS**。但仅校验二维两带方法符号/源因子，**不代表真实材料非零**，且合成模型并不要求时间反演对称 |

### Stage8o-real 关键数字（净化 JSON 可核对）
- AB，P3m1，6 Seitz 操作；GLIDE，Aem2，4 Seitz 操作；每一操作实现 441/441 网格双射。
- AB 非恒等晶体操作能谱误差 max 约 \(2.9–3.0\times10^{-4}\) eV、Gram residual 中位 \(1.8–2.4\times10^{-7}\)，最大约 \(2.4\times10^{-3}\)。
- GLIDE 面内镜面能谱误差 max \(3.1\times10^{-3}\) eV、Gram residual 中位 \(2.9–3.8\times10^{-7}\)，最差约 \(9.2\times10^{-3}\)。水平 glide 自身对 in-plane \(k\) 不变，能谱残差零，但 Gram 最差约 \(9.1\times10^{-3}\)。
- GLIDE 远离相邻带 \(>0.05\) eV 的态，面内对角 SPN 最大约 \(9.7\times10^{-3}\)，均值约 \(2.0\times10^{-4}\)，但未归因于唯一物理/数值机制。
- 原 Stage8n 的所谓 raw Wilson \(\sim10^{-3}\!-\!10^{-2}\) 实际已除以 \(\Delta k\Delta\lambda\)；更正的**原始 dimensionless 角度**约 \(3.9\times10^{-8}\)–\(5.3\times10^{-7}\) rad，B 差角 \(\sim7\times10^{-14}\)–\(1.2\times10^{-11}\) rad。机器 epsilon 裕量 \(\sim3\times10^2\)–\(5\times10^4\) 仍不能代表矩阵 SVD/det 的真实误差界。
- 因 AB C₃v 和 glide total null 没通过，**不要引用旧 \(10^{-9}\) 数值为 WSe₂ 真实响应**。不得将其与 toy model 的 \(10^{-6}\) 做未校准单位的三个数量级比较。

## 3. 最新 Stage8p：为什么重要、具体还差什么

**来源：**
- [stage8p_nonzero_benchmark.py](stage8p_nonzero_benchmark.py)
- [test_stage8p_nonzero_benchmark.py](../tests/test_stage8p_nonzero_benchmark.py)

Stage8o 的“合成验证”只用 symmetry-protected **零** 值：Wilson 和 Kubo 即使共同为零，也可能同时含有相反符号、错误 \(1/2\)、磁场因子或对称性伪消失。Stage8p 新增明确**非零**的可解析两带测试：

\[
H(k,\lambda,B)=\sin k\,\sigma_y+\frac{\lambda}{2}\sigma_x+(m-\beta B)\sigma_z,\quad
m=1.2,\ \beta=0.4.
\]
对应局部解析曲率
\[
F_{\lambda k}(k,B)=
\frac{(m-\beta B)\cos k}{4[\sin^2 k+(m-\beta B)^2]^{3/2}},
\]
并检验 \(\partial_B F\neq0\)。Wilson loop、\(\partial_\lambda P\) 的 resolvent 导数、中心差分 projector 三条路线在合成模型上比对非零值、符号、分母、source 1/2 及 band phase invariance。**请下一轮先复核 Stage8p 解析导数与实际基准测试是否覆盖足够 k、源步长、磁场步长和两种符号，不用再建立一个同功能 toy 方案。**

**严格限制：**此二带模型只是数学/实现基准，不是 TR-invariant 真实 WSe₂，不证明原有 96-band MMN 与 spin connection 正确；更不证明 material \(\Gamma\neq0\)。

## 4. 紧接着要做的第一件事情（按科研价值排序）

**目标：同一套现存真实 AB 96-band 数据，在精确同一个 k、同一个 \(\lambda,B\) 源定义下，拿到可审核的、非零的三路线对照。**

1. **核验 GitHub HEAD 和 CI**，阅读 Stage8o-real 已推送报告。运行 Stage8p 非零解析/回归，确认 resolvent、projector 和 Wilson 在**非零**模型上的数值一致，而不是仅凭两边同为零。
2. **修复 Stage8o 已有三方真实代码的输入约定：**\`stage8o_diag5.py\`、\`stage8o_kubo_cross.py\`、\`stage8m_response.py\` 必须统一 \(\lambda S_x\)、\(-g\mu_BBS_z\) 的数值单位。旧代码混用了 \(-\mu_BB\sigma_z\) 和 \(-2\mu_BB\sigma_z\)，以及 \(S_x\) 与 \(S_x/2\)；**这是新的真实三方比较尚未验收的重要代码问题**。
3. **纠正 Kubo/projector 的几何导数细节：**从真实邻接 MMN 得到不同 k 的 \(U(k')\)，构造 \(U^\dagger(k)M(k,k')U(k')\)。对 occupied projector 的 \(\partial_k P\) 应取真正的占据–空带块，并统一导数的反厄米连接与符号。核实所有 reciprocal boundary \(\mathbf G\)、实际方向度量和 Cartesian k 导数。不要用 \(U(k)\) 代替 \(U(k')\)。与真实 Wilson 的**局部 \(F\)** 比较应先固定 B，再比较 \(\partial_B F\)。
4. **先独立解析两带符号**，再只读使用原 Japan IMR EIG/MMN/SPN 在少量真实低对称 k 产生三方数据。原 Stage8o 已部署 \`stage8o_diag5.py\` 但完整 JSON 因登录网关截断取回失败。需要由 harness **直接在日本**保存结果为小型 sanitized JSON，只把摘要和 hash 提交 GitHub，避免把大量逐 k 记录送回聊天窗口。
5. **不要提前宣布物理 GO**：真实 AB tensor \(M=\gamma I\) 与 midpoint glide total=0 均应恢复，且三种估计器互相满足明确误差限，才讨论 WSe₂ 材料幅度。任何一项失败都要报告可归责的具体实现/输入差异。
6. **额外风险项（不急于加新作业）：**96-band window 的 band97 不可见，窗口顶部近简并，spin Gram 最大残差和真实晶体能谱非零残差尚未定量溯源；实际 sewing 还缺 nonsymmorphic phase/完整 MMN 的闭合验证。这些不能用 SPN Gram 中位很小一项自动消除。
7. **不做**旧稿中的生产级 Wannier、NEB、密 k 网格、真实 material layer-odd NMEE，除非先有可靠粗粒度 response 和单独授权的计算预算。

**首选下一阶段名称：Stage8p-real（非零基准迁移 + 三方真实数据核验）。** 本阶段原则上是**零新增 PBS/VASP 作业**，使用已有远端数据；若读取、内存或格式阻断，应记录恢复点，不得擅自启动新 SCF。

## 5. 关键路径

- [研究总索引](../RESEARCH_INDEX.md)
- [Stage8o-real 原始远端审查报告](STAGE8O_REAL_DIAGNOSIS_20261023.md)
- [Stage8o-real AB SPN/symmetry 机器证据](stage8o_AB_spn_symmetry_screen.sanitized.json)
- [Stage8o-real glide SPN/symmetry 机器证据](stage8o_GL_spn_symmetry_screen.sanitized.json)
- [Stage8o Kubo/core](stage8o_kubo_cross.py)
- [Stage8o 原真实三方对照脚本](stage8o_diag5.py)
- [Stage8m Wilson 基础引擎](stage8m_plaq.py)
- [Stage8m frozen-KS response](stage8m_response.py)
- [Stage8o 独立方法勘误](STAGE8O_INDEPENDENT_METHOD_AUDIT_20261023.md)
- [Stage8p 非零解析标定](stage8p_nonzero_benchmark.py)

### GitHub 工作边界

- 原有外部计算及 Stage8i–8o 源 JSON **只追加新版本与勘误，不覆盖或删除**；仅提交小型 sanitized 数值、程序和测试、来源 SHA。保留原 VASP/PBS 任务、POTCAR 与超算 SSH 信息在本地授权环境。
- 全部独立原研究仅写入 \`studies/bsf-nmee-connection/\`；论文侧 \`papers/jia-2026-nonlin-magnetoelectric-edelstein/research/README.md\` 仅做导航；Draft PR #3 不合并。
- 若进展需要重大计算资源，先给出**明确的新作业数、目的、失败回退、验收条件**，另行获得用户批准。用户不是人工复制命令、拷贝大型文件的中间人。

## 6. 交接后的输出契约

下一位执行者在读完本 handoff 后，不要从头复述原论文理论，也不要止步于“需要更多计算”。应直接交付：
- 一份 **Stage8p-real 中文科学报告**：真实三方对照的各个 k、B、λ 定义、数值、符号、相对误差及仍未通过的项；
- 新小型 JSON（包括来源 SHA、准确 k-index、矩阵维、spin/Zeeman normalization、gap guards、方法残差）和实际可复现的代码/测试；
- 真实 DFT 响应是否可被当作物理信号的 **GO/NO-GO**，并清晰地区分 **idea、方法、材料数值** 三个层次；
- GitHub commit、CI 状态，以及“是否需要额外计算”的具体理由。

**底线：Stage8p 合成非零基准已进入仓库且 CI 全绿；它不代表 Stage8m WSe₂ 材料数值已经恢复。下一个决定性证据是相同物理源下真实 EIG/MMN/SPN 的三方非零响应一致性及 AB/glide symmetry gates。**
