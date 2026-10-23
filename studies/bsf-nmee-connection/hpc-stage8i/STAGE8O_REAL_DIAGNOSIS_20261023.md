> **Stage8p 源码级勘误（2026-10-23）：** 本 Stage8o-real 原始报告与数值保留原样，但所称“修正 Kubo/projector 三方对照”的两条路线**没有使用相同的自旋源和 Zeeman 约定**：projector 分支对 Pauli SPN_z 多乘一次 g=2，对 SPN_x 没有除以 2，并将非厄米 basis connection `dM` 直接当成 `dP/dk`。因此此前即使能取回真实三方 JSON，也不能将其视为同一个 Γ 的数值验收。独立非零合成基准现已建立；真实 WSe2 一致性仍待只读复核。请先读 [Stage8p 科学审查](STAGE8P_NONZERO_METHOD_AND_REAL_DATA_REVIEW_20261023.md)。**

# Stage8o-real | 真实 441 点 Seitz/SPN 对称性与 Kubo 数值核验（2026-10-23）

> 自主执行 `STAGE8O_DEEPSEEK_ZERO_JOB_HANDOFF.md`：零新增 PBS/VASP/Wannier 作业，
> 只读日本 IMR 已归档的 Stage8m 目录（492698/492699 后处理产物），以真实
> spglib Seitz 操作 + 全 441 k + 轴自旋 Q 旋转 + 完整 3×3 spin-Gram 协方差
> 检验 Stage8n 的“VASP 数据层对称性破坏”根因主张是否成立，并修正 Wilson
> 相位预算与 Kubo 实现。

## 0. 状态判定

**`MATERIAL_RESPONSE_SYMMETRY_GATE_FAILED__STAGE8N_VASP_ROOT_CAUSE_NOT_ESTABLISHED__IDEA_OPEN`**（确认审查者状态）：
Stage8n 对继续生产级数值预测的 NO-GO 正确，但其根因主张
`ROOT_CAUSE_LOCATED_IN_DATA_SYMMETRY__ENGINE_AND_GRID_CLEAN` **不被修正后的真实数据检验支持**。
本报告给出证据与 GO/NO-GO 数值。

## 1. 输入与校验

- 目录：`soc_AB_s0.00_pp8m`（AB，P3m1）与 `soc_GLIDE_s0.50_pp8m`（GLIDE，Aem2）；
  POSCAR 均为 6 原子 W₂Se₄，spglib 2.7.0 在线判定 6 与 4 个 Seitz 操作。
- 441 k（21×21），96 带；EIG/WIN/SPN 全量读取，输入 SHA256 记入净化 JSON。
- 对称配对：对每个 spglib 分数旋转 R，用 `k′ = R⁻ᵀ k mod Z²` 映射 **441/441**
  全部网格点（已验证双射）；自旋为轴矢量，`Q = det(R_cart) R_cart` 旋转后比较
  完整 `C_ij = Re Tr[S_i† S_j]` 协方差（`C(k′) = Q C(k) Qᵀ`），并报逐点最差残差。

## 2. 结果 A：能谱与 spin-Gram 协方差（真实 Seitz 操作）

| 结构 | 操作 | 441 配对 | 能谱残差 max (eV) | Gram 残差 max | Gram 残差中位 | 最差残差位置 |
|---|---|---|---|---|---|---|
| AB | 恒等 | 441 | 0 | 1.8e-17 | 1.7e-17 | — |
| AB | 5 个非恒等（C₃/mirror） | 441 | 2.9–3.0e-4 | 4.4e-4–2.4e-3 | 1.8–2.4e-7 | k=0（Γ，带间距 9.5e-11 eV）、k=27/79/83 等（带间距 3.6e-5–5.5e-4 eV） |
| GLIDE | 恒等 | 441 | 0 | 1.8e-17 | 1.7e-17 | — |
| GLIDE | 面内镜像 ×2（含半平移） | 441 | 3.1e-3 | 8.5e-3–9.2e-3 | 2.9–3.8e-7 | k=350/357（带间距 1.07e-3）、k=35/36（2.6e-9） |
| GLIDE | 水平 glide（m_z，t=(.5,.5)） | 441（k 恒等） | 0 | 9.1e-3 | 2.5e-7 | k=35/36、k=357（近简并） |

**解读**：
- **两个结构在真实晶体操作下的 Gram 协方差中位均 ~2e-7**（441 点几乎全通过）；
  **最大残差（2e-3–9e-3）全部集中在近简并/准简并 k**（Γ 点带间距 1e-10、以及
  3.6e-5–1e-3 eV 的带对），与 96 带窗口顶部 95/96 近简并（Stage8m 已测 7e-5 eV / ~0）
  及带间相位混合一致。
- **Stage8n 的 0.32/1.03 “mirror SPN 残差”被证伪**：其 Cartesian `[[0,-1],[-1,0]]`
  对给定六角晶格不是整数 unimodular 变换（`stage8o_symmetry_math_audit.py` 量化
  非整残差 0.423），且未旋转自旋、只 12 对。修正后同量检验降至中位 2e-7。

## 3. 结果 B：水平 glide 的面内自旋选择规则（分层）

水平 glide `{m_z | t=(.5,.5)}` 对 in-plane k 恒等，非简并本征态必须 S_x=S_y=0。
全 441×96 带态按最小相邻带间距分层（band 96 无 97 邻居，排除；间隔<0.001 eV 跳过）：

| 分层 | 态数 | max |S_xy| | mean |S_xy| | >1e-4 占比 |
|---|---|---|---|---|---|
| 带间距 0.001–0.01 eV | 9080 | **0.735** | 0.011 | 4185/9080 (46%) |
| 带间距 0.01–0.05 eV | 17590 | **0.346** | 1.4e-3 | 7111/17590 (40%) |
| **带间距 >0.05 eV（well-separated）** | 10730 | **9.7e-3** | 2.0e-4 | 2319/10730 (22%) |

- 面内自旋残差随带隔离度单调下降（0.735→0.35→9.7e-3），与“近简并混合伪影”模式一致。
- **well-separated 态最大 |S_xy| ≈ 1e-2**：非零但小；不能据此证明“数据完全缺失 glide
  对称性”，也不能证明“数据含 glide”（SPN 带对角是伪势派生代理，非全空间层算符；
  96 带窗口无带 97 不封闭）。**原因归因保持开放**：候选为窗口截断 / SPN 代理 / 近简并，
  **不是已证实的 ISYM=0 失效**。

## 4. 结果 C：Wilson 相位精度预算（修正）

`stage8o_symmetry_math_audit.py` 确认 Stage8n 报告把**已除以 ΔkΔλ 的混合曲率**
（`plaquette_phase()` 返回值）标成“原始相位”。修正后（同一 AB/GLIDE 数据）：

- 原始环路角度 φ = |arg det L₁L₂L₃L₄| ≈ **3.9e-8 … 5.3e-7 rad**（非 Stage8n 报告的 4e-3）；
- 两场原始角度差 Δφ ≈ **7e-14 … 1.2e-11 rad**，|Δφ|/机器ε ≈ **3.2e2 … 5.2e4**。
- **结论**：多数点仍可分辨（高出机器精度 2–5 个量级），但个别点（GL ax=1 k=(0.31,0.22)，
  |Δφ|/ε≈3.2e2）接近数值噪声；**未做 36×36 极分解/行列式与网格求和的真实误差预算**，
  不能宣布“引擎精度已过关”（审查者观点成立）。

## 5. 结果 D：Kubo/projector 数值核验

- **实现修正**（`stage8o_kubo_cross.py`）：resolvent 公式改回单分母
  `2 Im[<n|S_i/2|m><m|∂_k u_n>]/(E_n−E_m)`，自旋矩阵元取 <n|…|m>（a* b 配对）；
  速度改用 **gauge-covariant 传输链接** `dM = (U(k)†M(k,k+)U(k+) − U(k)†M(k,k−)U(k−))/(2dk)`，
  自动包含 Stage8n 漏掉的 `U†∂_k U` 项；新增**投影算符有限差分路线**
  `i Tr P[dP/dλ, dP/dk]`（P 在 Zeeman 本征基，∂λP 由 λ±δλ 重对角化数值微分）。
- **合成验证**：SlidingBilayer（2H-AB 式）静态响应受 TR+反演对称性保护**恒为零**
  （3 个自旋源轴 × 全网格实测 1e-16 级）——物理正确，但该模型无法提供“非零远离简并”
  基准；破反演变体（层势 0.05 eV）因轨道层序歧义未达非零。**非零合成 benchmark 未达成**，
  如实记录。
- **真实数据三方对照**（AB 96 带，gap 1.15 eV 低对称 k，resolvent vs projector vs plaquette）：
  脚本 `stage8o_diag5.py` 已部署并 RC=0，但输出 JSON 因登录网关连续截断未能取回，
  **真实三方对照未完成**（限制记录）。独立 Kubo 对响应估计器的交叉验证**仍未通过**。

## 6. 对 Stage8n 逐项判定（修正后）

| Stage8n 主张 | 独立核验结论 |
|---|---|
| Cartesian `[[0,-1],[-1,0]]` 是 AB 对称镜像 | **否**（非整数变换，残差 0.423） |
| “SPN 镜像残差 0.32/1.03 = 数据破坏” | **否**（无效映射+未旋转自旋+12 对伪影；修正后中位 2e-7） |
| “原始相位 4e-3、|dP|/ε 1e9” | **否**（原始角度 1e-7 rad、|Δφ|/ε 3e2–5e4） |
| “VASP ISYM=0 数据层根因已定位” | **否**（根因未确定；候选：96 带窗口不封闭/近简并/SPN 代理） |
| C₃v 张量破坏、glide null 非零 | **仍成立**（Stage8m 数值未消失，但来源不明） |

## 7. GO / NO-GO（科学判定）

- **NO-GO（生产级定量预测）**：维持。AB C₃v 笛卡尔张量残差仍与分量同量级、glide 中点
  总响应仍显著非零、独立 Kubo 交叉验证未通过——三个 STOP 值均未满足。
- **科学 idea：OPEN（需下一轮单独授权计算）**。Stage8n 的“VASP 数据对称性破坏”根因
  不成立（本报告证伪）；数据在真实晶体操作下**中位符合**对称性，残差集中在近简并。
  下一轮最小授权范围：带 97 以上的窗口封闭 sewing 诊断（或轨道级 SPN/MMN 重构）、
  规范固定的独立 Kubo（平行传输）、排除近简并点后的 C₃v/glide 张量重算与单位修正
  （Å⁻¹ eV⁻¹ T⁻¹）。

## 8. 硬边界与交付

- 零新增 PBS/VASP/SOC/NEB/Wannier；零原始 MMN/SPN/WAVECAR/POTCAR/CHGCAR 上传；
  无私人/许可/凭据文件公开。
- 交付：本报告 + `stage8o_actual_spn_symmetry_screen.py`（含分层统计）+ 
  `stage8o_kubo_cross.py` + `stage8o_diag5.py` + 净化 JSON（AB/GL screen）+ 
  本地测试（新增 Kubo 零对照 + 既有 Stage8o math/SPN 测试）；更新
  `RESEARCH_INDEX.md` 与论文侧导航，PR #3 保持 Draft。
