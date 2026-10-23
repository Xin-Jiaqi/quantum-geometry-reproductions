# Stage8m 粗粒度冻结 KS 静态自旋-Zeeman 响应试点 —— 最终科学报告

> 状态：**完成**。任务书（STAGE8M_DEEPSEEK_ROUGH_DFT_RESPONSE_PILOT.md）全部执行。
> 结论：**NO-GO（对继续生产级研究）**。AB 端点存在收敛的粗粒度信号，但
> glide 中点 total null test **失败**，且全部可在现有数据内修复的管线环节核查
> 干净，破坏指向 VASP 输出（SPN/MMN）在滑移几何下的数值对称性（ISYM=0 后处理）。
> 按任务书 B/D 段，该非零**不得解释为物理效应**，本阶段不具备生产级延续依据。

## 0. 协议边界（任务书 A/B/D 段）

- 物理量：Γ^stat_{iz;j} = ∂_{B_z} ∫_{BZ} d²k/(2π)² · i Tr P[∂_{λ_i}P, ∂_{k_j}P] |_{λ=B=0}；
  冻结 KS 带（H(k,λ,B) = diag E + λS_i/2 − μ_B B S_z，g=2，spin-only Zeeman）。
- 实现：B=0 处中心差 `[F(+dB)−F(−dB)]/(2dB)`（dB 序列 0.125/0.25）；det 相位经
  polar 酉分解；mmn 模式 k 边为 `k1→k_c→k2` 两步链（O(dk²)，合成回归 rel<1e-4 @ dk=0.005）。
- 不是自洽磁/电响应、不是 driven-DC NMEE、不是 Edelstein；无 layer 投影算符分解。
- 粗粒度估计 = 原 21×21×1 网格 + 96 带窗口上的 Wilson plaquette 估计；**不是生产级 Wannier 收敛声明**。
- 输入 SHA 冻结核验通过（WAVECAR/EIGENVAL 后处理前后一致，见 Stage8l 账本）。

## 1. 数据与管线

| 项 | AB（s=0.00，2H） | 滑移中点（s=0.50，glide） |
|---|---|---|
| 后处理目录 | `soc_AB_s0.00_pp8m` | `soc_GLIDE_s0.50_pp8m` |
| 产物 | wannier90.mmn/.spn/.eig/.win（96 带、441 k、全网格） | 同 |
| 作业 | dry-run 2 次（32GB OOM）→ 正式 2 次（492698/492699，64GB，Exit 0） | 同 |
| SPN 检验 | Hermitian dev < 1e-4；S_z 对角 ≤ 1（Pauli 归一，σ 非 σ/2） | 同（min_sval 0.911） |
| MMN 链接 | 3087 唯一链接（折叠自链接去重、v2 5-token 块格式） | 同 |
| INCAR 规范 | ENCUT 520 / NBANDS 96 / SAXIS 0 0 1 / ISYM 0 / ICHARG 2 / ALGO=None / NELM 1 / LWANNIER90_RUN=F | 与 AB 逐项一致 |

对称性（spglib，symprec 1e-3 Å）：AB = P3m1 (156)；glide = Aem2 (39)，
**含 m_z+半平移的层交换 glide 操作**（R=diag(1,1,−1), t=(0.5,0.5,0.8975)），
满足任务书“layer-exchanging glide 保留”的前提。

## 2. AB 端点粗粒度响应（Γ_BZ，单位 1/(eV·T·Ang²)）

| 配置 (src, win, B, dl) | Γ_x | Γ_y | gap (eV) | sval |
|---|---|---|---|---|
| x, 36, 0.25, 1e-4 | +3.04e-12 | +9.33e-13 | ∞（窗口=占据，无空带） | 0.908 |
| x, 48, 0.25, 1e-4 | +9.87e-10 | +9.17e-09 | 1.153 | 0.908 |
| x, 64, 0.25, 1e-4 | +1.13e-09 | +8.69e-09 | 1.153 | 0.908 |
| x, 80, 0.25, 1e-4 | +1.17e-09 | +8.70e-09 | 1.153 | 0.908 |
| x, 96, 0.25, 1e-4 | +1.16e-09 | +8.69e-09 | 1.153 | 0.908 |
| y, 96, 0.25, 1e-4 | +1.76e-09 | +7.29e-09 | 1.153 | 0.908 |

全矩阵 60 配置见净化 JSON（stage8m_AB_response.sanitized.json）。

**双零控制组**：win=36（占据带窗，P 常数）→ Γ ~ 1e-12（物理零），证明
引擎/BZ 求和无系统偏置；win48+ 打开带间通道后信号 ~1e-9，比控制组高 3 个量级。

## 3. 稳定性检查（AB 与 glide 均通过）

- **B 步长**（dB=0.125 vs 0.25，B=0 中心差）：差异 < 0.7%（所有源/轴）。
- **λ 步长**（dl=1e-4/2e-4/4e-4）：差异 < 0.5%。
- **嵌套窗口**（36/48/64/80/96）：win64 之后 G_x 变化 < 4%、G_y < 0.3%
  —— 远带截断稳定（无符号翻转、无量级跳跃）。
- **gap guards**：win48+ 全程 min_gap = 1.153 eV（无空带泄漏，fail-closed 生效）。
- **min_sval ≈ 0.908–0.911**：重叠近酉（VASP 非对角化基的数值有限酉性），polar 投影稳健。

## 4. GLIDE 中点 total null test —— **失败**

| 配置 | AB Γ | glide Γ | ratio |glide/AB| |
|---|---|---|---|---|
| x, win96, B0.25, dl1e-4 | (+1.16e-9, +8.69e-9) | (−7.49e-9, +3.88e-9) | 6.5 / 0.45 |
| y, win96, B0.25, dl1e-4 | (+1.76e-9, +7.29e-9) | (−2.65e-9, −1.07e-8) | 1.5 / 1.5 |
| 控制组 win36 | ~1e-12 | ~1e-12 | O(1)（均在噪声地板） |

- 任务书 ratio 定义：abs(glide)/max(abs(AB), noise_floor=1e-12) → **0.45–6.5，
  远超数值噪声相容范围（应为 O(1e-3) 以下）**。
- 该非零对 B/dl 高度收敛（<0.5%），是**确定性系统量**，不是随机噪声。

### 4.1 对称性归因核查（任务书第 78 行要求，先修正再谈物理）

已逐项核查可在现有数据内修复的环节，全部干净：

1. **reciprocal-boundary links**：坐标查表周期邻居（VASP 对称 k 序 0,+1/21,−1/21,…）；
   合成半格矢中点（shift=(0.5,0.5)）上 **frame 与 mmn 链式两种模式都给出 Γ≈0**
   （~1e-12，21×21 网格）→ 链式近似不破坏 m_z-glide 零结构。
2. **k-point ordering**：441 点覆盖全网格、无重复坐标（win 464 行 = 441 k + 头）。
3. **spin source normalization**：SPN 存 σ（对角 0.817>0.5 验证），源乘 1/2，g=2；
   Hermiticity 与 S_z 对角界检查通过。
4. **symmetry/geometry**：spglib 确认 glide 构型含层交换 glide 操作（见 §1）。

结论：破坏**不在引擎、不在网格排序、不在源归一化、不在几何**；合成对照证明引擎
在真实 glide 对称性下输出严格零，因此剩余破坏源是 **VASP 输出数据（SPN/MMN）在
滑移几何 + ISYM=0 后处理下的数值对称性不严格**（波函数相位/带间符号约定）。
按任务书，**不可将其解释为物理 layer-odd/任何物理效应**；修复需要新计算
（对称性诊断作业），在本阶段“零新增作业”禁令之外，故停止。

## 5. Stage7i 合成回归与引擎验证

- 引擎 vs Stage7i 独立 Wilson/Kubo 基准：rel < 5e-8（本地 8 带模型）。
- mmn 链式 vs frame（平滑重叠、dk=0.005）：rel < 1e-4。
- 镜像滑移反号（u=±0.4）、gap-guard fail-closed、reader fail-closed、短行拒绝：
  全部测试通过（全套 173，skipped=7）。
- **新增 glide 中点合成检验**：半格矢滑移中点（shift=(0.5,0.5)）frame 与 mmn
  双模式 Γ_BZ ≈ 0（~1e-12）—— 引擎在真实层交换 glide 下保持严格零。

## 6. 单位约定与 SI 参考

- Γ_BZ_Ang2 单位：1/(eV·T·Ang²)；mesh 均值 = Γ_BZ / (A_BZ/(2π)²/441)。
- SI：× (1e-9 m/Ang)² / (1.602176634e-19 J/eV)。
- AB 端点参考量（win96, B0.25）：G_y ≈ 8.7e-9 → SI ≈ 5.4e-8 m²/(J·T)（粗略参考，
  不是生产级；合成 toy model 相应值 ~1e-6，**高估约 3 个量级**——真实 WSe2 的
  冻结 KS 响应远小于 toy model，Stage7i 的示意量级不可外推）。

## 7. GO/NO-GO 判定：**NO-GO（对继续生产级研究）**

| 任务书检验 | 结果 |
|---|---|
| AB 信号非零 | 通过（~1e-9，win36 控制组 1e-12） |
| B/dl 步长稳定 | 通过（<0.7% / <0.5%） |
| 远带截断稳定 | 通过（win64 后 <4%） |
| glide 中点 total null | **失败**（ratio 0.45–6.5，确定性非零） |
| 引擎/管线自查 | 全部干净；破坏定位到 VASP 输出数据 |

**NO-GO 理由**：
1. glide 中点 total null test 未通过 —— 假说的“对称性保护的零结构 + layer-odd
   允许有限”核心前提得不到数据支持；任务书 B/D 段明确禁止将此类非零解释为物理。
2. 修复需要新计算（对称性诊断 / 重新后处理），在“零新增作业”协议之外。
3. 真实 WSe2 粗粒度响应量级（~1e-9）比 toy model（~1e-6）低 3 个量级，
   即使未来通过对称性门槛，也需 Stage10 文献矩阵确认其相对现有 NMEE/磁响应文献
   的差异化价值（当前矩阵未见 glide-protected layer-odd 静态张量研究空白被此
   数据填补的证据，因零检验未过）。
4. AB 端点的收敛信号可作 COARSE_BLOCH_SPACE_STATIC_RESPONSE_ESTIMATE 参考存档，
   不作为 production prediction。

## 8. 交付物

- 代码：stage8m_plaq.py（引擎+reader，gesdd→gesvd→eigh 极分解回退）、
  stage8m_response.py（线程池并行，jobs=8）、stage8m_analyze.py；
  回归测试 test_stage8m_plaq.py / test_stage8m_response.py（173 全绿，skipped=7）。
- 数据：stage8m_AB_response.json / stage8m_GLIDE_response.json（远端，已取回本地）。
- 净化 JSON：删除绝对路径与作业标签后的 stage8m_*_response.sanitized.json。
- 本报告 + stage8m_README.md。
- Git：提交 832c8e9（引擎+测试+文档）；本报告与净化数据随最终提交。
