# Stage 8k：非简并态的自旋非守恒预审（2026-10-22）

**证据等级：仅基于已经上传 GitHub 的 Stage8i WAVECAR 衍生 JSON，未独立读取原始日本超算 WAVECAR，也未通过完整 PAW 算符验收。零新增 PBS 作业。**

本轮自动分析器：[stage8k_existing_spin_prescreen.py](stage8k_existing_spin_prescreen.py)。冻结的衍生结果：[STAGE8K_NONDEGENERATE_SPIN_CANDIDATES_20261022.json](STAGE8K_NONDEGENERATE_SPIN_CANDIDATES_20261022.json)。之前的 [Stage8j 勘误](STAGE8J_INDEPENDENT_VALLEY_SPIN_ERRATUM_20261022.md) 仍然有效。

## 为什么转向非简并态？

在严格无简并能级，若 [H(k), Sz]=0，则 H 本征态可以选成 Sz 本征态，因而真正的物理自旋 <Sx>=<Sy>=0。反过来，若非简并态的物理面内自旋**可靠地**非零，就否定了严格的 Sz 守恒。此判据比 Γ 点任意规范下 Kramers 简并自旋方向更有判别力。

但这仍不是证明了静态 NMEE：真实 Γ 需要价带-导带自旋与速度的非对角矩阵元、物理层算符、全 BZ 量子几何和带间虚跃迁。

## 已有五个位置每结构八条带的筛选

筛选范围是各结构五组之前记录的 k 点，band 33–40 共 80 个状态，不是全 BZ。门槛：相邻带距离至少 5 meV；旧伪自旋横向模长至少 0.025（spin-1/2 convention）。此相邻能级判断**只对当前有限能带集成立**，需原始文件验证。

| 源结构/实际 k | band | 旧 WAVECAR 衍生的面内伪自旋 | 最近相邻能带距离 | 证据 |
|---|---:|---:|---:|---|
| AB，(-1/3,+1/3)，Γ–M 线的 2/3 处，**非 K′** | 35 | **0.167801** | **7.84 meV** | 非简并候选，仍需排除近带混合与错误归一化 |
| AB，(0.476190,0)，靠近 M 但非精确 M | 35 | **0.188710** | **10.29 meV** | 较强面内伪自旋候选 |
| AB，(0.142857,0.142857)，Γ–K 线上 | 35 | **0.033579** | **71.25 meV** | 强度较低，但离近带简并更远 |
| 滑移中点，与 AB 同一近 M 点 | 35 | **0**（按源数值精度） | **20.64 meV** | 只供同 k 对比；未证明同一 orbital branch |
| 滑移中点，与 AB 同一 Γ–K 点 | 35 | **0**（按源数值精度） | **16.38 meV** | 同上 |
| 滑移中点，旧误标 Kp，即 Γ–M 内点 | 35 | 0.236570 | **0 meV**（band35/36 精确简并） | **简并态，排除出非简并检验** |

在现有 80 个样本内，筛选器找到 **AB 15 个、滑移中点 0 个**满足双阈值的伪自旋候选。**这不是说滑移中点在整个 BZ 内不含自旋混合**。同一带号不保证两堆叠下轨道本征态连续对应。

同时，原 spin-texture JSON 叫作 Lambda 的点在 (3/21,3/21)，而旧 EIGENVAL CBM 在 (4/21,4/21)。**不能用前者的自旋代替 CBM 自旋。**

## 下一步的关键“假阳性”门槛

1. **WAVECAR 读入正确性**：由于之前的 pymatgen wavecar band_energy 字段明显错误，两组自旋通道 plane-wave coefficient / G-vector 排列、record offsets、归一化、时反配对必须重新核验。简单的伪波函数 Euclidean norm 不为 1 不一定就表示 VASP 存储错误（PAW 使用广义重叠），但这使“物理自旋大小”不能仅由已有 pseudo 值确定。
2. **PAW 物理自旋**：目前缺真正 PAW 投影/增广算符，而 PROCAR 文件只有 21 字节；若没有可信办法从已有原始输出构建算符，应报告 BLOCKED_PAW_SPIN_OPERATOR，不能虚构。
3. **同 k 同带号并非同 orbital**：AB 与 glide 比较需可信 orbital matching 或 overlap tracking。额外滑移镜面对响应的真正影响必须通过算符、选态和规范做验证。
4. **非简并态与“真实静态 Γ”之间仍隔一关**：即使物理自旋非守恒，真实占据带到空带的 Sx 与 vj 交叉矩阵元仍可因选择定则或积分抵消而总响应为零。

## 已备好的自动交接

[Stage8k read-only DeepSeek assignment](STAGE8K_DEEPSEEK_READONLY_HANDOFF.md) 指定原始超算环境中的既有 WAVECAR/EIGENVAL 只读重分析，没有新 qsub，没有 VASP 改动，没有 Wannier 或 NMEE 作业。若检查失败，提交可信的障碍和现有衍生证据即可，不得要求用户手动复制文件或运行大量命令。

**当前状态：NONDEGENERATE_SPIN_PROXY_IDENTIFIED / PHYSICAL_SPIN_AND_RAW_WAVECAR_NOT_VALIDATED。** 这使 Stage9 的严格 Sz 守恒 no-go 有了更明确的否定性检验目标，但没有证实真实 WSe2 静态 NMEE 非零。
