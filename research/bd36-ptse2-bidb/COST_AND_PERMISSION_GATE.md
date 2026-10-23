# BD36 | 计算资源与授权门（COST_AND_PERMISSION_GATE.md, rev R1.2）

> 状态：**READ_ONLY**（2026-10-23, R1.2）。未获明确批准前，**禁止** qsub / sbatch /
> VASP 生产作业 / 远端修改 / WSL·Slurm 新提交。授权范围仅限只读检查、输入生成、
> 静态测试与本地非生产性验证（均已执行，见 LOGS/）。

## 1. 只读环境盘点（R1.1 实测/档案核查，2026-10-23）

| 平台 | 调度 | VASP | Wannier90 | 赝势 | .spn 路径 | 状态 |
|---|---|---|---|---|---|---|
| **IMR**（日本，PBS Pro）| PBS（P_030: 112核/460GB；CP_001）| 队列存在；**VASP 版本无法从本沙箱核实**（两跳登录需 pty，本会话沙箱禁 pty）| **待用户侧确认**（VASP 6.4.3+ 含 LWRITE_SPN）| 待确认 | **计划走 VASP6 原生 `LWRITE_SPN=.TRUE.`** | **未核实** → 授权前需用户交互确认 |
| **para**（并行云，Slurm）| Slurm（root.bingxing, 空闲）| VASP 5.4.4 **vasp_ncl**（VaspWannier 构建）；**strings 证实含 LWANNIER90/LWRITE_SPN/LWRITE_UNK/WANNIER90_WIN 关键字** | `wannier90.x` **= 3.0.0**（zwy，需 MKL+MPI 环境，登录节点不可 MPI-init）| POTCAR_REC 含 **Pt、Se** ✓ | VASP-written .spn 关键字存在，但 **HgI2 基准从未端到端跑通（无任何输出文件）** → **端到端 .spn 输出未验证** | **关键字已验证；端到端待就地验证作业** |
| WSL（win）| 本地 20核/7GB | VASP 6.6.0（vasp_ncl ✓）| postw90.x/w90chk2chk.x ✓（版本待查）| Se ✓；**Pt 缺失**（批准后补齐）| VASP6 LWRITE_SPN 路径 | 二进制齐全；内存 7GB 受限 |
| xingyi（广州）| Slurm | — | — | — | — | 组配额 98.6% 满，不推荐 |

**IMR 工具链核查结论（如实报告）**：本沙箱（两跳、密钥带口令、禁 pty）**无法实际登录 IMR**，
因此 IMR 的 VASP 版本、Wannier90、POTCAR(Pt/Se) 均为**待确认**状态。**授权申请中必须附带用户对
IMR 工具链的确认**（版本 ≥ 6.4.3、vasp_ncl、wannier90.x、POTCAR Pt/Se 路径）。

**R1.2 关键官方核实（VASP wiki, 2024-10 版）**：`LWRITE_SPN` **仅 VASP > 6.4.2 提供**
（wiki 原文 "Available for VASP version > 6.4.2"，且仅支持默认 SAXIS=0 0 1）→
**para 的 VASP 5.4.4 官方不支持 LWRITE_SPN**（其 VaspWannier 二进制含该关键字属自定义构建，
输出正确性未验证）→ **para 第 0 步的裁决意义进一步上升**：通过则 para 可用；不通过则
.spn 主路径只能是 **IMR（≥6.4.3）或 WSL（6.6.0）**。.spn 格式已按 wannierberri 生产解析器 +
wannier90 `conv_write_spn` 源码确认：头 `nbnd nk`、每 k 点 3×nbnd(nbnd+1)/2 行 "re im"、
**Bloch 本征态基**、压缩下三角（每三角元素 3 分量连续写）、对角虚部≈0；
Wannier 基算符经规范变换 S^W=U† S U（`read_spn.py` R1.2 已按此重写并过已知矩阵单元测试）。

**para 核查结论（R1.2 更新）**：wannier90.x=3.0.0（`.spn` 依赖 VASP 写出；官方文档该特性仅 >6.4.2，见上）；
VaspWannier vasp_ncl 二进制含 LWRITE_SPN 等关键字（strings 实证），但 HgI2 L1 基准
目录**无任何输出**（eig/amn/mmn/UNK/.spn 全无）→ **端到端从未验证**。因此 para 路径的
生产前必须通过一个**就地接口验证作业**（见 §3 方案 A 的第 0 步）。

## 2. 作业核时上限（R1.2 分阶段修正，INPUTS/caps.json）

- 单作业：**28 核 × ≤2 h（≤56 核时）**；收敛测试作业：**≤28 核时/个**；
- **Phase 0（para 第 0 步就地接口验证，1 作业）：56 核时**；
- **Phase 1（Round-1 最小集）：448 核时** = 3 结构 × 2 作业 × 56 + 收敛 4 × 28；
- **合计 504 核时**（= 56 + 448）。**R1.2 起不再使用“≤500”单一表述**——
  504 > 500，任何文档/检查器出现 500 一律判 FAIL（反例测试覆盖）；
- 静态检查器 `static_check_inputs.py` 强制校验分阶段上限（0 FAIL 才可进入授权队列）。

## 3. 首个可执行 SOC+Wannier+spin-matrix 作业方案（待单独授权）

### 方案 A：para（全链关键字已验证；需就地验证作业）
第 0 步（**接口验证作业，1 作业 ≤56 核时**，2253 基准 9×9×1 全流程缩短版）：
1. 静态 SOC（INCAR_soc, CHGCAR）→ 提取实际 k 列表（`make_win_kpoints.py`）；
2. `wannier90.x -pp`（win 的 kpoints 块已生成）；
3. VASP 接口（INCAR_soc_w90_vasp544, LWANNIER90+LWRITE_SPN, ISYM=−1）→ 检验
   `.eig/.amn/.mmn/.spn/UNK` 生成、`read_spn.py` 维度/Hermite 校验、`hr.dat` 插值误差；
4. 通过 → 进入 Round-1 三结构正式流程；不通过 → 报告失败诊断，改走方案 B 或 C。
第 1-3 步：Round-1（2253、+u 滑移、2252）registry 弛豫 → SOC 静态 → SOC+W90 接口。

### 方案 B：IMR（任务首选；**先经用户确认工具链**）
同流程，INCAR_soc_w90_vasp6 + PBS_draft.pbs（P_030）；版本确认后直接进入正式流程。

### 方案 C：WSL（需补齐 Pt 赝势 + 内存确认）
VASP 6.6.0 vasp_ncl + LWRITE_SPN 原生路径；7GB 内存对 6 原子 SOC 小胞通常可行。

## 4. 需要用户批准的事项（单一合并请求）

1. **IMR 工具链确认**（版本 ≥ VASP 6.4.3 / vasp_ncl / wannier90.x / POTCAR Pt+Se 路径）——
   若确认则首选方案 B；
2. 否则批准**方案 A（para）含第 0 步就地接口验证作业**（Phase 0 = 56 核时，见
   `PARA_STEP0_SCHEME.md` 独立方案）；
3. **Round-1 作业范围**：2253 基准 + +(0.08,0.19) 滑移 + 2252 单层
   （最小集 6 作业 + 收敛表，Phase 1 = 448 核时；合计 56+448 = 504 核时，分阶段执行）；
   不动 2255/2254/反向滑移；
4. **费米能方案**：响应阶段电子/空穴参考掺杂取导电态（σ 非奇异才换算 χ^j）。

未批准前：只读检查、输入生成、静态测试、本地干跑全部完成（本分支）；**零生产提交**。

## 5. 私密数据边界

- POTCAR、远程凭据、平台内部路径**不进入提交内容**（R1.1 已清理所有提交脚本内的绝对路径，
  改由本地 `run_env.sh` 提供，已 gitignore）；执行时由已授权注册表装配并本地记录 SHA256。
- 原始 BiDB SQLite 数据库不提交；main 与既有分支未动；Draft PR #99 仅含结构/输入/脚本/日志/文档。
