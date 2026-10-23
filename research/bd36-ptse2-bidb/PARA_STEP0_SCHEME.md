# BD36 | para 第 0 步：就地接口验证作业方案（PARA_STEP0_SCHEME.md, R1.2）

> 单独交付（R1.2 要求）。严格 READ_ONLY 期间**只提交方案，不提交作业**；批准后由
> 授权执行（1 作业，Phase 0 = **56 核时**：28 核 × ≤2 h，`sbatch` 一次）。

## 0. 背景与裁决目标（R1.2 官方核实后）

- VASP wiki（2024-10 版）确认 `LWRITE_SPN` **仅 VASP > 6.4.2 提供**，且只支持默认
  `SAXIS = 0 0 1`（我们的 INCAR 恰为默认值 ✓）。
- para 的 VASP 5.4.4（VaspWannier 自定义构建）**官方不支持** LWRITE_SPN；其二进制
  strings 含 LWRITE_SPN 关键字（R1.1 实证），但 HgI2 基准目录**从未跑出任何输出**
  （无 eig/amn/mmn/UNK/.spn）→ **端到端从未验证**。
- 第 0 步 = **一次 2253 基准的最小化全链作业**，裁决：para 能否产出可被
  `read_spn.py`（R1.2 重写版）解析的正确 .spn。
  - **通过** → 走方案 A 完成 Round-1（Phase 1 = 448 核时）；
  - **不通过** → 输出失败诊断（哪个文件缺失/格式错误），改走 IMR（≥6.4.3）或
    WSL（6.6.0）原生 LWRITE_SPN 路径，para 仅承担 SOC 静态部分。

## 1. 作业内容（单作业，2253_bilayer_baseline，9×9×1 全 BZ = 81 k）

| 阶段 | 输入 | 说明 | 预期产物 | 校验（脚本） |
|---|---|---|---|---|
| a | INCAR_soc + KPOINTS(9 9 1) | SOC 静态（LSORBIT, SAXIS 默认, ISYM=1 收敛用 / −1 接口用） | CHGCAR, OUTCAR（含实际 k 列表） | `make_win_kpoints.py OUTCAR` → k 计数/顺序/权重 |
| b | wannier90.win（kpoints 块由 a 生成） | `wannier90.x -pp`（3.0.0，需 MKL+MPI 环境） | *.win 展开（含 kpoints 块） | diff 与 a 的 k 列表一致 |
| c | INCAR_soc_w90_vasp544（LWANNIER90=T, LWRITE_SPN=T, ISYM=−1） | VASP 接口 | .eig .amn .mmn **.spn** UNK | `read_spn.py <seed>.spn 81`：头 `nbnd nk`、Bloch 基、Hermite、维数 (81,64,64,3) |
| d | （可选）wannier90.x 主运行 | 插值 hr.dat | hr 插值误差报告 | 插值误差 < 1e-5 量级 |

- **墙钟预算**：28 核并行，SOC 静态 ≤1 h + 接口 ≤1 h → 单作业 ≤2 h = 56 核时。
- **失败退出条件**（任一命中即停，不继续烧核时）：
  1. OUTCAR 无 k 列表 / k 数与 81 不符；
  2. 接口作业未生成 .spn（或尺寸 ≠ (81,64,64,3)）；
  3. `read_spn.py` 报格式错误（头解析失败 / 非实对角 / 行数不匹配）；
  4. wannier90.x -pp 或主运行 MPI-init 失败（登录节点已见 PMPI_Init 中止）。

## 2. 环境装配（批准后执行时）

- `module load`/`export`：MKL 2022.1.0 intel64、oneapi/mpi 2021.6.0 lib + lib/release；
  VASP 5.4.4 vasp_ncl 显式 export（R1.1 已记录 para 登录节点实测路径）；
- `run_env.sh`（本地装配，gitignore，不入库）提供绝对路径与 POTCAR（Pt/Se）SHA256 记录；
- 作业脚本：`INPUTS/2253_bilayer_baseline/sbatch_draft.sh` 为模板，第 0 步用缩短版
  （不改 INCAR 物理参数，仅降 NSW=0 与接口步）。

## 3. 与预算的关系（R1.2 分阶段）

- Phase 0 = 56 核时（本方案）；Phase 1 = 448 核时（Round-1 最小集）；合计 **504 核时**。
- caps.json 已改为分阶段字段（phase0_validation_core_h=56, phase1_round1_core_h=448,
  total_core_h=504）；`static_check_inputs.py` 对"单一 ≤500"表述判 FAIL（反例已测）。
- Round-1 中 2253 的接口作业（INCAR_soc_w90_vasp544）在 Phase 0 已做过同结构全链验证，
  不重复计入额外核时（仍按 6 作业最小集计 448）。

## 4. 通过/不通过后动作

- 通过：按授权范围提交 Round-1（2253 基准 + +(0.08,0.19) 滑移 + 2252 单层）；
- 不通过：本方案记录失败诊断；授权申请转为 IMR（用户确认 ≥6.4.3/vasp_ncl/POTCAR）或
  WSL（补齐 Pt 赝势）的 .spn 主路径；para 保留 SOC 静态角色。
