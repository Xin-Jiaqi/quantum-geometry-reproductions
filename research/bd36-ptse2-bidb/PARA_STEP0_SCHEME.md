# BD36 | para 第 0 步：就地接口验证作业方案（PARA_STEP0_SCHEME.md, rev R1.3）

> 单独交付（R1.2 要求）。严格 READ_ONLY 期间**只提交方案，不提交作业**；批准后由
> 授权执行（1 作业，Phase 0 = **56 核时**：28 核 × ≤2 h，`sbatch` 一次）。

## 0. 背景与裁决目标（R1.3 补丁源码核实，非 strings 判定）

- **补丁来源已确认**：para 的 "VaspWannier" = [Chengcheng-Xiao/VASP2WAN90_v2_fix]
  (https://github.com/Chengcheng-Xiao/VASP2WAN90_v2_fix) 补丁族（`mlwf.F`/
  `wave_high.F`，目标 `vasp.5.4.4.pl2_w90`，2022-04-06 版）。.spn 功能由该补丁
  2020-06-25 加入。**关键源码事实**：
  - `WRITE_SPN` **仅串行**：`IF (WRITE_SPN .AND. WDES%COMM%NCPU.NE.1) WRITE_SPN=.FALSE.`
    （静默跳过！）；`CALC_SPN_EXP` 注释 "This routine only works in serial"。
  - **LSPN_FMTED**（补丁 INCAR 关键字，默认 `.FALSE.`）：`FALSE` → `.spn` 为
    **unformatted (binary)**（头记录 + `nbnd nk` 记录 + 每 k 点 complex128 记录）；
    `TRUE` → **formatted ASCII**（`2ES26.16`，每三角元素 3 分量连续）。
  - 存储方向：`spn_n_loop DO N=1,NB_TOT / spn_m_loop DO M=1,N /
    spin_mat(NI,s)=SPN(M,N,s)` → 按 (m 外循环, n=1..m) 顺序存储**上三角元素
    S[n,m]（n≤m）**、3 分量依次排列（read_spn.py R1.3 已按此修正并经独立
    σ_y 符号测试验证：M[0,1]=−i 判别，旧方向会给 +i）。
  - VASP wiki（2024-10）确认原生 `LWRITE_SPN` **仅 VASP > 6.4.2**（VASP 6 原生
    恒为 formatted，自动在 win 加 `spn_formatted = .true.`；无 LSPN_FMTED）。
- para 的二进制 strings 含 LWRITE_SPN 等关键字与上述补丁一致（补充证据），但
  HgI2 基准目录**从未跑出任何输出**（无 eig/amn/mmn/UNK/.spn）→ **端到端从未验证**。
- 第 0 步 = **一次 2253 基准的最小化全链作业**，裁决：para 能否产出可被
  `read_spn.py`（R1.2 重写版）解析的正确 .spn。
  - **通过** → 走方案 A 完成 Round-1（Phase 1 = 448 核时）；
  - **不通过** → 输出失败诊断（哪个文件缺失/格式错误），改走 IMR（≥6.4.3）或
    WSL（6.6.0）原生 LWRITE_SPN 路径，para 仅承担 SOC 静态部分。

## 1. 作业内容（单作业，2253_bilayer_baseline，9×9×1 全 BZ = 81 k）
执行脚本：`INPUTS/2253_bilayer_baseline/phase0_step0.sbatch`（完整可执行；批准后
`sbatch` 一次，28 核 × ≤2 h = ≤56 核时）。

| 阶段 | 内容 | 预期产物 | 校验 / 退出条件 |
|---|---|---|---|
| a | SOC 静态（INCAR_soc：ISYM=−1 → **与接口同一完整 81 k 全 BZ**，CHGCAR 供接口） | CHGCAR, OUTCAR | EXIT-1/2：静态失败、OUTCAR 无 k 列表 |
| b | `make_win_kpoints.py OUTCAR` → win kpoints 块 | kpoints_block.txt | EXIT-3/4：提取失败、k 计数≠81 |
| c | 组装 win（去占位 kpoints 块 + 追加生成块）→ `wannier90.x -pp` | wannier90.win 展开 + kpt 文件 | EXIT-5：-pp 失败（MPI/MKL 环境） |
| d | 接口（INCAR_soc_w90_vasp544：LWANNIER90+LWRITE_SPN+**LSPN_FMTED=.TRUE.**+ISYM=−1+KPAR=1；**单 MPI rank 串行**——补丁约束） | .eig .amn .mmn **.spn** | EXIT-6/7：接口失败、OUTCAR 无 "SPN calculated"（串行未生效则静默跳过） |
| e | 真实验收 | — | EXIT-8/9/10：必需产物缺失（**UNK 非必需**，LWRITE_UNK=F）、`read_spn.py` 拒绝（维数须 (81,64,64,3)、Hermite）、接口 k 列表≠win kpoints |

- **墙钟预算**：静态 28 核 ≤1 h；接口 1 核串行 ≤1 h → 单作业 ≤2 h = 56 核时。
- **k 一致性（R1.3 修复）**：预处理与接口 **均 ISYM=−1**，同为 81 k 全 BZ；win kpoints
  块由预处理 OUTCAR 生成，接口实际 k 列表与其 diff 一致（EXIT-10）。

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
