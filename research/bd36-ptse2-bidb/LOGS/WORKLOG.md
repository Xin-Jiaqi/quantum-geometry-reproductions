# BD36 R1 worklog (2026-10-23)

- Read DEEPSEEK_TASK.md + PR #98 inputs (6 POSCARs + manifest + README).
- Built pure-numpy site-permutation symmetry audit; results match BiDB manifest
  (layer groups 72 p-3m1 / 69 p3m1; slides retain inversion, break C3+mirrors;
  permutation errors ~1e-8 A; slide shift measured = declared +-(0.08,0.19),
  layer gap preserved). See STRUCTURE_SYMMETRY_AUDIT.md + LOGS/audit_run.log.
- Environment gate (read-only): para = verified SOC+W90 chain (vasp_ncl,
  wannier90.x, Pt/Se POTCAR_REC, HgI2 L1 template); IMR versions unverifiable
  from sandbox (pty); WSL has Se POTCAR, Pt missing, VASP6+postw90 present.
- Materialized INPUTS for 6 structures (3 job types each + KPOINTS + POSCAR_SD
  + wannier90.win drafts + PBS/Slurm drafts), all NOT_RUN, hashed in manifest.
- Wrote PLAN.md / STRUCTURE_SYMMETRY_AUDIT.md / COST_AND_PERMISSION_GATE.md /
  INPUTS/CONVERGENCE_TABLE.md.
- No production compute submitted (permission gate).

## BD36 R1.1 corrections (2026-10-23, per user review of Draft PR #99)

1. POSCAR_SD selective dynamics: bilayers x,y frozen / z free (T F F); monolayer
   full relax (T T T). Post-relax verification added (verify_relaxed.py): inversion/
   C3/mirror re-audit + registry drift + slide residual; dry-run PASS (baseline D3d
   retained; slides Ci retained; slide residual ~1e-16 vs declared +-(0.08,0.19)).
2. Wannier projections unified: Pt:d + Se:p only; num_wann = 2*(5*nPt+3*nSe) =
   44 (bilayer) / 22 (monolayer); num_bands=64 >= num_wann.
3. VASP/W90 interface fixed: LWANNIER90=.TRUE. + wannier90.win as FILE; win
   atoms_frac generated per structure from that POSCAR (slide-correct, verified);
   win kpoints NOT hard-coded - make_win_kpoints.py generates post-SCF from the
   actual VASP k-list (ISYM=-1, 9x9x1 -> 81); no embedded WANNIER90_WIN, no
   spin=true. Platform variants: vasp6 (IMR) / vasp544 (para).
4. Spin matrix path: LWRITE_SPN=.TRUE. (VASP-written .spn). para read-only probe:
   wannier90.x == 3.0.0; VaspWannier vasp_ncl binary CONTAINS LWANNIER90/
   LWRITE_SPN/LWRITE_UNK keywords (strings), BUT the HgI2 L1 benchmark has NO
   outputs -> end-to-end .spn NOT validated -> in-place validation job required
   before production. IMR: version UNVERIFIABLE from this sandbox (no pty);
   LWRITE_SPN by VASP>=6.4.3 design, pending user confirmation. read_spn.py added.
5. Total spin (S_i = <w_n|sigma_i|w_m>, .spn) vs layer-projected spin
   (S_layer = P_layer S_z P_layer by Wannier centers) distinguished in PLAN.md
   §5; response via PT-degenerate-subspace projected traces (gauge-safe).
6. Caps fixed: per job 28 cores x <=2 h (<=56 core-h); round-1 total <=500 core-h
   (minimal set 448); committed scripts cleaned of internal absolute paths
   (env placeholders + local run_env.sh, gitignored).
7. static_check_inputs.py: 17 checks/structure + budget; 6/6 structures PASS
   (0 FAIL), exit 0. make_win_kpoints.py + read_spn.py + verify_relaxed.py added.

No production compute submitted. PR #99 updated with per-item PASS/FAIL.

## R1.2（2026-10-23，READ_ONLY）修正与验证

七项修正（commit 摘要）：
1. **选择性动力学语义纠正**：VASP 语义 T=自由、F=固定。双层 registry-preserving
   由 `T F F`（写反：实为固定 z、放开 x,y）改为 **`F F T`（x,y 固定、z 自由）**；
   单层保持 `T T T`。同步 build_inputs.py / 全部 INPUTS README / static_check_inputs.py
   （A2 语义文档检查，独立于生成器）。
2. **read_spn.py 重写**（真实格式）：头 `nbnd nk`；每 k 点 3×nbnd(nbnd+1)/2 行
   "re im"；**Bloch 本征态基**（非 Wannier 基）；压缩下三角（每三角元素 3 分量连续，
   wannier90 conv_write_spn + wannierberri 交叉确认）；对角虚部≈0 校验；
   Wannier 规范变换 S^W(k)=U† S(k) U。已知矩阵单元测试（4 带 σ_x/σ_y/σ_z、置换 U、
   独立 numpy 对照、Hermite）全过。
3. **make_win_kpoints.py**：修复跳首 k 点 bug（k 行从 index 3 起，官方 3 行头）；
   C1 计数/C2 格式/C3 权重(和=nk)/C4 Γ 首点/C5 分数坐标（拒 Cartesian，|k|<1）；
   VASP 官方显式示例（Cartesian）→ 拒绝；81 点全 BZ 夹具 → 通过；3 点 Γ 回归 → 不跳点。
4. **verify_relaxed.py**：解析带 Selective dynamics 头的真实 CONTCAR（audit 亦同步
   SD-aware）；每原子**笛卡尔**位移分解为面内 |dr_xy| 与面外 dr_z；模拟弛豫夹具：
   z-only（反演保持）→ PASS，面内漂移 0.187 Å → FAIL（>1e-3 门）。
5. **SOC INCAR**：ISPIN=2 + 每原子 3 分量 MAGMOM（0 0 0.02）保持；官方建议
   **LASPH=.TRUE.**（PAW 球内非球面梯度，5d 过渡金属）与 **GGA_COMPAT=.FALSE.**
   （恢复非立方晶格全对称，磁性/各向异性推荐）加入全部 4 个 INCAR 模板。
6. **预算分阶段**：Phase 0（para 第 0 步）=56，Phase 1（Round-1 最小集
   3×2×56+4×28）=448，**合计 504**；caps.json 字段化；检查器对"≤500"单表 判 FAIL。
7. **检查器独立性与反例测试**：A2 SD 语义文档检查、LASPH/GGA_COMPAT、分阶段预算；
   `run_tests.py` 15/15 PASS，其中 5 个反例（T F F、缺 LASPH、win spin=true、
   500 预算、GGA_COMPAT=TRUE）均被检查器识别 FAIL。

官方核实：VASP wiki LWRITE_SPN 仅 >6.4.2（para 5.4.4 官方不支持 → 第 0 步裁决意义上升）；
.spn 格式经 wannierberri 生产解析器 + wannier90 conv_write_spn 源码确认。
新增 PARA_STEP0_SCHEME.md（第 0 步独立作业方案）。
