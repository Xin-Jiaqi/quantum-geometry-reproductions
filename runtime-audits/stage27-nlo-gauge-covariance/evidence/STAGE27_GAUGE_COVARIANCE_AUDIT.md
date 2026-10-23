# Stage 27 — Campaign 参数静默丢失修复 + 源规范/胞规范协变几何审计

> 状态：**两项完成**。(A) 修复 Stage 26 review item A（campaign-plan 静默丢弃
> `--numerical-baseline` / `--second-derivative-indexing`）；(B) 独立审计
> Haldane 源规范投影算符几何 vs Wannier 胞规范中心协变导数几何（review item B），
> 在完全统一的 Hamiltonian / Wannier 中心 / k 网格 / 展宽 / 张量约定下完成。
> 所有负面结果与历史 Legacy 基准保留（未覆盖任何历史谱文件）。
> PR（stacked on Stage 26 branch）：见发布记录。

## (A) Campaign 参数静默丢失修复

**缺陷**（review 确认）：`cli.py::_campaign_kwargs` 不含 `numerical-baseline` 与
`second-derivative-indexing`；`campaign.py::create_campaign` 无对应参数且
`RunOptions(...)` 构造不转发。CLI 接受标志但冻结 plan 静默使用
`reviewed-hermitian-v1/cartesian`——历史回归 campaign 不可复现，违反 fail-closed 约定。

**修复**（nlo_isolated，全部通过 41/41）：
- `cli.py::_campaign_kwargs`：加入两个字段（CLI→create_campaign）。
- `campaign.py::create_campaign`：新增 `numerical_baseline`（默认
  `REVIEWED_NUMERICAL_BASELINE`）与 `second_derivative_indexing`（默认 None →
  运行时 build_run 自动解析），转发进 `RunOptions`；`asdict(options)` 自动落入
  冻结 plan 的 `jobs[].options` 与 `_job_fingerprint`。
- **fail-closed**：reviewed 基线 + 显式 legacy 索引在 create_campaign 阶段即抛
  `ValueError`（build_run → NLOConfig.validate 交叉检查），不会静默降级。
- 新测试（test_campaign.py 3 个）：显式 legacy 标志存活到冻结 plan 的 options 且
  指纹与 reviewed 计划不同；非法组合在计划期被拒；CLI `_campaign_kwargs` 含两标志。
- review item C：`_generalized_derivative(indexing=...)` 未知字符串现在显式
  `raise ValueError`（不再静默回退 cartesian），新测试覆盖。

## (B) 源规范 vs 胞规范协变几何审计（`validation/stage27_gauge_covariance_audit.py`）

**动机**：Stage 26 正确撤回了 Stage 25 无条件的 `SM75 = −(π/2)·raw_dressed` 声明。
其投影算符 oracle（σ = −π·R_model_plain）只关于指定的 plain 连接，不是
diagonal-TBA 物理响应的证明。review 要求：用 `τ` 相同的规范协变导数
`D_a P_cell = ∂_a P_cell + i[P_cell, τ_a]` 重做，并把残余定位到具体对象。

**严格统一约定**（全部比较共用）：
- H_src：解析 2 带 spinless Haldane（M=0.4, φ=π/2, t1=1.0, t2=0.2, NN=(3,2), NNN=(3,2)）；
  H_cell = U†H_srcU，U = diag(1, e^{−ik·δ1})，δ1=(1,0) Å。H_cell 与合成包 HR 核
  逐点一致（G0: 7.9e-16）。
- τ_a = i(∂_a U)U†（笛卡尔 k 单位），与包 wannier90-wcenter.dat 的中心差逐点一致
  （G0b: diag(0,1,0) Å, 偏差 1.7e-9）。
- k 网格 5×5×1（收敛 9/13）、η_sp=0.05、η_band=0.005（扫描 0.05…0.002）、
  分量 xxx/xyy/yxy、(bc) 对称化、f_nm·δ 展宽、−2π/Nk·Im C_sym——与 Stage 26
  oracle 完全一致。投影算符 FD：h1=1e-4（一阶/协变）、h2=3e-3（二阶）。

**结果**（`STAGE27_GAUGE_COVARIANCE_AUDIT.json`，全部 PASS）：

| 检查 | 内容 | 值 |
|---|---|---|
| G1 | D_a H_cell == U†∂_a H_src U | 9.6e-9 |
| G2 | D_a P_cell == U†∂_a P_src U | 1.2e-8 |
| G3 | σ(P_src, ∂) == σ(P_cell, D)（规范协变） | 7.0e-5 |
| G4 | σ(P_cell, ∂) ≠ σ(P_cell, D)（plain 是错误对象） | 1.38 |
| G5 | σ(P_cell, D) == **−π·R_dressed**（NLO 生产 raw = 协变响应） | 3.2e-4，常数 3.1406 |
| G6 | σ(P_cell, ∂) == −π·R_plain（Stage 26 plain 对，规范依赖） | 6.3e-4 |
| G7 | R_dressed == R_src（生产 raw ≡ 源规范 raw） | 1.4e-15 |
| G8 | σ(P_src, ∂) == −π·R_src | 3.5e-4，常数 3.1407 |
| G5b | η_band 稳定性（max 0.0059 @η=0.05，否则 ≤6.2e-4）；k 网格 5→9→13：3.2e-4→2.0e-4→1.7e-4 | 收敛 |

**配对表**（所有 σ 几何 × −π·raw 组合，诚实报告）：
- 匹配对（比值 ≈ 1.000）：σ_src↔R_src/R_dressed、σ_cell_cov↔R_src/R_dressed、
  σ_cell_plain↔R_plain。σ_src 与 σ_cell_cov 都匹配生产 raw（R_dressed），因为
  R_dressed ≡ R_src 且规范协变。
- 混合对（比值 ≈ 0.53 或 2.16）：σ_src↔R_plain、σ_cell_plain↔R_src/R_dressed、
  σ_cell_cov↔R_plain——任何跨几何配对都不成立。

**关键物理结论**：
1. **规范协变建立**：源规范（平坦连接，普通 ∂）与胞规范协变导数 D 给出相同的
   投影几何响应（G3, 7e-5，FD 极限内）；D_aP_cell = U†∂_aP_srcU 恒等式数值成立
   （G1/G2 ~1e-8）。
2. **生产响应的定位**：NLO diagonal-TBA 生产 raw **恰好是**规范协变（=源规范）
   响应：σ_cell_cov = −π·R_dressed，测量常数 3.1406≈π，残差 3.2e-4（FD 极限）。
   Stage 26 的"负面发现"被本地化：它比较了错配的规范对象（plain-∂ 胞投影 σ 对
   dressed raw——配对表比值 2.16，且 −π/2 从来不是正确常数）。
3. **plain 胞规范对是自洽但规范依赖的伪像**：σ_cell_plain↔R_plain 内部一致
   （G6），但与源/协变响应相差 1.38（G4）；普通 ∂_aP_cell 不是物理匹配对象。
4. **−π/2 关系不普适且未假设**：测量常数对所有匹配对都是 ≈π（3.1406/3.1407），
   不是 π/2；Stage 25 的 −π/2 声明明确不作假设（π/2 与 π 的差异属于 σ 组装
   （−2π/Nk·Im）与 raw 组装（iabc+iacb）之间的约定因子问题，SI/2D 归一化另属
   未决项）。

**保留的未决项**：SI 绝对前置因子、2D 片层归一化（Lz×0.1）、bath/散射协议、
DFT/W90Flow 认证、FDR SM Eq. (75) 原文不可得（`SKIPPED_NO_FDR_ARCHIVE`；
投影形式为 arXiv:2409.16358 Eqs. (3)+(6) 的独立重构，非 SM 逐字文本）。

## 发布记录

- 分支：`fix/campaign-indexing-flags-stage27`（基于 Stage 26 分支
  `fix/shift-current-cartesian-derivative-stage26`，stacked PR，diff 仅含 Stage 27
  变更；未触碰科学仓库与用户 dirty 工作区）。
- PR：见 GitHub（title 含 Stage 27 标识），base = Stage 26 分支；未合并、未推 main。
- 历史 Legacy 光谱（nlo-run-output-legacy*）与 Stage 26 谱文件未改动。
