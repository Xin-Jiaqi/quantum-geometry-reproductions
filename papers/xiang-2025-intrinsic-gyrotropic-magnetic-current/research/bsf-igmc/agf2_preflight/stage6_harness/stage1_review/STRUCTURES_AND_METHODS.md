# AgF₂ 双层 B1/B2 结构与方法审查（Stage 1）

> 范围：Zhu et al., *PRL* **135**, 056801 (2025)（DOI [10.1103/dmzg-ck2t](https://journals.aps.org/prl/abstract/10.1103/dmzg-ck2t)）AgF₂ 双层滑移铁电体系 —— 文献与计算输入审查。
> 结论先行：**BLOCKED_SOURCE / BLOCKED_PROTOCOL**。作者弛豫后的 B1/B2 原子坐标**不公开**（PRL 数据可用性声明："data are available from the authors upon reasonable request"；SM 不含任何原子坐标）；本阶段按 HARNESS_TASK A.3 的授权，基于**公开可靠的 C2DB 单层**（2AgF2-1，CC BY-SA 4.0）构造了 B1/B2 候选结构（**明确标注：构造的待弛豫候选，非作者弛豫结构**），并完成对称性/磁对称性验证。任何 DFT 结果均未运行（`data_origin=NOT_RUN`）。
> 修订记录（按 PR #8 审查意见）：①**移除全部版权文件**（原文 PDF/TeX/原图/SM 全文仅作本地私有归档，见 `sources/README.md`）；②**设置核验统一为单一全局坐标变换**（S, δ）同时映射全部群操作（拒绝逐操作原点位移）；③**磁层群标签按母层群消歧**（uni 33 同时对应 p11a' 与 pb'11，按已验证母群 LG12 选 pb'11，给出轴设置/原点/Seitz 对应）；④**k 网格改为显式取整 14×9×1 初始候选**并记录 Gamma/MP 与 ISMEAR/SIGMA 收敛计划；⑤**来源记录强化**（C2DB 行校验哈希、转换算法、许可原文、工具包来源统一为用户仓库 v0.18.0）；⑥3.0 Å 层间间隙定位为**未收敛初猜**并给出多初猜弛豫方案；⑦作者坐标索取计划（双来源分支比较）。
> 二轮修订（按二轮独立审查）：⑧**轴矢量连接变换修复（BLOCKER 1）**：实自旋 m′=(−1)^τ·det(R)·R·m；Type-I = 层间 Néel 反平行（幺正 m₀₀₁，面内 S_x→−S_x）、Type-II = Néel 平行（反幺正 m′₀₀₁，S_x→+S_x），替换旧极矢量规则与颠倒的赋值；12/12 逐原子映射表 + 独立 stdlib 回归测试；⑨**Git 历史清理（BLOCKER 2）**：分支重写为基于 main 的单一干净提交，含版权 blob 的旧提交不再可达、本地对象库清除；⑩k 网格奇偶表述修正（14×9×1 中 14 为偶数，Gamma-centered 可混合奇偶）；⑪脚本来源路径可配置（`GTOT_TOOLKIT_SRC`/`AGF2_C2DB_ROW_JSON`/`AGF2_PRIVATE_ARCHIVE_DIR`）与 offline/local_data_required 分工。

---

## 1. 结构来源与处理流程

| 文件 | 内容 | SHA-256 |
| --- | --- | --- |
| `structures/C2DB_2AgF2-1_monolayer/POSCAR` | C2DB 单层 AgF₂（UID 2AgF2-1，LG17 p2₁/b11） | `42468f68fc87912a365df505ab7e73952df88fab574be74c31bb077c0686fa7c` |
| `structures/C2DB_2AgF2-1_monolayer/AgF2_monolayer.cif` | 可视化 CIF | `e55ffb7a8bd2747e91cd2a77d5d0ff72700a095a2fb2ecd086336f2a65360a6d` |
| `structures/constructed_B1_candidate/POSCAR` | **构造候选** B1（未弛豫） | `b3c1b942df69b5442e52d0d38fe839b4b61feff3619b2adb4736faa82501c122` |
| `structures/constructed_B2_candidate/POSCAR` | **构造候选** B2（未弛豫） | `ab2f3d84af0a074079824b1f0b94d928d16609e889ac7b1ba3fd9a50e2c69087` |
| 各目录 `provenance.json` | 来源与构造记录 | 见 `results/STRUCTURES_MANIFEST.json` |

**单层来源（原文给定之外的公开权威数据）**：C2DB（Computational 2D Materials Database，[cmr.fysik.dtu.dk/c2db](https://cmr.fysik.dtu.dk/c2db/)），UID `2AgF2-1`，注册表 `C2DB_all_data.csv` 第 4194 行 ↔ 数据库系统 4195（ASE sqlite：`project-cy/c2db.db`，2024-09-26）。单元格 a=3.82702 Å、b=5.87967 Å、c=31.52656 Å（z 真空），Ag₂F₄（2 个 Ag、4 个 F），PBE/GPAW 弛豫。已在用户 2022-11-30 旧库交叉验证（system 10217，晶胞与能量 E=−19.3594 eV 一致）。许可证：CC BY-SA 4.0（见用户资料库 `DATA_LICENSES.md` 注册表）。

**双层构造规则（论文给定）**：B = S + ÔS，Ô = {m₀₀₁ | t₀}（SM §1.2 约定：算符符号遵循 Bilbao，m₀₀₁ = 法向沿 [001] 的镜面，绕层间中面 z_m = Lz/2 作用；t₀ 为平面内分式平移）：
- B1：t₀₁ = (0.5, −0.16)；B2：t₀₂ = (0.5, +0.16)
- 顶层 = 底层经 m₀₀₁（绕 z_m）镜像后叠加平面内平移 t₀；两层共用一个平面内晶胞（论文：与单层相同横向晶胞）
- **假设（未发表，明确标注）**：层间 F–F 间隙 g = 3.0 Å（vdW 接触估计）→ d(Ag–Ag)=4.969 Å；真空 20 Å（论文协议）；POSCAR 物种块顺序 Ag×4 + F×8（VASP-5 规范）

**构造验证（全部通过，单一全局坐标变换核验）**：
1. 单层 2D 层群 = **p2₁/b11（LG17）**：spglib 检测的 4 个操作在**同一全局原点平移 δ=(0,0,0)** 下同时映射到工具包 LG17 标准操作集合（E、{−I|0}、{2_x|(½,½,0)}、{m_x|(½,½,0)}）；直接原子映射 4/4 精确成立；spglib 3D 分类 P2₁/c（#14）= LG17 的 3D 实现。4 个操作与 C2DB/论文的连接算符 {m₁₀₀|(1/2,1/2)}、{2₁₀₀|(1/2,1/2)} 直接对应（C2DB 原点 = 工具包标准原点）。[注：z 按计算胞周期比较（真空区对称），与 spglib 3D 分类一致]
2. 双层 B1、B2 的 2D 层群 = **pb11（LG12）**：spglib 检测的 2 个操作在**同一全局原点平移 δ=(0.75,0,0)（≡ −0.25 mod 镜面核）** 下同时映射到工具包 pb11 标准操作集合（E、{m_x|(0,½)}）；直接原子映射 2/2 精确成立；spglib 3D 分类 **Pc（#7）**，与论文"p2₁/b11 → pb11"一致 ✓。
3. 连接算符：B2 ≡ {m₀₀₁ | (0.5,0.16)} B1（mod 晶格），与 SM 式 (3)（t₀₂ = −N⁻t₀₁ + t₀，N⁻=m₀₀₁）一致 ✓。

## 2. 磁对称验证（工具包 + spglib 2.5.0，`results/magnetic_symmetry_report.json`）

| 体系 | 自旋构型（输入假设） | 磁性层群 / 磁点群 | 与论文一致性 |
| --- | --- | --- | --- |
| 单层 | 反铁磁排列：Ag(0,0)↑x̂、Ag(½,½)↓x̂（共线，沿 [100]） | **p21'/b'11（uni 86 = P2₁′/c′）**，磁点群 **2′/m′**，含反幺正镜面 **m′₁₀₀** ✓ | C2DB 声明 MSG P2₁′/c′ ✓；预检源磁群 m′_x（反幺正镜面存在）✓ |
| B1/B2 Type-I（层间 Néel 反平行，0% 应变基态） | 顶层镜像位点磁矩与底层**相反**（Néel 矢量反平行；m₀₀₁ 连接的原子自旋相反） | **pb'11（母群 pb11，type III，MPG m′）**，磁点群 **m′ = {1, m′₁₀₀}**（与 Type-II 相同） | 论文："B1(B2) has the magnetic point group m′ = {1, m′₁₀₀}" ✓；**设置对应**：uni 33 同时对应 p11a'（母群 p11a）与 pb'11（母群 pb11）——按已独立验证的母群 LG12 消歧选 **pb'11**；在 δ=(0.75,0,0) 全局原点平移下，spglib 磁性操作 {E|0}、{m_x|(½,½,0)}·T 映射为工具包 pb'11 标准操作 {E|0}、{m_x|(0,½,0)}·T（反幺正镜面 ⊥ a、b/2 滑移）✓ |
| B1/B2 Type-II（层间 Néel 平行，4% 应变基态） | 顶层镜像位点磁矩与底层**相同**（Néel 矢量平行；m₀₀₁ 连接的原子自旋相同） | 同上（内部群同为 m′；连接算符 m′₀₀₁） | 论文 Type-II 连接算符为 m′₀₀₁ ✓（SM："the magnetic atoms connected by N⁻ = m₀₀₁ have the same spin"） |
| 连接算符（磁，轴矢量规则） | — | Type-I：{m₀₀₁|t₂} **幺正**——实自旋为轴矢量，m′ = (−1)^τ·det(R)·R·m，R=diag(1,1,−1) → 面内 S_x→**−S_x**（连接位点自旋相反）✓；Type-II：**反幺正 m′₀₀₁**（τ=1）→ S_x→**+S_x**（连接位点自旋相同）✓ | 论文 Fig. 2(b)/(c) ✓；12/12 原子逐位精确映射（`results/magnetic_symmetry_report.json` per_atom 表） |

**物理注记（非阻塞，含二轮审查修复）**：单层完整磁点群 2′/m′ 含**幺正反演** → 单层 IGMC C 与 AHE σ_xy 在对称性下为零（与论文"单层铁电被反演禁止、AHE 为零"一致）；预检文档以 m′_x 子群给出 C 的允许分量，属子群必要性判据。**双层的价值在于反演被堆叠打破**（pb11 无反演），C/σ_xy 可在双层非零并与 P_z 联动切换——这是后续 Stage 4+ 的物理主线。
**二轮审查修复（轴矢量连接变换，BLOCKER 1）**：上一版把 Type-I/II 的磁矩赋值与变换写反（用极矢量 z-翻转规则得到假阳性 PASS）。按论文定义（主文 Fig. 2(b)/(c) 与 SM §1）：**Type-I = 层间 Néel 矢量反平行**、由幺正 {m₀₀₁|t₂} 连接（轴矢量规则使面内 S_x→−S_x，连接位点自旋相反；开关 k 空间自旋极化、不开关 AHE）；**Type-II = Néel 矢量平行**、由反幺正 {m′₀₀₁|t₂} 连接（S_x→+S_x；开关 AHE、不开关自旋极化）。内部磁点群 m′={1,m′₁₀₀} 对两类均成立（与论文一致），与外部跨态连接核验分开报告。新增独立回归测试 `agf2_preflight/stage6_harness/tests/test_magnetic_connector_axial.py`（纯 stdlib，CI 可运行）：非共线探针 (0.3,0.7,−0.5)、12/12 逐原子映射、以及"交换 τ 或沿用旧极矢量规则必须失败"的判别力断言。

## 3. 参数核对（详见 `PARAMETER_VERIFICATION.md`）

**原文给定（SM §2 逐字提取，作为 source-derived baseline）**：VASP+PAW、PBE、GGA+U（Dudarev）U_eff(Ag-4d)=4.0 eV [Derzsi et al. PRB 105, L081113 (2022)]、DFT-D3（Grimme 2010）、SOC 开、ENCUT 500 eV、k 间距 0.02×2π Å⁻¹、真空 20 Å、离子力 <10⁻² eV/Å、电子收敛 10⁻⁶ eV、Berry 相极化、CI-NEB 势垒、Wannier90+WannierTools。
**仍然缺失**：PAW POTCAR 版本（SM 未给出；本机 WSL 检查：无 vasp 二进制、无 POTCAR）、整数 k 网格、ISMEAR/展宽、MAGMOM 初值、偶极修正与否、作者弛豫后的精确晶格常数/层间距、作者单层坐标。
**操作选择（operator choices，非作者参数，需批准+收敛记录）**：整数 k 网格（初始候选 14×9×1）、ISMEAR/SIGMA、MAGMOM 初始化、偶极修正、层间 F–F 间隙（初猜 3.0 Å，未收敛）。

## 4. 计算计划与资源估计（未执行，仅供评审）

| 计算 | 体系/原子数 | k 网格（估计） | 备注 |
| --- | --- | --- | --- |
| 单层参考 | 6 原子 | 14×9×1（显式取整，见下） | 校验单层 LG/磁群与能带 |
| B1、B2（Type-I，0% 应变） | 12 原子 | 14×9×1 | 弛豫 + 极化 + 势垒（CI-NEB 8–16 像） |
| B1s、B2s（Type-II，4% 双轴应变） | 12 原子 | 14×9×1（应变后重新取整） | 应变晶胞（+4% a,b），弛豫 + 磁基态比对 |
| 层间距初猜扫描（几何试算，CONDITIONAL-GO） | 12 原子 | 12×8×1（粗） | ≥2–3 个不同初猜层间距，非 SOC、固定面内晶格；检查最短跨层原子间距、残余力、磁序与连接对称 |
| 极化/势垒 | — | — | Berry 相极化；NEB 沿滑移路径（B1↔B2） |
| IGMC/Wannier（Stage 4+） | 12 原子 | 更密网格 | Wannier90 投影 + 自旋/位置算符链审计 |

**k 网格显式取整策略（审查意见第 4 条）**：Δk = 0.02×2π Å⁻¹；N_i = floor(2π/(a_i·Δk)) + 1：
- a=3.827 Å → N_a = 13.06 → **14**；b=5.880 Å → N_b = 8.50 → **9**；z 方向 1（真空）。
- 初始候选 **14×9×1**（Gamma 中心；14 为偶数、9 为奇数——Gamma-centered 网格可混合奇偶，奇偶性作为操作选择记录）；做 k 收敛扫描（12×8、14×9、16×10），记录 Gamma vs Monkhorst–Pack、偶数/奇数网格、ISMEAR/SIGMA 的选择与收敛数据后冻结。

**资源估计（量级，假设 16–64 核小集群）**：12 原子 + U + SOC + 500 eV + ~126 k 点（14×9），每 SCF 约 2–10 min；每个结构弛豫 50–100 离子步 → 每结构约 200–500 核时；6 结构 + 层间距扫描 + NEB + 极化 + Wannier 总计约 **2,000–5,000 核时**（1–2 周墙钟，32 核并行）。**前提**：需要 VASP 许可与 POTCAR（本机无），以及 Slurm 集群（本阶段按指示不使用）。

## 5. 未决问题（BLOCKED 项清单）

1. **作者弛豫 B1/B2 坐标：不公开**（PRL 数据可用性声明；SM 无坐标）。→ ①候选结构已被用户**有条件接受为"待弛豫试算的起始几何"**（不视为作者弛豫结构）；②按审查意见第 7 条：建议向通信作者索取 B1/B2 relaxed POSCAR/CIF、层间距、PAW/DFT+U、磁矩初始化与 CI-NEB 初终态；收到后**与 C2DB 派生候选作为两个不同来源的分支**并列比较（不覆盖旧 provenance）。
2. **POTCAR 版本与锁：缺失**（SM 未给出；无 VASP 安装）→ `potcar_lock_sha256 = null`，审计给出 `BLOCKED_SOURCE_PROTOCOL`。
3. **整数 k 网格/ISMEAR/MAGMOM/偶极修正：操作选择** → 按 §4 收敛并冻结，记录在收敛记录中。
4. **物理输入未经用户批准**（`physics_approved_by_user=false`，按审查意见**不得自动改为 true**）→ 待用户在计算资源/结构方案明确批准后再改；本阶段 `NOT_RUN` 与 `BLOCKED_SOURCE_PROTOCOL` 是正确的科学状态。
5. **层间 F–F 间隙 3.0 Å 是未收敛初猜**（非冻结材料参数）→ 在合法 VASP 环境可用且预算获批后，先做 §4 的多初猜非 SOC 弛豫，B1/B2 配对采用同一参数，在 SOC/极化/NEB 之前冻结协议。
6. **外部磁连接核验（二轮 BLOCKER 1）**：已按轴矢量规则修复并重验（Type-I 幺正 m₀₀₁ ↔ 层间 Néel 反平行；Type-II 反幺正 m′₀₀₁ ↔ Néel 平行；12/12 原子逐位映射、判别力测试通过）——原 `BLOCKED_AXIAL_SPIN_TRANSFORM` 状态解除；内部磁点群 m′={1,m′₁₀₀} 两类相同、独立报告。**Git 历史清理（二轮 BLOCKER 2）**：分支已重写为单一干净提交（基于 main，仅含最终审查文件）；本地对象库已清除旧 blob；GitHub 侧对旧 PR 引用 blob 的彻底移除需仓库所有者向 GitHub 提交敏感数据移除申请（本仓库无法代为执行）。

## 6. 交付文件索引

- `results/agf2_preflight.stage1.json` —— 预检门（对称性 + 输入合同；结构哈希全部通过；仅 POTCAR 锁阻塞）
- `results/audit_receipt.stage1.json` —— receipt 审计（**BLOCKED_SOURCE_PROTOCOL**，非 PASS）
- `results/monolayer_symmetry_report.json`、`results/bilayer_symmetry_report.json`、`results/layer_groups_2d_report.json`、`results/magnetic_symmetry_report.json`
- `results/SOURCES_MANIFEST.json`、`results/STRUCTURES_MANIFEST.json`、`results/agf2_protocol_frozen.json`
- `stage6_receipt.stage1.json` —— 依模板填写（`data_origin=NOT_RUN`）
- `sources/README.md`、`sources/.gitignore` —— 私有归档政策（版权文件不入库）
- 复现脚本：`scripts/export_c2db_monolayer.py`、`scripts/build_bilayer_candidates.py`、`scripts/gtot_common.py`、`scripts/verify_layer_groups_2d.py`、`scripts/verify_bilayer_symmetry.py`、`scripts/verify_magnetic_symmetry.py`、`scripts/verify_monolayer_symmetry.py`、`scripts/write_source_manifest.py`、`scripts/write_structures_manifest.py`

### 6.1 复现性说明（二轮审查非阻塞项，脚本分工）

- **离线可执行（offline，空白 checkout 可直接运行，输入全部在仓库内）**：`build_bilayer_candidates.py`、`gtot_common.py`、`verify_layer_groups_2d.py`、`verify_bilayer_symmetry.py`、`verify_monolayer_symmetry.py`、`verify_magnetic_symmetry.py`（后四者需要 pinned 工具包，路径可由环境变量 `GTOT_TOOLKIT_SRC` 指定，默认 `/tmp/gtot/toolkit-v0.18.0/src`）、`write_structures_manifest.py`。
- **本地数据依赖（local_data_required，不可从空白 checkout 重建）**：`export_c2db_monolayer.py` 需要 C2DB 本地库行（`AGF2_C2DB_ROW_JSON`，默认 `/tmp/agf2_2AgF2-1_structure.json`）；`write_source_manifest.py` 需要仓库外的私有归档目录（`AGF2_PRIVATE_ARCHIVE_DIR`，默认工作区中仓库根兄弟目录 `private-agf2-stage1-archives/`）。仓库内冻结的 POSCAR/CIF/provenance/清单均由哈希钉定（`results/STRUCTURES_MANIFEST.json`、`results/SOURCES_MANIFEST.json`），空白 checkout 只能**复核哈希**，不能声称重新生成原始 C2DB 行或版权文件本身。
- 独立回归测试 `agf2_preflight/stage6_harness/tests/test_magnetic_connector_axial.py`：纯 stdlib，CI（stage6 工作流）与本地均可运行，无需工具包/numpy。
