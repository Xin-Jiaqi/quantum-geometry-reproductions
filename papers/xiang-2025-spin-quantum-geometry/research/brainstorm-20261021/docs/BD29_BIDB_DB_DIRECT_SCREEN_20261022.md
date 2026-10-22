# BD29 | 真实 BiDB SQLite 筛选：Pd2Te2 意外成为最优群论试验台

**2026-10-22。已独立读取用户提供的完整 BiDB ASE/SQLite 数据库（约 48 MB），只读、结构和 metadata 均已实际检查。没有提交 DFT 或 HPC。**

## 数据状态与提取方法

数据库通过 PRAGMA quick_check = ok；11,184 个结构系统，992 单层和 10,192 双层。字段 slide_stability 中严格标注 Stable 的有 2,736 个，另有 Stable-flat 30 个（这些统计不要误同论文另一个 2,586 条数据集的筛选定义）。BiDB 的 stoichiometry 写成 AB2/ABC 等通式，故筛选代码使用 species 中原子序数/数目并按最大公约数归一化，避免 CrI2 与 CrI3 等假阳性。条件 number_of_layers=2 且 slide_stability=Stable。

数据来源：[BiDB 官方字段文档](https://cmr.fysik.dtu.dk/bidb/bidb.html)及 [Pakdel et al., Nature Communications 15, 932 (2024)](https://doi.org/10.1038/s41467-024-45003-w)。官方字段单位：binding_energy_gs 为 meV/Å²；ehull 为 eV/atom；gap_pbe 为 eV。2DHub 表明数据库为 CC BY-SA 4.0。只公开本次选定 UID 与少量派生统计，不公开用户上传的 48 MB 数据库原件。

| 化学式 | Stable 双层数 | 代表原胞 | 关键统计 |
|---|---:|---:|---|
| VTeS | **2** | 6 原子 | 两个 Stable 均 p3m1/无反演、BiDB 层间 FM，PBE 金属 |
| VSeTe | **8** | 6 原子 | p-3m1（有反演）与 p3m1（无反演）均存在，8 条层间 FM、金属 |
| AgVP2Se6 | **17** | 20 原子 | 14 条层间 FM、3 条 AFM，PBE gap 约 1.08–1.22 eV |
| MoSSe | **7** | 6 原子 | p3m1 和 p-3m1 双层 |
| MnPS3 | **7** | 20 原子 | BiDB 记录层间 AFM，不可替代真实层内 Néel 序 |
| WSe2 | **4** | 3 条 6 原子、1 条 12 原子其他多晶型 | 6 原子含有/无反演两类 |
| CrI3 | **3** | 16 原子真实 Cr2I6 双层 | 2 条层间 FM、1 条层间 AFM |

数据库中的层群是**结构**层群，不是带 SOC 与真实磁序的 magnetic layer group；层间 FM/AFM 不保证层内真实磁基态。不同 UID 的取向操作可能不同，不可以把不同堆叠都当成**同一相对方向上的纯滑移路径**。

## 更高性价比的意外候选：真实 m_z 层交换对称性

BD27/BD28 早先的六方群论例子严格假定 C3z、竖直镜面 mv 与纯水平层交换镜面 mz 同时成立。原先优先考虑的 VTeS 与 WSe2 的某些稳定结构不满足 mz。真正满足该结构类型的 BiDB Stable 双层，按 layer_group=p-6m2（D3h）筛选**只有 7 条**：

| BiDB ID | 家族 | 原胞 | gap PBE | DB magnetic | 来源单层 ehull（eV/atom） |
|---:|---|---:|---:|---|---:|
| **980** | **Pd2Te2** | 8 | 0 | 0 | **0.13369** |
| **975** | **Pd2Se2** | 8 | 0 | 0 | **0.18294** |
| 2903 | NBa2 | 6 | 0 | 0 | 0.16450 |
| 3555 | ScSe2 | 6 | 1.590 | 1 | 0.04652 |
| 1350 | Cu2S2 | 8 | 0 | 0 | 未逐项核对 |
| 1463 | Ni2Se2 | 8 | 0 | 0 | 未逐项核对 |
| 1478 | Ni2S2 | 8 | 0 | 0 | 未逐项核对 |

**Pd2Te2 完整来源 UID：** Pd2Te2-51a4897edf43-2--1_0_0_-1--0.67_-0.33（id 980）。非磁性标记、p-6m2、无反演、PBE gap=0、binding_energy_gs=47.124629 meV/Å²，slide_stability=Stable。其来源单层 dynamically_stable=Yes，**但是 ehull=0.134 eV/atom**：这仅支持用作计算/理论模型，不能宣称热力学稳定或易合成。

**Pd2Se2 完整来源 UID：** Pd2Se2-a3734aabd910-2--1_0_0_-1--0.67_-0.33（id 975）。类似非磁性 D3h 金属、8 原子，binding_energy_gs=50.10556 meV/Å²，单层 ehull=0.183 eV/atom。

### 实际原子位点几何置换测试

从 ASE SQLite 的 numbers、positions、cell 二进制字段重建真实原子坐标，严格进行元素保持的一对一原子匹配、二维 Bravais 平移模商。对 p-6m2 的纯 horizontal mz 映射，确实检查上下层交换，而非仅复制 BiDB 层群标签：

| BiDB id | C3z | mv | **纯 mz** | 最大成功匹配残差 |
|---:|---|---|---|---:|
| **980 Pd2Te2** | 8/8 PASS | 8/8 PASS | **8/8 PASS** | 0.001291 Å |
| **975 Pd2Se2** | 8/8 PASS | 8/8 PASS | **8/8 PASS** | 0.001044 Å |
| 2903 NBa2 | PASS | PASS | PASS | 0 Å |
| 3555 ScSe2 | PASS | PASS | PASS | 0 Å |
| 2684/2685 WSe2 | PASS | PASS | **mz FAIL** | 首次失败位置约 1.916 Å |
| 8055 VTeS | PASS | PASS | **mz FAIL** | 约 0.673 Å |
| 7964 VSeTe | PASS | PASS | **mz FAIL** | 约 2.005 Å |

完整的可复验代码在 [Group Theory Toolkit Draft PR #32 的 BiDB symmetry gate](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/bidb_geometry_layer_swap_gate.py)。该代码不含数据库本体。FAIL 只排除检查的**纯 horizontal mirror**，不等于排除所有复合 Seitz 操作。PASS 仅是结构几何约束，**不是**磁空间群、SOC/Kubo 验收。

## 材料与理论两条线重新排序

1. **BD27/BD28 非磁性构型空间响应理论：Pd2Te2 首选，Pd2Se2 负/对照。** 8 原子、真实水平层交换、p-6m2、PBE 金属性，无须先以磁性 SCF 和载流子门控起步。需要证明真正的 current-induced spin susceptibility 非零；Pd/Te 较大 SOC 只是选材动机。单层 e-hull 较高，所以定位为**模型验证**而非实验材料预言。
2. **磁性推广：VSeTe 优先于 VTeS**（8 vs 2 个 Stable 记录，均 6 原子且金属）。Janus 的层交换、FM 磁化轴与 antiunitary operations 必须另外识别。不能把非磁性 D3h 规则直接套在 p3m1 FM 的结构上。
3. **较复杂极性磁性：AgVP2Se6** 有 17 个 Stable FM/AFM 层间状态和较小磁交换差，但 20 原子、低对称磁群与竞争非共线磁态使其不宜作为第一项 Kubo 材料计算。
4. **已发表 LEE baseline**：MoSSe/MoTe2/WTe2 的公开工作可用来复现响应算符；WSe2 为 SOC 强的 6 原子备选，但目前检查的 stable p-3m1/p3m1 结构**没有纯 mz**，不可冒充 D3h 试验台。

## 下一步只做两个低成本 Gate

**A. 物理 Hamiltonian 认证：** 以 Pd2Te2 id980 的原子/轨道基建立能体现原子 SOC 的最小源支撑模型，先过完整 orbital-site 对称性协变。只有确认非零 spin-current Kubo/Boltzmann 的模型才值得研究配对 F_minus(delta)、F_plus(delta) 纹理。无 DFT 自动授权。

**B. 公式同层级竞争核对：** 已有 2026 LEE PRL 直接覆盖 LEE 存在与固定堆叠的允许张量；2023 BSF 直接覆盖层群构型。我们潜在的高层次理论要研究**整个堆叠空间上的算符依赖零点/最低开启阶/纹理**。某个其他材料中已有 Edelstein，并不自动否定该普适分类；但一个更一般的已发表定理若已经蕴含这些具体结果，则会形成真正的直接竞争。仍需逐式核对 LEE SM。

**不需要用户再提供别的数据库或额外导出 CIF。** 已在私人工作环境导出选定 Stable 行 CSV 与 11 个代表结构 POSCAR/manifest（仅提供用户，不上传公开 GitHub）。MnPS3 R1 CLOSED；R2-P0 READ-ONLY；无新 DFT、Slurm、WSL、Japan HPC。
