# Stage 8l | Stage8k 独立 WAVECAR 证据审计与最省资源的下一门槛（2026-10-22）

> **零新增 VASP/PBS 作业。** 本报告独立交叉核验 commit 39639ac/41f427e 中的解析器、归档 JSON 与前期 Stage8i JSON。远端原始 WAVECAR（约 4.9GB）由 DeepSeek harness 在 Japan IMR 读取；本轮 ChatGPT 仅核对 GitHub 中的源程序和机器证据，未声称再次打开 HPC 原始 WAVECAR。历史 JSON 保留不覆盖；原 Jia PRL 复现不变。

## 1. 真实进展：严格 Sz 守恒的独立必要条件已受否定性检验

harness 已从头解析 VASP 6.4.3 ncl WAVECAR，报告记录跨度/总文件大小一致、目标 k 的 WAVECAR 能量与 EIGENVAL 逐带相等、grouped-up/down 格式和真实 K/−K 时间反演关系通过。

在 PAW 变换对 spin 不混合、广义重叠也分 spin block 的通常前提下，如果物理哈密顿量真正满足 [H,Sz]=0，则相应广义本征值问题分成互不耦合的 up/down 子块。非简并本征态只能位于其中一个块，伪波函数的横向 spinor overlap 必为零。可信非简并横向伪自旋非零，因此足以否定严格 Sz 守恒这一全局 no-go 的前提。

**术语纠正：PAW 变换一般不酉。** 不能像先前报告所说用“酉等价”证明；应通过自旋对角的广义哈密顿量和重叠矩阵的 block structure 论证。数值伪自旋不等于完整 PAW 物理自旋幅度。即便已否定严格守恒，也不能推出静态 NMEE 非零。

## 2. 七条明确的 band35/36 归档错位

把新独立证据 STAGE8K_INDEPENDENT_SPIN_EVIDENCE_20261022.json 中误名为 E_band35_eV 的字段及对应 S_grouped，与旧 STAGE8I_SPIN_TEXTURE_20261022.json 中明确标注 band35/band36 的能量与三个 spin 分量交叉比对，发现 **7 个非简并来源位置中，新文件保存的实际是 band36 而不是 band35**：能量偏差在 4 μeV 内，spin 分量差在 5×10^-6 内。第八个中点 Γ–M 内点 band35/36 精确简并，单凭 E 无法唯一定义带号和该单带自旋；不得强行匹配。

| AB 实际 k 点 | 新文件保存、实际属于 band36 的横向伪自旋模长 | 旧 band35 的横向伪自旋模长 | 两带间隔 |
|---|---:|---:|---:|
| (−1/3,+1/3) Γ–M 内点，非 K′ | **0.088912** | **0.167801** | 7.84 meV |
| (0.476190,0) 近 M | **0.164810** | **0.188710** | 10.285 meV |
| (0.142857,0.142857) Γ–K 内点 | **0.037573** | **0.033579** | 71.252 meV |
| (1/3,1/3) 真 K | 0 | 0 | 87.543 meV |

这不是 WAVECAR 原始记录解析失败的证据，更像生成摘要时将 Python 0-based 35 写作物理 band35（应为 band36）。独立归档的三个 AB **band36** 状态也具有明显非简并横向伪自旋，因此主要微观结论仍获支持；但报告宣称的“三条 band35 原始独立核验”需先重导出真正 zero-index=34 的 spin vector 才算得到对应的机器证据。

零作业核验程序 [stage8l_crosscheck_band_ledger.py](stage8l_crosscheck_band_ledger.py) 与独立回归测试 [test_stage8l_band_evidence_consistency.py](../tests/test_stage8l_band_evidence_consistency.py) 将阻止这种数据错位被再次接受。不得修改旧 JSON 以伪装从未出错。

## 3. 两处工程验收改进

- 原自研 reader 未拒绝不支持的 rtag，原单测通过捕获自己抛出的 AssertionError 可能错误通过。已收紧只支持本次 VASP ncl complex64 格式 rtag=45200，补充 spin/recl/nplane 安全检查，并纠正误把 grouped spinor 写成 interleaved 的源码注释。
- 原分析脚本硬编码远端 Python import 路径且允许 nearest k 自动替代不在 21×21 网格的 exact M / Lambda。已改成相对路径导入和单独报告是否为 exact requested k，防止把附近采样点冒称高对称点或真正 CBM。

以上变化只影响后处理，不修改任何已完成的自洽波函数或能量。

## 4. 更省算力的下一门槛：已有 WAVECAR 后处理，优先于重做两个 SOC SCF

VASP 官方 [LORBIT 页面](https://vasp.at/wiki/LORBIT) 明确指出 LORBIT orbital projection 是 Kohn–Sham orbitals 的后处理，已收敛 WAVECAR 可以设置 ALGO=None 与需要的 LORBIT 直接重启；[ALGO 页面](https://vasp.at/wiki/ALGO) 也说明 None/Nothing 为只读既有轨道的 postprocessing。

**建议方案，尚未授权任何新作业：**
1. 对现有许可的 VASP6 非共线 SOC、相同 spinor WAVECAR、KPOINTS、W+Se PAW 和原 INCAR 作现场预检。在新目录复制/以只读方式接入旧波函数，别在原始作业目录覆盖输入输出。首先只对 AB 做最小 LORBIT=11 投影试算，建议初值 ISTART=1, ALGO=None, LORBIT=11, LWAVE=.FALSE., LCHARG=.FALSE.；必须核验站点版 VASP 非共线 SOC 与 CHGCAR 重启兼容、是否需要 ICHARG=11、原始文件保全和 PBS 资源配额。
2. 如果 AB 能获得真实完整 PROCAR，再决定是否在单独授权范围内对 glide 做同样的后处理。出现异常要追查前两次 Exit_status=255 的 MPI/I-O 原因；**两次完整 SCF 从默认方案降为最后备选**。
3. PROCAR orbital/site projection 能辅助 AB 与 glide 的 band-character matching，但它**仍不是完整且唯一的空间层 PAW 自旋算符**，interstitial 及 onsite PAW 增广需要另外处理。不得从 PROCAR 对角投影推导 off-diagonal spin–velocity Kubo 响应。
4. 若需要检验物理自旋-速度带间机制，可以先用现存 WAVECAR/PAW onsite 数据设计少数 k 和少数 occupied/empty bands 的 operator-complete feasibility 原型，再考虑全 BZ Wannier 插值。

目前新任务并不包括新的 qsub，也不自动批准 Stage7j Wannier 或真实 NMEE 计算。

## 5. 科学边界与下一关键可证伪命题

- 已有 Stage8h：受限弛豫的 WSe2 half-slide 结构保留 nonsymmorphic layer-exchanging glide，total equilibrium/static Γ 被禁止，layer-staggered Γ 仅被对称性允许。
- Stage8i/8k：正的 **采样** SOC 能隙，非简并伪 spinor 横向混合支持 AB 的严格 Sz 守恒失效；但完整 PAW quantitative spin、正确 band35/36 归档、实际 spin and velocity matrix vertices 以及任何真实 WSe2 static Γ 均未完成。
- Stage10：hidden layer/ferroelectric Rashba/spin quantum geometry/intrinsic insulator NMEE 已有先例，只有新物理的真实材料数值与可靠机制才可能构成论文主张。

**当前状态：PSEUDO_NONCONSERVATION_SUPPORTED__ARCHIVE_BAND_LABEL_ERROR_IDENTIFIED__PAW_LAYER_OPERATOR_PENDING__MATERIAL_STATIC_NMEE_NOT_CALCULATED。**
