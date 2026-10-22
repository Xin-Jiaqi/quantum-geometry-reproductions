# Stage8l | DeepSeek harness 单文件自主交接（仅使用已有 WAVECAR，零 PBS 作业）

**研究分支：** research/bsf-nmee-stage7-sw-gate-20261020，仓库 Xin-Jiaqi/quantum-geometry-reproductions。

本阶段最重要的目标：**把 Stage8k 的 band35/36 证据编号修正，并首次直接筛选已有 real WSe2 spinor WAVECAR 中的占据–空带伪自旋非对角顶角**，为真正物理 PAW + Wannier/Kubo 决定是否值得花资源提供高质量预审。

用户只需口头授权执行这份 GitHub 文档，**不承担下载 4.9GB 文件、逐条 SSH 命令、脚本执行、复制结果和手工 git 提交**。

## 必读来源

- 本目录 [Stage8l 证据更正和低成本计算路线](STAGE8L_EVIDENCE_RECONCILIATION_AND_EFFICIENT_NEXT_GATE.md)。
- Stage8i 原始 [自旋纹理 JSON](STAGE8I_SPIN_TEXTURE_20261022.json)、Stage8k [机器证据 JSON](STAGE8K_INDEPENDENT_SPIN_EVIDENCE_20261022.json)、[raw 解析器](stage8k_independent_reader.py)、[后处理代码](stage8k_spin_analysis.py)。
- 新增 [band ledger 自动校验](stage8l_crosscheck_band_ledger.py) 和 [pseudo occupied-empty spin-only vertices](stage8l_pseudo_interband_gate.py)。
- Stage7j 中已经建立的 orbital gauge、Wilson source-connection 与物理 spin vertices 判据；Stage9 中的严格 Sz 守恒 no-go。

## 自动执行（所有步骤自行完成）

1. 在本地 fetch 最新研究分支并保留不可变 HEAD SHA，读取最新 8l 与历史 8i–8k 证据。**不修改历史报告的数值或覆盖旧 JSON**。
2. 运行 band-archived-JSON crosschecker：新 Stage8k key_results 的命名 E_band35_eV/S_grouped 有 **0-based index+1 错位**，在七个非简并记录实际为 band36，另一个 glide 精确简并不许单带标注。应从现有远端真实 WAVECAR 只读重新分别导出 physical 1-based **band35 (Python index34)** 和 **band36 (Python index35)** 的 E、normalized pseudospin、完整最近能带间隔、k-index、文件 SHA、相位测试；用新版本化 CORRECTED JSON 留下原件并链接 erratum。
3. 通过已有合法 SSH 读取 Japan IMR Stage8i **原来就有的** WAVECAR/EIGENVAL（两个作业 491186/491187.cmspbs1），从现场 resume state 定位实际目录；原始文件 SHA 与 GitHub 归档比对。**绝不**把 WAVECAR、POTCAR、SSH keys 或许可证内容复制到公开 GitHub、聊天或临时外部服务器。
4. 在授权 HPC 的只读数据上使用新增 stage8l_pseudo_interband_gate.py，AB 和 glide 各做一次读取性质筛选，无新 VASP。默认 k 取真实 Γ–M 2/3，近 M 网格点 (10/21,0)，Γ–K 内点 (3/21,3/21)，真 K (1/3,1/3)；窗口 occupied spinor bands **33–36 (1-based)** 与 empty spinor bands **37–44 (1-based)**。严格对比 WAVECAR enocc 与 EIGENVAL 全 band，确认 k-grid EXACT 命中，不允许把最近 k 当精确高对称点。
5. 计算现有伪波函数的跨隙 spin matrices Sx_{vc}, Sy_{vc}，每个非简并 band pair 的幅值及 selected-window aggregate Frobenius norm。必须区分 **pseudo-plane-wave Euclidean-normalized spin matrix ≠ physical PAW spin matrix**；检查独立相位旋转只改变复数矩阵元相位不改变幅值。对于精确简并处，单个 band pair 的幅值具有简并子空间 gauge 依赖，应使用完整子空间块指标并写明 pseudo metric 不满足 PAW physical gauge 的局限。
6. 逐项科学解释：非简并横向 pseudo spin 已支持 AB strict Sz nonconservation，但即便跨隙 pseudo Sx 非零，也**仅**证明有一个 spin-vertex channel，不证明 spin×velocity product/Berry curvature 积分非零。滑移 midpoint 的 total static tensor 仍因真实 constrained geometry glide 严格为零；layer-staggered static Gamma 被允许但没有数值。
7. 仔细处理原解析器检查缺口：raw WAVECAR 格式 rtag、nplane、spin channel 和 record length 不匹配时 FAIL CLOSED；交叉验证新旧脚本的 energy/band indexing 与 selected k. 需要证据时提供 sanitized SHA 与数值比较，不要从旧论文叙事中猜结果。
8. GitHub 仅提交独立研究的新增代码/测试、**正确带号**的净化 JSON、中文审查和 RESEARCH_INDEX/paper-associated navigation，保留 Draft PR #3，冻结 Jia 原 PRL 复现完全不变。
9. 给出未来**单独待批准**的最低成本 PAW orbital projection proposal：VASP 官方允许基于现存 WAVECAR 的 ALGO=None + LORBIT=11 重启后处理。若需要可先做 AB 后做 glide、在新目录保留源文件、核实完整 PROCAR 与合理退出。**不要把它等同于完整物理 PAW layer spin 或 Kubo 顶角，更不要默认完整重做两个 SOC SCF。**

## 提交/权限界限

- 本次 **新增 PBS/VASP/QE/NEB/Wannier 作业为零**。
- 不允许修改或删除真实 Stage8g/h/i 计算结果；不允许为了给出正结论调整 band label 或挑选 k。
- 如远端 WAVECAR 无法被解析/存在变更，报告一个具体阻断项及自动恢复状态，不让用户手工传文件。
- 报告必须区分 (a) band index 更正是否完全通过，(b) 实际跨隙 pseudo spin vertices 是否非零及数值稳定性，(c) PAW spin/momentum operator 准备情况，(d) 是否值得下一轮独立授权少量投影后处理作业。
- 最终交付：中文一页科学结论 + commit/CI + 正确 JSON + 缺口/下一门槛，不做空泛 status。
