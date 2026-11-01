# 王健组近期文献覆盖核查

核查日期：**2026-11-01**。对象是深圳大学/香港大学的Jian Wang及其合作者，不混入同名作者。窗口以2025–2026为主，回溯2023–2024相关基础工作。

**结论：尚未全部接入。** 本次发现并核实的2025–2026 PRL/PRB集合共13篇（12篇已发表、1篇已接收）：其中4篇的指定核心模型可用，9篇没有论文级验收。另有已接入的Curvature–Flux预印本，因此仓库当前是**5个核心模块**。这些数字是本次检索集合，不是“全组全部论文”的穷尽统计。

## 已接入

| 工作 | 来源 | GitHub可用范围 |
|---|---|---|
| NMEE | [PRL 136, 226302 (2026)](https://doi.org/10.1103/r3x2-kb93) | E×B自旋磁化、指定AC/静态/印刷协议 |
| Spin Transport | [PRL 135, 146303 (2025)](https://doi.org/10.1103/14mp-263q) | 一/二阶spin polarization/current及完整有限频率内核 |
| Quantum geometric map of magnetotransport | [PRB 113, L201406 (2026)](https://doi.org/10.1103/kd12-4ycj) | 指定spin/orbital PHE/MNHE和多带几何 |
| IGMC | [PRL 134, 116301 (2025)](https://doi.org/10.1103/PhysRevLett.134.116301) | Rashba/Dirac传导与位移电流、频率/SI桥 |
| Curvature–Flux | [arXiv:2604.09725v3](https://arxiv.org/abs/2604.09725) | 几何、通量和带间GHE/KME；截至本次仍按预印本登记 |

“已接入”不等于全图或通用材料完成；具体限制见[当前状态](../RELEASE_STATUS.md)。

## 直接相关但未接入的五篇

| 工作与出处 | 缺少的独立能力 | 现有工具为何不够 |
|---|---|---|
| [Quantum Geometric Fluctuation-Dissipation Relation for Nonlinear Transport，PRL 137, 026301 (2026)](https://journals.aps.org/prl/abstract/10.1103/6qcc-wdtx) | 电场线性阶DC噪声与频率积分shift/injection photocurrent关系，Haldane核验 | 平均电流/自旋响应不等于电流关联函数及其频率积分恒等式 |
| [Third-order intrinsic anomalous Hall effect in altermagnets，PRB于2026-10-19接收](https://journals.aps.org/prb/accepted/10.1103/8nkd-ksnd) | E³本征Hall、spin-group/SOC选择规则、Berry曲率四极矩及Lieb/V₂Se₂O模型 | 一/二阶内核与普通磁点群不能自动提供三阶响应和spin-group判定 |
| [Shot Noise of Photocurrents in Magnetic Quantum Materials，PRL 135, 166702 (2025)](https://journals.aps.org/prl/abstract/10.1103/nsky-zmzp) | magnetic shift/injection DSN、PT选择规则及相应模型 | 几何张量可复用，但还没有光致噪声观测量与专门核验 |
| [Intrinsic T-odd spin Hall effect in altermagnets and intrinsic-extrinsic classification，PRB 111, L241404 (2025)](https://journals.aps.org/prb/abstract/10.1103/79nf-1g49) | 该文指定T-odd通道、Drude/本征分解、模型与极限顺序 | Spin/Cui已具备部分spin geometry/vorticity，但尚无该文的独立协议验收 |
| [Unconventional Hall Effect in Gapless Superconductors，PRL 135, 166304 (2025)](https://journals.aps.org/prl/abstract/10.1103/fql8-f3tl) | BdG、准粒子流→横向超流、热力学/四端对照 | 普通Bloch带和正常态NEGF不等于超导电流与粒子-空穴规范一致性 |

交替磁论文旧题为[“as a transport fingerprint of altermagnets”](https://arxiv.org/abs/2604.26665)。本次找到APS接收页及DOI **10.1103/8nkd-ksnd**；不凭数据库“published”字段把接收日期当成正式发表日期。

## 更广义输运线的四篇遗漏

| 工作与出处 | 当前缺口 |
|---|---|
| [Hidden Dissipation in Topological Insulators: Near-field Radiation from Local Shot Noise，PRL 137, 046303 (2026)](https://journals.aps.org/prl/abstract/10.1103/drhy-gkr1) | 有限频率局域shot noise与电磁环境耦合、近场辐射耗散；全局电导/噪声不足以替代 |
| [General method for calculating transport properties of disordered mesoscopic systems based on NEGF，PRB 111, 035409 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.035409) | 该文无序平均方法与指定观测量核验；已有SCBA/随机有限模型不自动等价 |
| [Scaling behavior and emergent conductance plateau … domain wall structure，PRB 111, 035428 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.035428) | 畴壁无序模型、标度与电导分布验收 |
| [Gauge invariant quantum transport theory for non-Hermitian systems，PRB 111, 075135 (2025)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.075135) | 复势系统完整多端电流公式与规范检验；Hermitian Ward工具不等于该文实现 |

以上均经APS作者与深圳大学/HKU单位确认。它们不是本轮自动开工的新任务。

## 相关前序工作，不能因概念相同就算覆盖

- [PRL 130, 036202 (2023)，Quantum Fluctuation of the QGT](https://doi.org/10.1103/PhysRevLett.130.036202)：几何涨落/内禀Hall signatures尚无论文级入口。
- [PRB 107, 075411 (2023)，Third-order intrinsic anomalous Hall effect with generalized semiclassical theory](https://doi.org/10.1103/PhysRevB.107.075411)：完整三阶内核未验收。
- [PRB 109, 115121 (2024)，Linear displacement current solely driven by the quantum metric](https://doi.org/10.1103/PhysRevB.109.115121)：电场驱动的metric位移流不能与磁场驱动IGMC混称。
- [PRB 109, 075419 (2024)，Intrinsic in-plane MNHE in tilted Weyl semimetals](https://doi.org/10.1103/PhysRevB.109.075419)：通用几何部件已有，指定Weyl模型及原文通道未单独接入。
- [PRB 110, 245406 (2024)，Equivalence of semiclassical and response theories for second-order nonlinear AC Hall effects](https://doi.org/10.1103/PhysRevB.110.245406)：现有动态方法可作为基础，尚未建立该文逐通道等价性基准。
- [PRB 109, L201404 (2024)，Gapless superconducting state and mirage gap in altermagnets](https://doi.org/10.1103/PhysRevB.109.L201404)：未接入BdG谱/配对物理。
- [Nature Communications 15, 2012 (2024)，Photocurrent fluctuation by shot noise](https://www.nature.com/articles/s41467-024-46264-1)：未接入噪声内核，非PRL/PRB。
- [Light-induced thermal noise anomaly，arXiv:2412.12662](https://arxiv.org/abs/2412.12662)：未接入热噪声/共振光学协议，本次未确认期刊发表。

## 检索与判定方法

1. 从现有文献清单和Zotero出发，检索APS、arXiv、Crossref；专用academic-search MCP未挂载，采用原站与公开HTTP接口。
2. Crossref按`Longjun Xiang Jian Wang`、`Fuming Xu Jian Wang quantum`检索2023–2026 APS记录（每次100条）；另一个宽查询遇到HTTP429，未作为“无更多论文”的证据。
3. OpenAlex用Jian Wang的ORCID `0000-0002-4519-5003`补查2025–2026，发现新增方向后回APS核对。ORCID聚类仍含3条明显无关的等离子体/医学/工程记录，已排除。OpenAlex漏掉已知新论文且混用接收日期，因此不单独据此宣称穷尽或确定发表日期。
4. 按DOI去重、将arXiv与已发表/接收记录对应，再核查仓库目录、实际API、脚本及测试。文章标题被引用、含有相似函数或Zotero已收藏，均不算接入。
5. 五篇已接入工作的正文+SM/附录及独立物理核验沿用本轮逐篇证据。新发现条目这次完成的是书目/核心观测量与代码覆盖核对，**尚未完成其全文和SM的逐式审计**，不能将目录登记称为新能力。

机器清单：[literature-coverage.json](../physics/literature-coverage.json)。当前判断不依赖测试数量或笼统的PARTIAL标签；后续若选某篇接入，只需为其核心物理建立独立入口与验收，不预先扩成材料平台。
