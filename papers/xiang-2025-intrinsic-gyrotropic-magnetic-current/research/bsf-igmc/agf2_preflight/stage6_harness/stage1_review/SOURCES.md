# 来源清单（Stage 1 文献与计算输入审查）

> 下载日期：2026-10-20。哈希为 SHA-256。机器可读清单见 `results/SOURCES_MANIFEST.json`（schema `stage1_source_manifest_v2`）。
> **版权与再分发政策**：本仓库**不包含**任何出版社/arXiv 的原文 PDF、TeX 源、原图或全文 SM 提取文本。上述文件仅为**本地私有归档**（`sources/` 目录，`.gitignore` 全部忽略），用于溯源哈希与逐字参数提取；公开仓库中只保留 URL、SHA-256、协议短摘录与本说明（PR #8 审查意见第 1 条）。

## A. 论文与补充材料（Zhu et al., PRL 135, 056801 (2025)）

> 以下文件保存在**仓库外**的本地私有归档目录（会话工作区中仓库根的兄弟目录 `private-agf2-stage1-archives/`），**绝对路径不外泄、不提交**；仓库中只留文件名、URL、SHA-256 与协议短摘录（PR #8 审查意见第 1 条：版权文件不得进入公共暂存区，`sources/` 目录仅保留本政策说明）。

| 私有归档（不入库） | 说明 | 链接 | SHA-256（本地） |
| --- | --- | --- | --- |
| `zhu2025_prl_published.pdf` | PRL 正式版 PDF | https://journals.aps.org/prl/abstract/10.1103/dmzg-ck2t | `9ed58200…` |
| `zhu2025_prl_sm.pdf` | 正式版补充材料 PDF（§1.2 约定、§2 计算方法、Fig. S1–S4、Table S1–S4） | 同上（Supplemental Material 链接） | `37ff3366…` |
| `zhu2025_arxiv_v2.pdf` | arXiv 预印本 v2（含更完整推导文本） | https://arxiv.org/abs/2502.17095 | `59c53499…` |
| `zhu2025_arxiv_v2.tar.gz` | arXiv TeX 源包（**不得整体再分发**） | 同上 | `704aab4f…` |
| `zhu2025_prl_sm.txt` | SM 全文 pdftotext 提取（私有工作副本；仅短摘录取用） | — | manifest |
| `zhu2025_src/` | 解压出的 main_text-clean.tex、Figure_1–3.png、00README.json（原图/TeX，**不入库**） | — | manifest |

完整哈希列表见 `results/SOURCES_MANIFEST.json`（含每项的 `availability=private_local_only`）。

**关键事实**：
- DOI `10.1103/dmzg-ck2t` 真实存在（Crossref 验证：PRL 135, 056801 (2025)，2025-08-01 出版）。
- PRL 数据可用性声明：**"The data that support the findings of this Letter are not publicly available. The data are available from the authors upon reasonable request."** → 作者 B1/B2 原子坐标、INCAR/CHGCAR 等**不公开**（建议按审查意见第 7 条向通信作者索取）。
- SM §2（计算方法，**短摘录**，见 `results/agf2_protocol_frozen.json`）：VASP+PAW、PBE、GGA+U（Dudarev）U_eff(Ag-4d)=4.0 eV [Derzsi et al. PRB 105, L081113 (2022)]、DFT-D3（Grimme 2010）、SOC、ENCUT 500 eV、k 间距 0.02×2π Å⁻¹、真空 20 Å、力 <10⁻² eV/Å、电子 10⁻⁶ eV、Berry 相极化、CI-NEB、Wannier90+WannierTools；辅助：Bilbao、pymatgen、spglib、PyXtal、VASPKIT、VESTA。
- SM §1.2（约定）：层群标准设置遵循 **Spglib**；算符符号遵循 **Bilbao**（m₁₀₀ = 法向沿 [100] 的镜面；m₀₀₁ = 法向沿 [001]）。
- 正文关键数值（用于未来复现比对）：P = ±1.66 pC/m、势垒 17.1 meV/f.u.、Type-I 比 Type-II 低 0.26 meV/f.u.、最大自旋劈裂 51 meV、易轴 [100]、磁点群 m′={1, m′₁₀₀}、双层层群 pb11、单层层群 p2₁/b11（No. 17）、单层连接算符 {m₁₀₀|(1/2,1/2)} 与 {2₁₀₀|(1/2,1/2)}、Type-II 在 4% 双轴拉伸应变下为基态。

## B. C2DB（公开第三方数据，CC BY-SA 4.0）

| 条目 | 说明 | 独立检索途径 / 位置 |
| --- | --- | --- |
| 单层 AgF₂ UID `2AgF2-1` | LG p2₁/b11（No. 17），label exfoliated02-21，PBE/GPAW 弛豫 | 网站：https://cmr.fysik.dtu.dk/c2db/ （按 UID 检索 2AgF2-1）；注册表 `C2DB_all_data.csv` 第 4194 行（id=4194）↔ ASE sqlite system 4195（用户本地库 `project-cy/c2db.db`，2024-09-26） |
| 原始行校验 | 注册表行的 SHA-256 与转换算法记录 | `structures/C2DB_2AgF2-1_monolayer/provenance.json`（含 source_row_uid_sha256、转换脚本参数与输出哈希） |
| 数据许可 | CC BY-SA 4.0（C2DB 数据许可声明；另见用户资料库 `10_材料基础文件/DATA_LICENSES.md` 注册表） | https://creativecommons.org/licenses/by-sa/4.0/ |
| 旧库交叉验证 | 2022-11-30 库 system 10217（同晶胞、E=−19.3594 eV，确认条目稳定） | 用户本地 `project4/c2db-2022-11-30.db` |
| C2DB 论文 | Sødequist & Olsen, arXiv:2401.05992 / *APL* **124**, 182409 (2024) | https://arxiv.org/abs/2401.05992 |

> C2DB 结构为**第三方 PBE/GPAW 弛豫结果**，不是作者的 AgF₂ 结构；仅作为"构造候选"的母层（用户已条件接受其为**待弛豫初始结构**，不视为作者弛豫结构）。

## C. 工具与库（pinned）

| 工具 | 版本 / commit | 来源 | 说明 |
| --- | --- | --- | --- |
| group-theory-operations-toolkit | **v0.18.0**（commit `40a68e9ee2235a765e768f07f0f11fb99729c5c2`） | **https://github.com/Xin-Jiaqi/group-theory-operations-toolkit**（用户仓库；tag v0.18.0） | 对称性/磁层群分析专用（pinned 克隆；verify 脚本经环境变量 `GTOT_TOOLKIT_SRC` 指定其 src 目录，默认指向本机 pinned 克隆） |
| spglib | 2.5.0 | Python 环境 | 3D 空间群/磁性空间群分类 |
| Python | 3.12.4 | 本机 /Library/Frameworks/… | 运行环境 |
| pdftotext | poppler（/opt/homebrew/bin） | — | SM 文本提取（仅工作副本） |

> 注意：本机 `~/Documents/GitHub/group-theory-operations-toolkit` 为 v0.15.0 且工作树脏，**未用于 pinned 分析**；所有分析使用上述 v0.18.0 pinned 克隆。

## D. 许可与脱敏声明

- 所有源文件仅作研究评审用途；SM 短摘录仅用于协议参数逐字冻结（`results/agf2_protocol_frozen.json`）。
- **GitHub 上传脱敏**：不提交 POTCAR、OUTCAR 等 VASP 原始输出与受许可输入；不提交出版社/arXiv PDF、TeX 源、原图、全文提取（`sources/.gitignore` 强制）；只提交评审文档、构造结构文件（公开 C2DB 派生）、哈希与脚本。
- 结构文件哈希（POSCAR/CIF/provenance）见 `results/STRUCTURES_MANIFEST.json`。
