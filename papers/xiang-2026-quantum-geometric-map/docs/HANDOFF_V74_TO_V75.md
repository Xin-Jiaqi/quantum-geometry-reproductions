# v74 → v75 物理复现交接 | 2026-10-23

原文：Xiang et al., PRB 113, L201406 (2026)，整体科学状态 **PARTIAL**。本轮直接基于仓库 v49/v50 已用原文 SHA256 核验的 SM Eq. (4)–(20)/(31) 约定，**没有重新上传有版权的原始 PDF，也没有声称作者已修改公式**。

## v74 新增且经独立验证的内容

原 v50 已能提取 SM Eq. (15)/(20) 的 \(1/\eta_E\) 极点；v74 没有重复宣称首次实现。新模块使用解析 Fermi Fréchet 导数而非 k 有限差分，精确区分：电场带内输入极点、E→B 路径其余项及完整 B→E 路径；保留两个输入展宽和输出之和。

在 XY 主文旋向、beta=.125,mu=.3,T=.03,nr=64,nphi=88,cutoff=2.5,etaE=.003,etaB=.001 的源自旋 Dirac 模型中，\(\eta_E \chi\) 的 BE/E B 两路分别为 +.005619090327/−.008125646806，合计 −.002506556479；直接 Eq20 收缩 −.002506260111，印刷 Eq31 +.002506260111。已有的源 Eq20↔Eq31 符号矛盾 **仍在**。仅电场带内通道贡献 −.002506387722，附近其余项和 BE 是大项抵消，不可省略。

etaB/etaE=0.1,1,10 三条联合趋零路径以及非零 omegaE 负控制均纳入冻结数据。**不能从离散 4 个调节量点推断严格极限的全部路径交换性。** 静态热平衡 B 的准备方式和真正磁电输运耗散依然是独立未闭合问题。原文使用 selected 1/etaE→tau 的现象学替换，并非对任意所有矩阵通道施加同一个散射率。

## 路线图与验收接口

- **VERIFIED：**两路 Liouvillian 解析磁场导数和电场极点拆分、独立 Eq20 收缩、Fréchet 对完整 rho(k) 数值差分、磁旋向/电子空穴对照、频率与三条展宽路径、原 v50 有限差分的独立交叉（仅在全仓库 CI 实际通过后记录）。
- **CONDITIONAL：**理想 Dirac 两带模型与 cutoff、温度和表观的绝热 pole，仍未定义真实物理碰撞算符。
- **UNRESOLVED：**原文 SM Eq20/31 的相反符号，Main Eq10/SM Eq62 helicity，来源 Eq15 的物理微观 tau，静态 B 准备态与输运判据、完整周期体系 q→0/omega→0 顺序及横向 magnetization current，Fig2 SI 2pi 疑问；独立 orbital MNHE 的全 Wannier 连接、位置场导数、三带和 Eq49/52 争议。

## v75 应优先解决什么

1. 不是继续编写一个同等两带 tau-pole 代码。应在真实原文 index / spin Zeeman / full-current 约定下，给出 **SM Eq20 与 Eq31 的原始代数核查表**，明确哪些可能由符号/磁场方向/电荷约定解释，哪些不可能。若不能闭合，保存两个 literal source branches 并准备精确向作者询问的问题。
2. 真正的磁场静态热平衡准备态、正则化 etaE 与宏观 Drude tau 的关系，要由明确可守恒碰撞模型选出，不能把有限 q toy model 和 v50 Dirac 几何的 pole 数值一致误认为微观散射来源。
3. spin PHE 与 orbital MNHE **分两条线**完成原文量纲、三带和完整 magnetic Wannier 数据。进入材料前严格验证算符、单位和晶格到体相的极限。
4. 用独立作者原始曲线/补充材料文本比较，而非修正系数去拟合 Fig2(b)，保留 v52 曾识别的约 2π SI 归一化疑问。

源码：`src/qgmap/selected_electric_pole_v74.py`。45 新物理测试 `tests/test_selected_electric_pole_v74.py`，冻结值独立重算一项 `tests/test_selected_electric_pole_v74_frozen.py`，额外六项 v49/v50 原始源对照 `tests/test_source_v50_bridge_v74.py`，数值冻结 `results/frozen/v74_selected_electric_pole.json`，复算脚本 `scripts/reproduce_v74.py`，科学报告 `docs/SELECTED_ELECTRIC_POLE_SOURCE_AUDIT_V74.md`。

GitHub PR 与 full-suite CI、merge 后独立 main CI、准确 SHA 等只应在实际完成后写入交接和发布记录；本地 v60+ 历史项目不包含完整 v49/v50 源码，新增 6 项交叉测试在独立交接包中应显式 SKIP，**绝不能写成 52/52 standalone PASS**。

进度估计仍保持 GitHub 工程 90–95%、有限模型能力 80%、原文整体复现 30–40% 的项目管理粗略数值，而非公式或实验复现比例。


## v74 final certified GitHub release | 2026-10-23

- [PR #88](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/88) **squash merged**, final tested PR HEAD `cadbfb47a9c2881f5c88023d5a0742f5f1588aac`, merged scientific source SHA `1e20ae4c937db8ab5f861d7d3fbfff8dd227416f`.
- Final PR-head [Actions 39878472657](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39878472657): **1105/1105 PASS**, `Ran 1105 tests in 158.408s`, `PARTIAL` original-paper status; all six v49/v50 source cross-validation tests explicitly executed PASS, not skipped.
- Independent post-merge exact-sha main [Actions 39878902592](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39878902592): **1105/1105 PASS**, `Ran 1105 tests in 126.641s`, original-paper `PARTIAL`.
- 45 standalone v74 physical science gates plus one frozen independent regeneration **46/46 PASS**; six additional original-v49/v50 comparisons were **SKIPPED** in the independent v60+ standalone package because authentic original v49/v50 modules reside in the complete GitHub repo; these six PASSED in complete remote CI. No local full historical-suite PASS claim is made.
- Final portable ZIP: **124 files**, **123/123 SHA256** matches on clean extraction; v71–v74 standalone subtests **179 discovered, 173 passed + 6 skipped**, frozen JSON separately re-generated with exact equality. The standalone historical snapshot is not asserted to contain all full-repo source-audit modules.
- **VERIFIED** source-Eq15/20 selected electric intraband pole decomposition, differentiated-magnetic Liouvillian and large B->E/E->B cancellation across three finite-regulator trajectories, original v50 source finite-difference cross-oracle. **CONDITIONAL** ideal Dirac source protocol and selected phenomenological (1/\eta_E\to\tau) prescription; universal continuum limit not proven. **UNRESOLVED** original printed Eq20/31 spin sign, main/SM spin helicity, microscopic conserving collision/transport, full orbital magnetic Wannier MNHE and Fig2 SI. Overall original Xiang source remains **PARTIAL**.
- Any after-CI commits in the release record/handoff are documentation-only and do not alter the source/test/frozen values accepted on the exact tested merge commit.
