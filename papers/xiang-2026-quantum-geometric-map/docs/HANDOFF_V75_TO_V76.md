# v75 → v76 理论闭合交接 | 2026-10-23

原论文 Xiang et al., PRB 113 L201406 (2026) **整体科学状态 PARTIAL**。基线 v74 PR #88 已在 GitHub 合并且远端独立运行 1105/1105 PASS，详见附带 v74→v75 Handoff。本轮在此基础上继续推进，但不能把局部物理定理误称整篇论文完全复现。

## v75 最重要的物理结果

1. **符号问题成为严格的源公式 no-go：** 采用原文 v49 SHA256-pinned SM Eq.(20)、v_nm=i eps_nm r_nm、Z_nm 的共同约定，必得 +tau sum Z_nm df_n；印刷 Eq.(31) 为相反负号。全局反转 spin/B/position/change of units 只会共同翻转两条路径，消不掉相对负号。随机厄米四带直接收缩 -0.203967264526674，印刷分支 +0.203967264526674。作者原意仍 OPEN；不可单方面宣称勘误。
2. **对微观 tau 给出受限的严格闭合：** 选定 Dirac 导带固定能量 E，等能壳 r(E,phi)=E/(1+beta cos(phi))、DOS 权重 ~(1+beta cos(phi))^-2，标量杂质 spinor overlap 1+cos(delta phi)。构造明确 elastic Born-style 核 W(i→j)=gamma w_j[1+cos delta]. 检查 detailed balance、粒子数守恒、唯一平衡零模、正谱和解析/时域/独立加权谱解。在直到无倾斜时 gamma tau_tr=2 而 gamma tau_q=1；beta=.125 为 1.961508，beta=.5 为 1.492820。不能把选定电场 pole eta_E 不加定义地认作 quantum gamma；必须选 transport collision mode。
3. **反例：** 倾斜 beta=.125、同一半径 r=.3 下不同 phi 能量跨度 2 beta r=.075，所以固定半径均匀 angular scattering **不是弹性过程**。正确模型必须用变形等能壳，并区别 density-of-states weights。
4. **硬边界：** 仅单导带弹性能壳、声明的短程 Born-style 角核，尚未构造与源双频率全矩阵 EB Liouville 一致的碰撞超算符；源原始作者旋向、Eq.(31) 印刷号、Fig2 2π、orbital MNHE 和三带项仍 OPEN。整体 PARTIAL。

新增单模判据：电场源只有满足 Ls=lambda s，标量 tau=1/lambda 才严格成立。beta=0 时模残差约 1.25e-17，beta=.125 时 0.0043664，beta=.5 时 0.0168601；有限倾斜的 tau_eff 是电流投影，不等于所有响应共用的 tau。

## v75 数值与文件

- `src/qgmap/elastic_spin_collision_v75.py`：守恒/细致平衡等能壳弹性散射、稳态与 AC、独立带权谱/时间解、SM 相对符号不变性。
- `tests/test_elastic_spin_collision_v75.py`：51 named physical tests。
- `tests/test_elastic_spin_collision_v75_frozen.py`：1 independent frozen regeneration gate。
- `tests/test_source_v50_collision_bridge_v75.py`：4 optional tests 需要 GitHub 原仓库 v49/v50 源码。独立 v60+ snapshot 中应标记 SKIP，严禁称 56/56 本地 PASS。
- `scripts/reproduce_v75.py` 和 `results/frozen/v75_elastic_collision_and_source_sign.json`；完整详细物理报告 `docs/ELASTIC_COLLISION_TAU_AND_SOURCE_SIGN_V75.md`。

## v76 最优研究任务

优先将 *物理可守恒的能壳碰撞* 实际嵌入原论文双频率 Liouville 的对角占据空间，同时保持真实平衡静态 B 的准备过程、两种场插入次序以及源 E pole；用同一个杂质/碰撞模型给出所有响应通道的一致 \(\tau\)，不要单纯用不同核的任意 gamma 把已有数值凑到 Eq.(20)。必须把 per-shell elastic kernel 扩展到各能量与 valence/conduction、谱/非对角相干投影和配对顶角。对源 Eq20/31 保留双分支验证；作者纠正前 no-go 仍在。随后解决 v52 的绝对 SI 2π 和 orbital Eq49/52、三带几何/磁位置补全。

本轮 GitHub 发布的确切 PR、HEAD、CI、merge 和 postmerge 只有在真实远端验收完成后补记；项目管理完成度不是物理证明比率。


## v75 merged release verification | 2026-10-23

- [PR #90](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/90) **squash merged** as `73d79b602f997ac6ab62c19c5b8e5abe879d5bb5`; final tested PR HEAD `d401721fcea0c995caa3ee58863a19e9ed887f61`.
- [PR final HEAD Actions #39882521001](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39882521001): **1161/1161 PASS**, 177.902 s, source article explicitly `PARTIAL`; all four historical v49/v50 bridge gates individually **executed and passed**.
- [Independent post-merge main Actions #39883014821](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39883014821): **1161/1161 PASS**, 100.137 s on exact merge SHA `73d79b602f997ac6ab62c19c5b8e5abe879d5bb5`; original paper `PARTIAL`.
- v75 has **51 new physical science tests + one independent frozen-number regeneration + four original-source v49/v50 cross-checks**, all **56** executed within complete remote test suite. Standalone v60+ archive includes **52 PASS + 4 intentionally SKIP**, not 56 standalone PASS. Named local v71–v75 module groups separately **225/225 PASS**. All eight new source/test/data/science Git blob SHA matched exactly the tested local files before the documented final HEAD.
- The exact scalar transport relaxation reduction `tau=1/lambda` requires the electric source to be a genuine collision eigenmode `L v=lambda v`. No-tilt Dirac: `gamma*tau_tr=2`; tilted source eigenmode residual nonzero (beta=.125: 0.00436638, beta=.5: 0.01686008). This is an *explicit single-conduction-energy-shell elastic model*, **NOT** original article full conserving spin EB relaxation.
- Eq20/31 common-source relative sign incompatibility is verified but author intent/physical measured sign is **unresolved**; no correction is asserted. Original bulk PHE/MNHE paper science remains **PARTIAL**. These post-merge updates are documentation-only and do not alter any code, tests or frozen numerical results.
