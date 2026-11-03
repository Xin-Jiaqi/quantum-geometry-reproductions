# AgF₂ B1_g3.0 15-step PARTIAL：独立科学审查与只读诊断任务

**审查日期：2026-10-22。** 关联 [PR #38](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38) head `7582260897ebc716215c770de3b80b5929d8e8ba`（审查开始时）；旧 PBS Job 490620.cmspbs1 与原 PR #25 均保持独立。

**本文件是对现有 GitHub 计算证据的独立审查，不是新 VASP 运行，也不授权任何新的有偿算力或本地 VASP 任务。** 已核对 `B1_RELAX_RUN_RECEIPT.json`、`B1_RELAX_SCIENTIFIC_ASSESSMENT.md`、`G1_FORCE_AND_GEOMETRY_AUDIT.md`、`INCAR.relax`、`KPOINTS.relax`、`POSCAR.relaxed` 和 main 的原始 B1_g3.0 POSCAR。原始 OUTCAR/vasprun/OSZICAR 仍在研究者私有 WSL，本次审核未亲自访问，需保留此证据层级差别。

## 1. 科学裁决

**状态：PARTIAL / 结构未收敛 / 磁性种子保持 / 电子初步有隙。**

- 15/15 ionic steps 的电子迭代均报告 EDIFF=1e-6 到达，证明电子自洽流程在当前 12×8×1、PBE+U(Ag-4d)=4eV、D3(0)、非 SOC 模型下可执行。
- 初始与末态 TOTEN 分别 -35.17580、-35.34932 eV / Ag4F8 12-atom cell，能量降低 0.17352 eV。仅表明给定计算路径上存在更低能构型，不证明末态局部极小值。Step 14 的能量升高表明 line search 过程需要重查。
- Type-I Ag 磁矩按 POSCAR 顺序约 [-0.566,+0.566,+0.566,-0.566] μB、净磁矩为 0，在这次优化中没有坍塌。**这只是固定磁性分支的自洽解，尚不是磁基态。**
- 初末 occupation-defined gap 约 1.494→1.445 eV，但 k=12×8×1 与高斯展宽 0.05 eV 未做收敛，仅能称“本次采样网格下初步有隙”。
- 离子力最小值是第 6 步的 max|F|=0.051 eV/Å，末步为 0.497 eV/Å，约为本轮 0.02 eV/Å 判据的 **24.9 倍**。无任何已证明的 force-converged geometry。
- Ag–Ag z 分隔由 4.969→4.044 Å，最近跨层 F–F 从 3.681→3.052 Å；第14步短暂到 2.644 Å。**只能得出界面显著压缩且优化器路径曾靠近较短接触的结论**，不能据此称“正在接近局部极小值”“真实平衡界面略短”或“结构已部分收敛至极小值附近”。2.644 Å 高于事先定义的 2.2 Å 最短间距硬停止阈值，但结合力振荡而软停作为保守执行动作有合理性，**不是化学键/反常接触的最终物理判断**。
- 连续降低的总能量和增大的力并不矛盾；梯度大的下坡构型仍可能远离任何稳态。尤其第15步同时具有全轨迹最低能量和 0.497 eV/Å 力，说明不能用能量最低来替代离子稳定性。

### 1.1 独立纠正横向滑移单位

现有材料记录 `t_y=-0.3400 → -0.3391`（分数坐标），那么

\[
|\Delta t_y|=0.0009,\qquad b=5.8796657458\ \mathrm{\AA},\qquad
|\Delta y_\mathrm{relative}|=0.0009\,b\approx\mathbf{0.0053\ \AA}.
\]

旧科学报告写“\(\Delta\approx0.005=0.03\ \AA\)”的 0.005 指标不能用在滑移分数坐标上；它与 0.0009 的记录存在换算不一致。进一步对比仓库初始及第15步 POSCAR，两个 Ag 层**整体 y 方向的绝对平移**分别约 -0.02635 Å 和 -0.03150 Å，差为 **+0.00515 Å**，与 0.0009×b 符合。**0.03 Å 接近各层自身的共同绝对平移，不是层与层的相对滑移。** 正确陈述应为“这条起点及该优化路径中层间横向相对移动约 0.005 Å，远小于 z 向 0.93 Å 的压缩”。

小滑移仅说明**本条计算轨迹**中没有明显横向位移，不能单独证明 B1 具有稳定的滑移势阱、正曲率或有限势垒。确证局部滑移稳定性需要离子充分收敛、对滑移方向做受控位移扰动/能量曲率检查，并避免只沿某一高对称路径停留在鞍点。

### 1.2 重点防止两种过度解释

- **“CG 在 soft interface 模上振荡”是工作假设**，不是已确证的 Hessian 本征软模诊断。存在另一种可能：VASP IBRION=2 的 line-search trial/corrector 步幅不合适、晶格/内部形变耦合、非简谐势面或数值表现。必须先读原始 `trial:`、`trialstep`、`gam`、`g(F)`、离子受力及能量。
- **“B1 堆叠稳定/局部极小”证据不足**，因为在力未收敛时，横向漂移小可能只是这条初始路径受到高对称或力抵消影响；建议保留为“no appreciable lateral drift along this trajectory”。

## 2. DeepSeek 的下一项执行任务：严格只读，无新 VASP

> 直接在当前已有的 WSL 原始 15 步计算目录完成以下分析。可以查找已有 OUTPUT，但**禁止启动、重试或修改任何 VASP/PBS 作业，禁止改 INCAR 或已提交的结果文件**。全程无需用户向 ChatGPT 复制任何中间消息。完成后将新的证据性 Markdown/JSON 直接提交到本 PR #38 的原分支，再发送同一个 PR 链接即可。

1. 找到并校验当前 15 步的原始 OUTCAR、vasprun.xml、OSZICAR；对照这次 receipt 的时间和输入散列。若任何源文件不在，标记 NOT_VERIFIABLE，不用研究者复制粘贴。
2. 从 OUTCAR 提取全部 CG `trial`/`trialstep`/`gam`/`g(F)`/`ort` 以及所有 trial/corrector 标记，结合每一步 TOTEN、energy without entropy（若可）、12 原子的完整 `TOTAL-FORCE` 和离子坐标。**注意 VASP 的试探步与修正步，不要按能量曲线直接推出物理势能有双阱**；给出“哪个步骤开始 line search 过冲、是否因 POTIM 选择或方向更新”的有证据诊断。
3. 将力分解为：两层整体 COM 相对 z 运动、面内相对 y 滑移、各层内部 F/Ag 扭曲。每步计算跨层最近 F–F 笛卡尔最小镜像距离、两个 Ag 层 z-centroid 之差、各层 y-COM；**采用跨 PBC 的最小位移**而不是简单在 wrapped fractional coordinate 求差。报告出现近接触前后的各原子最大受力位置与方向。
4. 审计第 6、13、14、15 步和一次 step 14→15 的能量与受力：哪一项导致 -35.25755→-35.34932 eV 的能量回落？step 15 的 0.497 eV/Å 主要落在哪些原子和坐标分量？若无 Hessian/动态矩阵，不准断言“软模本征矢”。
5. 独立交叉核对当前 `POSCAR.relaxed` 是否逐原子对应最后有效的 OUTCAR `TOTAL-FORCE` 坐标及其 SHA256；保留它作 **unconverged snapshot**，不可当作供发表的结构。原记录说明当时 CONTCAR 未初始化，若真，说明用 OUTCAR 恢复几何是必要的，但需检查最后一个完整的 ionic step 及是否有 predictor-only 坐标。
6. 给出**下一次最小算力计算的参数建议，但不得启动**。只有在看过 CG 的 `trialstep` 和原子力分解后，才能判断是否应尝试较小 POTIM、改变 IBRION 或先做受控层间距离线扫描。不要默认从 force=0.497 的 step 15 自动续跑；可比较 step 6 的较低力与 step 13 的较低能，但也不能擅自挑一个作为“平衡结构”。
7. 将结果推至**同一个 PR #38 分支**，建议文件 `B1_CG_LINESEARCH_FORCE_DECOMPOSITION.md`、`B1_CG_READONLY_AUDIT.json`，并附命令、程序版本、非敏感来源哈希。敏感原始文件继续私有保存。只读分析可标记 PASS/PARTIAL/BLOCKED，不需要新计算授权。

### 2.1 对之后真正继续弛豫的门槛

确认 line-search 问题与 force 分量、确定合适初始快照和稳定执行参数，才讨论**另行授权单条** WSL continuation 或受控 z 向扫描（不会自动转用日本 PBS）。下一条数值计算应该只回答“是否存在一个真正受力收敛的 B1 局部极小值”，仍然不进入 B2、NEB、SOC 或 Wannier。

## 3. 官方程序方法依据（非本项目数值来源）

- VASP Wiki [Structure optimization](https://vasp.at/wiki/Structure_optimization)：IBRION=2 是试探/修正线搜索，可通过 stdout/OUTCAR 的 `trialstep` 评估步长；IBRION=1、2、3 各有适用条件，不能未查日志就机械换算法。
- VASP Wiki [IBRION](https://vasp.at/wiki/IBRION)：RMM-DIIS/CG/damped MD 的定义与适用性。
- VASP Wiki [STOPCAR](https://vasp.at/wiki/STOPCAR)：LSTOP 在下一离子步停止，LABORT 在下一电子步停止；本项目需核对实际 STOPCAR 内容及最后一步完成性。

**结论：** 这次试算真实推进了材料的电子/磁性可行性，但**尚没有提供单个稳定磁性双层结构**。先把已有数据的 line-search 与结构投影审计做透，再决定是否继续占用计算资源。
