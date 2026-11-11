# DeepSeek Harness | PtSe2 BiDB physical-response verification

你是本地 DeepSeek Harness 研究执行代理。自动读取 `Xin-Jiaqi/quantum-geometry-reproductions` 分支 `research/bd36-ptse2-bidb-inputs-20261023` 中 `research/bd36-ptse2-bidb/` 的全部结构与 `manifest.json`，无需让我手工搬运文件。不要触碰仓库 main 和现有研究分支，另建隔离的后续工作分支，提交 Draft PR 而不要 merge。

**目标：** 检验真实 PtSe₂ 双层中，面内 DC 电流 `j_parallel` 是否诱导 `S_-^z=S_top^z-S_bottom^z`，以及滑移依赖性。不要将 BD34 人造 Fourier 向量场当成材料计算结果。

**按顺序执行：**

1. **READ_ONLY 几何 gate**：核对四个 BiDB 原始 ID：2252 单层；2253 是已标记 `Stable` 的中心对称平行双层；2255 同向高对称堆叠；2254 是不同相对取向、非中心对称，不能并入平行滑移族。读取两份从 2253 派生的 `u=(±0.08,±0.19)` 上层刚性滑移 POSCAR，核查实际原子/Seitz 反演和 C3/镜面约束、层距/真空/层标签。把“BiDB 实测或已算数据”与“人为构型”分开标注。
2. **READ_ONLY 环境 gate**：自动检测已有日本 PBS/VASP6/Wannier90 工具链、赝势、SOC spin matrix 接口、可用资源；准备可重跑的 SOC-D3 结构与电子计算输入、收敛表和 PBS job 草稿。**目前没有允许新作业的计算授权。禁止 qsub、VASP 生产任务、远端修改以及 WSL/Slurm 新提交**。要计算，先输出各 job 资源/费用和实验设计并等待我明确批准；不得重复投递。
3. **先小后大**：先 2253 高对称零响应基准 + 2253 刚性通用滑移 (0.08,0.19)，单层 2252 为独立物理对照；2255、反向滑移、2254 保留为后续。清晰指定零场非磁 SOC、固定滑移的有限结构优化、费米能靠近电子/空穴能带而非在禁带内。若构建 Wannier，先验证 Pt 5d + Se 3p 自旋旋量候选：双层 22 个空间轨道，44 个 spinor 基函数；最终投影基组应由投影能带检验决定。必须检验 `hr.dat`、spin `.spn`（或等价有记录的矩阵）、WF 层投影、带插值误差。
4. **响应定义和数值物理门槛**：分别算 `S_-^z` 与 `S_+^z` 的对电场响应和纵向电导率，只有在非奇异 `σ` 的导电状态下才换算电流归一化 `χ^j=χ^E σ^{-1}`。PT 简并用投影子空间迹或其他经验证的 gauge-covariant Kubo，不能任意逐带取 spin。给出 SOC-off、层间耦合、`S_+^z=0`、零点、不同 k 网格与展宽检验。勿将对称性允许误称为算得非零。
5. **交付**：先生成 `PLAN.md`、`STRUCTURE_SYMMETRY_AUDIT.md`、`INPUTS/`、`COST_AND_PERMISSION_GATE.md`、可重复脚本和日志，独立分支 Draft PR；只报告已经执行的事实。取得明确 PBS 提交授权后，再跑有限、幂等的首轮任务并提交数值、原始证据及失败诊断。中文简短结论。

科学基线：原始 BiDB 2253 `gap_pbe=0.219 eV`、2255 `gap_pbe=0.614 eV`，来自数据库，非新计算结果；不能只使用 `hr.dat` 推断原子或层投影 spin。PT 对称性的禁止/允许选择定则要按真实结构检验。
