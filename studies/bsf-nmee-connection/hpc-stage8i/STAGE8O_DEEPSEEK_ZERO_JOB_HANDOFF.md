# Stage8o-real | DeepSeek harness 独立零新作业验收任务（单文件入口）

**目标：** 在**已存在**的真实 WSe₂ VASP2WANNIER90 EIG/SPN/MMN 与原有结构上，检验 Stage8n 宣称的 VASP ISYM=0 “数据层根因”究竟成立与否，进而拯救或否定当前 coarse static Gamma 数据。**这是科学实证诊断，不是重复写计划。** 用户无须上传原始文件、复制多条终端命令或手动 git 归档。

## GitHub 来源
仓库：`Xin-Jiaqi/quantum-geometry-reproductions`，分支 `research/bsf-nmee-stage7-sw-gate-20261020`。
请先 fetch 最新 HEAD（含已新增的 Stage8o 修正），阅读：
- `hpc-stage8i/STAGE8O_INDEPENDENT_METHOD_AUDIT_20261023.md`：独立审查，逐点指出 Stage8n 统计与物理定义上的错误；
- `hpc-stage8i/stage8o_symmetry_math_audit.py`：六角晶格/轴自旋/相位数值检查；
- `hpc-stage8i/stage8o_actual_spn_symmetry_screen.py`：真实全原子 spglib-Seitz 操作 + 全 441 k + full 3×3 spin-Gram 检验；
- `hpc-stage8i/STAGE8N_DIAGNOSIS_REVIEW_20261023.md` 与 `stage8n_diag1.py / stage8n_diag4.py`：原始待更正主张，保留历史不覆盖；
- `hpc-stage8i/stage8m_plaq.py` / `stage8m_response.py`：真正 Wilson 返回值/单位/坐标索引，`stage8n_diag3.py` 及 Kubo 代码。

## 一、自主执行（**零新增 VASP/PBS/QE/NEB/Wannier 任务**）

1. 在已有 WSL/日本 IMR 授权环境 fetch 最新分支；通过 Stage8m 已归档的 manifest/resume-state 自行找出 AB 和 glide 的真实 POSCAR/WIN/EIG/SPN/MMN（既有 492698/492699 作业目录），核对所有原始 SHA。禁止让用户传 4.9GB WAVECAR；禁止任何私人文件/许可证/SSH 凭据公开。
2. 自主在日本有现成 spglib 的运行环境中运行**已提交** `stage8o_actual_spn_symmetry_screen.py`，给出两个结构的输出 JSON 和实际原始输入 SHA。先核实真实 POSCAR 全部六个原子、spglib 的 fractional rotation + translation、晶格设置（实际 z 坐标/原点可能有 fractional 非整平移）。
3. 对每个真实 Seitz 操作必须列出：**441/441** 精确 reciprocal k 对；完整 band energies 在这些 k 对下的最大 sorted-eigenvalue residual；正确轴自旋 \(Q=\det(R)R\) 旋转后的 **3×3 Gram covariance** 相对残差。**切勿使用无效 Cartesian [ [0,-1],[-1,0] ] 映射，不能忽略 reciprocal lattice folding，不能比较未经旋转的单一 Sx²。**
4. 对真正的 horizontal glide \(R=\mathrm{diag}(1,1,-1)\) 注意它对 in-plane k 是恒等映射。检验全网格（除去所有近简并）的 spinor band diagonal \(|S_x|,|S_y|\) 是否接近零，明确分布、极值、对应 k/band/gap；同时在已知 Stage8k 旧采样 near-M/Γ–K 的点做交叉比对。只有全部严谨检查完成才能讨论 VASP 数据是否明显缺失该对称性。
5. 对 96-band window 的未封闭风险单独标记：原文件无 band97；Gram residual 可能是截断而非 VASP 软件对称性失效。若 EIG band spectrum 的晶体对称残差已大，则可独立讨论数值求解/Hamiltonian；若 EIG 通过而 Gram 失败，要继续分开“缺带、SPN convention 与真正 sewing”。
6. 纠正 `stage8n_diag4.py` 的 phase precision budget：直接从同一四 Wilson links 求**角度** \(\phi=\arg\det L_1L_2L_3L_4\)；导出 raw \(\Delta\phi\)（未除以 \(\Delta k\Delta\lambda\)），量化 SVD 极分解/36×36 determinant/网格平均真实误差，而不是直接比较已经除以面积的 F 值和 machine epsilon。
7. 修正 static Gamma 的**量纲**以及 \(k\) 轴斜基到笛卡尔转换，不再沿用 `1/(eV T Å²)` 标签。数值必须说明 spin operator/source convention 与没有 charge prefactors 的限制。
8. 不允许用未经检验的 Kubo 跨方法一致性作“引擎已排除”证明。原先只有 1/10 个局部 k-axis 对照 rel≈.04，其他最大 275×。先在合成**非零**且远离简并的模型实现独立 projector/解析 Kubo 对照，严格检验 spin-matrix conjugation、\(\partial_k U\)、同一 MMN overlap gauge、source derivative 的符号。然后若可行，在真实 96-band 局部 k 上测试。
9. 必须提供两个实际归一化的 STOP/GO 值：**AB C3v tensor residual**（对真实 Cartesian M=γI），**glide total tensor null residual**；只有对称性通过且 Kubo/Wilson 非零响应在误差内一致，才允许重新评价 WSe2 的量级和物理新意。
10. 若原脚本的依赖、POSCAR/spglib 格式或 k-grid 与原数据不符，优先定位并修复**新增 Stage8o 程序**本身，使用 synthetic tests 独立检查；不能为了“看起来完美”随意放大 symmetry tolerance。保留程序/真实输出 SHA 和改动原因。

## 二、交付与科学判定

自动提交安全的新 Stage8o-real evidence JSON、独立脚本修复、非零 benchmark 测试、中文审查报告；更新 `studies/bsf-nmee-connection/RESEARCH_INDEX.md` 和论文侧**仅导航** `research/README.md`，保持 PR #3 Draft，原 Jia PRL 复现冻结。

**只能在证据支持下选择状态：**
- `SYMMETRY_RECOVERED__RESPONSE_STILL_UNVALIDATED`：修正 k/spin 变换后真实 SPN 检查符合，但原 Wilson Gamma 数值尚无独立 Kubo；
- `EIG_SYMMETRY_FAIL`：实际能谱在真实晶体操作下不匹配，查结构、k-point、SOC 运行与真对称性；
- `SPN_SYMMETRY_FAIL__96BAND_WINDOW_UNCLOSED`：能谱过关，完整 Gram 仍不合，原因可能为带截断/PAW/SPN/读取器，不直接指定 ISYM=0；
- `GAUGE_KUBO_INCONSISTENT`：对称性过关但独立非零静态响应不一致；
- `COARSE_RESPONSE_METHOD_VALIDATED__LAYER_ODD_NOT_COMPUTED`：只有原始 SPN/能谱对称性、AB C3v+glide null、非零 Wilson/Kubo、单位和步长/窗口等全部核验通过才可使用；
- `BLOCKED_INPUT_ACCESS_OR_FORMAT`：真实文件或格式无法自动获取，有明确错误日志与恢复状态，不编造 PASS。

**硬边界：**零新增 PBS 作业，零全新 SOC SCF、零 NEB、零生产级 Wannier、零原始 MMN/SPN/WAVECAR/POTCAR/CHGCAR 上传。所有能自主完成的工作自行完成，不把用户当手工命令/传输中间人。

最终中文报告需要明确：**Stage8n 的“VASP 对称性破坏根因”是否真的成立？若否，哪一个环节首先出错？此 WSe2 idea 是否值得下一轮单独授权的材料计算？**
