# Stage8p | 本地 DeepSeek harness 自主零新作业：非零三方法门槛与真实 WSe2 局部响应验收

本任务只需用户给 harness 一句“继续执行 GitHub 最新 Stage8p”，**不需要用户手动 SSH、下载 WAVECAR、传输数 GB 波函数、逐条运行命令或复制输出**。所有文件定位、代码运行、整理、修改、测试、推送都由 harness 在已有授权 WSL/日本 IMR 上自主处理。

## 已备好代码与必须先理解的错误

仓库：`Xin-Jiaqi/quantum-geometry-reproductions`；分支 `research/bsf-nmee-stage7-sw-gate-20261020`；不要合并 Draft PR #3。

阅读：
- [Stage8p 独立实证评审](STAGE8P_NONZERO_METHOD_AND_REAL_DATA_REVIEW_20261023.md)
- [非零解析基准](stage8p_nonzero_benchmark.py) 与 `../tests/test_stage8p_nonzero_benchmark.py`；
- [统一源场/投影 Kubo 核心](stage8p_projector_kubo_core.py) 与 `../tests/test_stage8p_projector_kubo_core.py`；
- [直接读取现存 WSe2 96-band EIG/SPN/MMN/WIN 的八项局部三方法核验](stage8p_real_zerojob_three_way.py)；
- 历史 Stage8o-real 报告、`stage8o_diag5.py`、`stage8o_kubo_cross.py`；
- 冻结 Stage8i WAVECAR 和 Stage8m 两次真实 VASP2WANNIER90 后处理的原始 SHA manifest。

必须注意 **Stage8o 所谓三方实现本来就不一致**：projector 路线把 `g=2` 再乘到完整 Pauli SPN_z 上，并把 λ 自旋源写成 σx 而非 σx/2；还把 `dM` 直接当作 `dP/dk`，邻居右端用中心 k 的 U。原始报告没有声称三方通过，是正确的；不得把其网络取回问题当作唯一阻碍。

## 零新作业具体操作

1. 自主核验最新研究分支、文件、安全边界和两套真实后处理目录；对已有 WAVECAR/EIGENVAL/SPN/MMN/WIN/POSCAR 核对 SHA（既有 manifest）。保留旧 Stage8m/n/o JSON，所有新文件使用版本号。
2. 在有 NumPy 的本地环境先执行 `stage8p_nonzero_benchmark.py` 与两组 Stage8p 回归，要求解析非零 \(\Gamma(k=0.4)\approx0.07862712\)、projector/resolvent/Wilson **三个都非零、符号和模量一致**；带基 MMN 形式经任意邻居 band 相位变换后不变。不允许以全零 toy benchmark 代替。
3. 在已有日本 IMR `soc_AB_s0.00_pp8m` 与 `soc_GLIDE_s0.50_pp8m` 后处理产物上运行 `stage8p_real_zerojob_three_way.py`，分别使用新输出路径。此脚本选择**真正** 21×21 网格点 (8,2)、(7,9)，两种 λ 源 x/y、两种 k 轴，0.25 T 双侧 B 中心差；用 **同一** `H=diagE+lambda SPN_source/2 - mu_B B SPN_z` 对比 resolvent 的 analytic dP/dλ、finite-projector 的 dP/dλ 和 Stage8m 旧 Wilson estimator。按源 SHA、gap、min link singular values、MMN link nonunitarity、邻居 projector 幂等误差审计。
4. 重新找回旧 `stage8o_diag5.py` 在日本已 RC0 写出的远端小 JSON，若确实存在只取小型 JSON 证据。登录 stdout 被网关截断不能等同于本地文件不存在。不要用户代为搬运任何东西，也不要将大型许可/原始数据入 Git。
5. 回答物理问题：新两条“resolvent / projector fd”在**非零**真实 WSe2 k 是否一致？它们与 Wilson 是否相同符号/数量级？若它们一致但 Wilson 失败，优先排查两步 k-chain、可数值观测的 MMN 有限窗口非酉、reciprocal boundary 和源-field gauge 的差异。若两条本身不一致，检查源强度/FD 精度/邻居 Hilbert map。**两条共用 MMN 仍不算完全独立数据；没有完整 PAW layer spin、96/97 sewing 和全 BZ 结论。**
6. 实际 glide \(\{m_z|t\}\) 在所有 in-plane k 保留不变，而对 total transverse physical spin odd；不要拿总自旋 SPN 对角或全 Gram median 自动推出材料 layer-odd Gamma 数值。不要因为 band-internal near-degeneracy 就排除一个完整 occupied 36-band projector（占据子空间物理上可以规范协变）；但 96 带窗口是否封闭需要另检。
7. 选择清晰 STOP/GO：`NONZERO_METHOD_PASS_REAL_WILSON_KUBO_FAIL`、`REAL_SOURCE_BASIS_MISMATCH`、`REAL_WILSON_KUBO_LOCAL_AGREEMENT_ONLY`、`BLOCKED_RAW_INPUT`。即使全部局部对照通过，也仍不能称现实 WSe2 NMEE 预测成功，因为 AB C3v/full BZ 与 glide total null 尚失败。

## GitHub 交付

请真正执行，不要仅汇报计划。提交新 sanitized JSON（两堆叠）、源 SHA、局部三方法结果、脚本修订、回归、中文物理报告；更新 `RESEARCH_INDEX.md` 和论文侧导航 `research/README.md`；把关键结论回复 Draft PR #3（保持 Draft/OPEN，不合并），不得修改冻结论文模块。

所有大数据、POTCAR、WAVECAR、MMN/SPN、凭据与私人路径留在合法远端。**授权范围：零新增 PBS、VASP、SOC、NEB、Wannier 计算或新 k-mesh**；仅已有数据的只读 Python 后处理。

最终报告必须明确：**Stage8o 原三方法比较究竟为什么没有可靠通过？修正后真实局部响应是否取得正负号一致的非零独立数值？如果不一致，下一步最值得修复哪一个具体算符或链接，而不是再做大规模 DFT？**
