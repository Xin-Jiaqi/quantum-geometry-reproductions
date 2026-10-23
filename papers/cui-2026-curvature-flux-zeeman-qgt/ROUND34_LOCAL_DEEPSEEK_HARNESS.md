# Round34 | 本地 DeepSeek harness 一次性交接指令

> **执行对象：本地 DeepSeek harness。** 这是独立任务文档。用户不希望反复复制粘贴指令，你应自行阅读本仓库当前文档、运行脚本、诊断错误并形成一份最终报告。**不要向用户要求逐个粘贴文件内容或冗长的操作结果。**

## 任务目标

解决 Cui 论文项目的真实物理自旋数据来源阻断。尽量取得**同一次、同一 Wannier90 seed 与相同 k 顺序**的非共线 SOC 电子结构、Wannier 变换和 Pauli 自旋源文件，验证：

\[
V(k)=U_{\rm dis}(k) U(k),\quad
S_a^W(k)=V^\dagger(k)S_a^{B}(k)V(k),\quad
H^W(k)=V^\dagger(k)\mathrm{diag}(\epsilon_{nk})V(k).
\]

用 \`hr.dat\` 的独立 Fourier Hamiltonian 交叉核验 \(H^W\)。这是原生文件**数学与部分来源匹配**测试，不是原文材料 SI KME 的完成证明。

GitHub 仓库：\`Xin-Jiaqi/quantum-geometry-reproductions\`。工作目录：
\`papers/cui-2026-curvature-flux-zeeman-qgt\`。
技术细节、已实现代码、官方格式来源：[ROUND34_NATIVE_SPN_SOURCE_CHAIN.md](ROUND34_NATIVE_SPN_SOURCE_CHAIN.md)。

## 安全、权限和成本（硬性）

1. **只读现有数据优先**。不更改用户已有 QE、Wannier90、VASP 或 EPW 计算目录；所有转换、测试和报告放到**单独的本地目录** \`local_evidence/round34/\`。它不应被 git 跟踪。
2. **禁止自动上传任何原始计算、用户文件、伪势、运行日志、凭据或大文件**到 GitHub/第三方。仅生成本地哈希和脱敏摘要。要发送数据必须经过用户明确授权。
3. 不使用付费云计算/集群作业，不执行高开销正式 DFT。若后续需要新算例，优先选现成的小规模**官方非共线 SOC Wannier90/QE 示例**，严格限定约 4 CPU 核、单作业 20 分钟、临时磁盘 2 GB；超出则停止并登记原因。不可用时**不得偷偷安装不明来源软件、伪势或获取额外权限**。
4. 明确区分 \`FOUND_EXISTING_SAME_RUN\`、\`GENERATED_BOUNDED_OFFICIAL_EXAMPLE\`、\`BLOCKED_MISSING_INPUT\`、\`BLOCKED_NO_TOOL\`、\`BLOCKED_GAUGE_MISMATCH\` 和 \`BLOCKED_RESOURCE_CAP\`。不得因文件名/sha256/用户输入的 material_id 直接宣称真实 DFT 自旋已认证。
5. **单次汇报**。能自行解决的路径/格式/测试失败请先自查；除缺少权限、资源超标或无可用源文件这种真实阻断，不要反复打扰用户。

## 自主执行步骤

**A. 获取代码与原始输入候选。** 使用本机已有 Git 仓库，先 \`git fetch\` 并查看 \`main\` 的 Round34 源码。如果没有仓库但网络允许，可克隆 \`https://github.com/Xin-Jiaqi/quantum-geometry-reproductions.git\`。若用户当前有未提交更改，**不要重置、清理或覆盖它们**，在独立只读 worktree/目录运行。

以用户现有的科研项目/计算目录为候选（只在这些有限目录搜索；**不得递归搜索整个家目录或根目录**）。优先找：
- \`seedname.spn.fmt\` 或格式化 ASCII 的 \`seedname.spn\`（注意原生非格式化 \`.spn\` 可能是机器相关 Fortran 二进制）。
- \`seedname_u.mat\`、存在解纠缠时的 \`seedname_u_dis.mat\`、可选 \`seedname_v.mat\`。
- \`seedname.eig\`、\`seedname_hr.dat\`、\`seedname_r.dat\`、\`seedname.win\`、\`seedname.wout\`、\`seedname.nnkp\`。
- 上游 QE \`pw2wannier90.x\` 输入/输出、非共线+SOC 的 SCF/NSCF/投影证据；与原始相同的 Wannier seed、能窗、k 点列表、\`write_spn\`/ \`spn_formatted\`/ \`write_u_matrices\`/ \`write_hr\` 等实际版本设置。

调用已提交脚本（**不是运行网络扫描**），使用你实际找到的有限项目根目录，可以重复 \`--root\`：
~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python scripts/collect_round34_local_sources.py \
  --root "/一个已确认的本地计算项目目录" \
  --out local_evidence/round34/receipt.json
~~~
使用真实本地路径替换示例占位符，不要把上面的中文字符串原样执行。脚本默认限制扫描深度/数量、拒绝递归扫描整个 home、不跟随目录软链、只读文件且不访问网络。可以按照你发现的候选目录多次添加 \`--root\`，一次输出报告。

**B. 现成数据可用时优先核对。** 不要只检查文件存在。使用仓库已提交的 \`src/flux_lab/native_spin_chain.py\` 与 \`src/flux_lab/wannier90_raw_gate.py\`：
- 若 \`.spn\` 是机器相关二进制且本机存在**官方** \`w90spn2spn.x\`，在隔离副本转换为 \`.spn.fmt\`；如果没有官方转换器，只登记 \`BLOCKED_NO_TOOL\`，禁止自行猜测 Fortran 记录布局。
- 原生 \`spn\` 读取 \(\langle\psi_n|\sigma_{x,y,z}|\psi_m\rangle\) 的上三角；原生 \`u/u_dis/v\` 必须按 **Fortran 列优先**读取。检查三维 k 坐标顺序、bands/wannier 维数、isometry、\(V=U_{\rm dis}U\) 和 Pauli/投影自旋的厄米性。
- 只有 \`eig\` 和同 seed 的 \`hr/r\`、真实 \`win/wout\` 等设置完整、\`wsvec\` 不需额外未知补偿时，做独立 H(k) 交叉比较。遇到 unexpanded \`wsvec\` 不要关闭检查，而要记录 \`BLOCKED_GAUGE_MISMATCH\` 与需要的适配器。
- 保存每份文件 **SHA256、大小、来源类别、原始计算软件版本、运行设置摘要**，同时注明何者从原始材料获得、何者只是一份合成/官方示范。禁止将 \`spn\` 的单位角动量误写为 \(\hbar\sigma/2\)；本仓库使用 dimensionless Pauli \(\sigma\) 约定。

**C. 没有现成完整真实数据时。** 若已安装 QE/ Wannier90 和小规模官方 SOC 样例及合法可用的伪势，选一个最小可重复的 spinor 教学样例，在隔离目录生成完整对应链。**只执行资源上限内的单个非昂贵示范**，不要开始 Cui 原文的正式材料 DFT、DFPT、EPW 或付费远程作业。这个阶段结果标为 \`GENERATED_BOUNDED_OFFICIAL_EXAMPLE\`，绝不改写成源论文材料复现。如果样例或工具不存在，直接给出带证据的 \`BLOCKED_MISSING_INPUT\` 或 \`BLOCKED_NO_TOOL\`，不要发起大安装。

**D. 独立验收与输出。** 至少执行：
~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m unittest discover -s tests -p "test_native_spin_chain.py" -v
python -m unittest discover -s tests -p "test_local_round34_collection.py" -v
python scripts/audit_native_spin_chain.py
~~~
最后一个审计使用**明确合成的固定输入**，不能当成新产生的真实物理材料数据。现有真实输入额外用源解析/一致性函数测试，结果和合成审计分开保存。不要因为此时没有本地 \`numpy\` 而暗中改系统环境，先查找本地项目环境；必要安装需要按本地权限规定执行。

完成后在本地保存**且只需要用户转交一个脱敏 JSON 文件**：
\`papers/cui-2026-curvature-flux-zeeman-qgt/local_evidence/round34/receipt.json\`，
另附最多一页 \`SUMMARY.md\`。报告应包含来源状态、文件角色是否齐全、原始哈希摘要、\`spn\` 是否 ASCII、\`Udis U=V\` 残差、HR/eig Hamiltonian 偏差、spin norm/厄米检验、不同算法结果、被跳过或阻断的环节、是否运行新小算例与其 CPU 时间/来源、后续最小缺失文件。

**禁止**将真实源文件、大量日志或本地绝对路径自动写入公共仓库。若报告会公开提交，先生成路径脱敏版并征得用户许可。

## 最终回报格式（只报一次）

写一段简明结论及指向本地 \`receipt.json\` 的路径，给出：
- \`状态: FOUND_EXISTING_SAME_RUN / GENERATED_BOUNDED_OFFICIAL_EXAMPLE / BLOCKED_...\`
- \`物理自旋数学投影: PASS / FAIL / NOT_RUN\`
- \`源文件物理真实性: AUTHENTICATED_BY_INDEPENDENT_SAME_RUN_EVIDENCE / UNVERIFIED\`，并列出具体来源证据（不得凭 hash 自称认证）
- \`真正的 Cui 材料 KME: PARTIAL\`
- \`最小下一步: ...\`

不要让用户来回复制粘贴过程日志。
