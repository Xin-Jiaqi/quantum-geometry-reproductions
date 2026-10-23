# Stage8p 独立理论审查：Stage8o-real 的 Kubo 三方比较尚缺共同算符定义（2026-10-23）

> **阶段：零新增 PBS/VASP 作业，独立复核真实 commit `ac487c00253d87b6bc62d180218bb81aba77c975` 的科学报告、净化 JSON、`stage8o_kubo_cross.py`、`stage8o_diag5.py` 与 Stage8m 静态响应代码。** 无权将已经编写、但尚未在日本真实数据上执行的 Stage8p 后处理脚本视作“真实材料 NMEE 已复现”。不覆盖 Stage8o 历史证据，冻结 Jia 等人 PRL 复现。

## 1. Stage8o-real 有效进展

- 真实 P3m1 AB 的 6 个、Aem2 中点的 4 个 Seitz 操作，全 441/441 reciprocal k 精确配对；axial-spin 3×3 Gram 协方差残差中位约 2×10^-7（此前 Stage8n 无效 Cartesian 映射导致的 0.32/1.03 不能作为数据根因证据）。
- Gram 最大残差 AB ～2.4×10^-3，glide ～9.2×10^-3，在小带间隙的 k 点更大；但 **Gram 是单位变换不变量，若 96-band 对称封闭，窗口*内部*的精确/近简并及相位规范本身不应造成 Gram 残差**。需要检验窗口跨越 band96/97、实际计算精度与矩阵表示是否正确，而不能把“简并导致 Gram 大残差”当作独立物理解释。
- Glide 全采样 well-isolated (邻带间距>0.05 eV) 的单态横向 SPN 对角代理最大约 9.7×10^-3；严格 glide 的 exact 非简并真实物理自旋应为零，因此**不是 symmetry null gate 已通过**。哪怕 Gram 中位很小也不能反过来证明所有 operator sewing 精确。
- 当前真实 WSe2 的 AB C3v 响应各向同性与 glide 总量 null 均未通过验收；科学 idea 仍然 OPEN；不应默认重新启动 VASP。

## 2. 关键的新发现：Stage8o 所谓三方对照在进入数据层前已经定义不一致

直接核验以下源代码：
- `stage8m_response.py:build_hamiltonian`：`H=diagE+lambda*SPN_x/2-mu_B*B*SPN_z`；SPN 代表 Pauli sigma，正确 spin-only Zeeman for g≈2。
- `stage8o_diag5.py:H_B`：`H=diagE-mu_B*B*SPN_z`，与上式相同。
- **但** `stage8o_kubo_cross.py:projector_gamma_k` 使用默认 `g=2.0`、`H_B=diagE-g*mu_B*B*SPN_z`，又将 `H(lam)=H_B+lam*SPN_x` 而非 `+lam*SPN_x/2`。因此 projector 与 resolvent/Wilson 路线分别多乘了 **磁源约 2 倍、自旋源约 2 倍**，对 B 导数的比较可产生系统倍率错位。
- 该函数还把反厄米的 Bloch basis connection `dM=<u|∂k u>` **直接当作 Hermitian 的 ∂kP** 放进投影交换子。这在数学上错误。对精确自洽的占据投影导数，应有 `dPk=[D,P]`（在同一正交态基底的导数连接 D 下），或直接用从 *真实 MMN 链* 映回中心基的左右 k 投影差分。
- `stage8o_kubo_cross.py:dM_covariant(U,Mplus,Mminus,dk)` 用同一个中心 `U(k)` 同时变换邻近 `M(k,k±)` 的右端，而实际右端属于**不同 k 的 Bloch band basis**，必须使用 `U(k±)`，或完全避开带本征向量的相位规范，用 center-basis mapped neighbor projectors。`stage8o_diag5.py` 的 resolvent 路线使用 `U(k±)`，但 projector 路线并没有，导致同名“covariant”实际表达式不相同。
- Stage8o `stage8o_diag5.py` 的 `phi_plus/phi_minus` 在乘以 dk·dl 之前先 `abs()`，丢掉 loop phase 方向符号，`|phi_plus−phi_minus|` 不一定等于真实签名相位差。
- 这些问题意味着 **即使 Stage8o 的 gateway 曾取回三方对照 JSON，其三列数值目前也不是同一个可比较的数学量**。Stage8o 对尚未验证事实的保留非常恰当；但修复路径应首先是算符和基底一致性，不是新 VASP 计算。

## 3. 做到了之前没有做到的事情：可解析的**非零响应**三路线基准

新增 `stage8p_nonzero_benchmark.py`，使用明确的两带 Hamiltonian
\[
 H(k,\lambda,B)=\sin k\,\sigma_y+{\lambda\over2}\sigma_x+(m-\beta B)\sigma_z,\quad
 m=1.2,\;\beta=0.4.
\]
这是**方法测试专用的 toy**，不是 TR 对称真实 WSe₂，也不能证明 NMEE 物理新意。其解析零场静态源响应：
\[
F_{\lambda k}(B)=\frac{(m-\beta B)\cos k}
  {4\,[\sin^2 k+(m-\beta B)^2]^{3/2}},
\quad
\Gamma(k)=\partial_B F|_{0}
 =-\frac{\beta\cos k\,(\sin^2 k-2m^2)}
 {4(\sin^2 k+m^2)^{5/2}}.
\]
在 k=0.4，解析 \(\Gamma\approx0.07862712\)，**明显不为零**。解析、occupied projector 双重有限差分、带间 resolvent 和 Wilson 环路必须在同符号、同单位、同源系数下重合，不能“都为零所以测试通过”。
- `stage8p_projector_kubo_core.py` 进一步将上述基准的 Hamiltonian 转成**实际 VASP 需用的 band basis + neighbour MMN links**。对同一中心 k 的 source derivative 使用完整 interband resolvent 与直接 projector finite difference 两套实现；k derivative 使用 `M(k,k±) P(k±) M(k,k±)†` 回投中心 Hilbert space，避免错误把 dM 本身当 ∂kP。
- 合成非零回归测试覆盖正负符号、独立 band phase 改变、必要的 M-link 非酉 leakage 指标、source normalization。**此测试只证明了独立程序对应数学定义，在真实文件上仍需实测。**
- `stage8p_real_zerojob_three_way.py` 已提供 Japan IMR 已有 MMN/SPN/EIG/WIN 的只读本地运行入口：两处实际 21×21 fractional-grid k × 两个 source 轴 × 两条 k 轴，使用同一 \(\lambda\,\sigma_i/2-\mu_B B\,\sigma_z\) 构造 analytic-resolvent vs finite-projector vs existing Stage8m Wilson 数值；返回精确 k 坐标、源 SHA、delta B、非酉链接与投影幂等误差，不声称生产级响应。大文件留在远端，只提交 JSON 摘要。

## 4. 实施策略：先取回最小真实数据，再决定新 DFT

1. DeepSeek 在现有 Japan IMR 授权系统找回先前 `stage8o_diag5.py` 输出 JSON 的远端路径和 size/hash；若文件确实存在，由 harness 自主读取远端小型 JSON，不经过会截断登录网关的长 stdout。
2. 用 Stage8p 新程序分别对 AB 与 glide 的**现存**真实四种文件执行，零新 qsub，得到所有 8 个 local cross-checks，每条给出三路线 \(\Gamma_k\)、error、\(\|M^\dag M-I\|\) 与邻近投影幂等误差。必要时先做一处 k 验证，停止在相位/链接异常，不用稀疏结果直接宣布全 BZ 张量物理可信。
3. 数值意义：即使 resolvent 与 FD-projector 一致，也因为二者**共用同一 MMN 定义的 k 投影导数**，并非完整独立数据来源；再与真正 gauge-covariant Wilson 环路一致才进一步过关。若仍不一致，优先诊断有限 96-band MMN 的**非酉/泄漏**、双步链 vs 中心差导数、dk 太粗、k 边 reciprocal-shift 相位/物理空间算符等。
4. 同时利用已经正确匹配的实际 glide Seitz 操作，按**occupied subspace** 直接计算局部 projector sym constraint，而不是笼统地排除 near-degenerate k：**带内简并只改变单带 gauge，不应改变被全体占据态定义的 projector**。如需 96/97 完整闭合才单独论证新增 NBANDS 计算，不能由小带间距直接认为必须重算全部 SOC。

**门槛：**新的真实材料 \(\Gamma\) 暂未被验收；production quantitative NO-GO，idea OPEN；不新增 SCF/PBS/VASP、NEB/Wannier。Stage8m/8n/8o 旧净化 JSON 和原始科学报告原样保留，并附上指向本报告的更正入口。

**最新科研状态：`NONZERO_SYNTHETIC_RESPONSE_VALIDATION_AVAILABLE__STAGE8O_THREE_WAY_SOURCE_MISMATCH_CORRECTED_IN_NEW_CODE__REAL_THREE_WAY_PENDING__MATERIAL_GAMMA_NOT_ACCEPTED`。**
