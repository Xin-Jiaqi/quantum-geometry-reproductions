# 公式与数值审计

依据：用户提供的 17 页 PRB 正式版，DOI 10.1103/xr5v-gkv5；主要定位 p3 Eqs. (10),(11),(15)–(17)，p9 Eq. (A9)，p12 Eq. (B17)，p14 Eqs. (D9)–(D11)。公开 arXiv v1 是不同排版版本，不替换正式版的系数。

## 原文表达的实现边界

`printed_D10` 的轨道部分完整保留：

$$
\widetilde G^{B;az}_{n,L}=-\frac12\epsilon_{zbc}\partial_b g_n^{ca}
-\operatorname{Re}\sum_{m\ne n}\frac{M^z_{mn}R^a_{nm}}{\epsilon_m-\epsilon_n}
+\operatorname{Re}\sum_{m\ne p,\,p\ne n}\frac{\mathcal M^z_{mnp}R^a_{pm}}{\epsilon_m-\epsilon_p}.
$$

这里 e=hbar=1 只用于无量纲 q、b 的代数；实际 SI 换算单独进行。M_mn=sum_p M_mpn，M_mnp=-(v_mn+delta_mn v_p)×R_np/2。最后一项的条件没有 m≠n，因此不会因为只有两个能量块就自动删除。

metric 项在周期 BZ 上分部积分为 +int[(partial_x f)g_ya-(partial_y f)g_xa]/2。带隙内零温时该积分为零，其余两项仍可非零。CuMnAs 块求和是本实现的持续简并子空间延拓，原文没有提供完整非阿贝尔算法；测试验证其 U(2) 不变性，并拒绝会局部分裂的简并。

`density_A9` 在固定轨道基下独立采用原文 A9 的 +omega 驱动：

$$
D_L=\{H_y,(\rho_0)_x\}/4-\{H_x,(\rho_0)_y\}/4,
\qquad D_S=i[M_S,\rho_0]/2,
$$

$$
(\rho_B)_{nm}=\frac{(D_L+D_S)_{nm}}{i(\epsilon_n-\epsilon_m)},\qquad
G_a=\int\operatorname{Tr}(\rho_B R_a).
$$

只在不同能量块间求逆，不把精确双重简并的小分母人为加展宽。rho0=f(H) 的导数使用矩阵函数的差商；等能量块使用 f'。独立参考算法直接有限差分 H、f(H)、谱投影，使用 R_a=i sum_s (partial_a P_s)P_s，再求解同一驱动，避免对本征矢相位求导。

## OPEN_SOURCE_DISCREPANCY：轨道幅度不等价

按当前实现，CuMnAs 带隙内得到：D10 的 moment 与 multiband 两项相等，轨道总和为 A9 的两倍；两条路线的自旋项一致。该现象同时通过随机 U(2) 旋转和独立有限差分检查。对通用非简并三带模型，D10 与 A9 仍不同，且无法用同一个 1/2 因子校正所有项。

因此，本目录分别验证两个明确的计算对象：印刷公式的忠实求值，以及其密度矩阵驱动的独立求解。**不指定哪一个已经是唯一正确的实验幅度，不静默改写原文，不以原图幅度反向调整系数。** 后续需要核对 A9/A12 到 B15/B17 的重排或作者代码；本轮核心机制接入不依赖先解决整篇推导。

这一差异也不属于整体 cosine/phasor 因子：整体约定变化会同时改变轨道、自旋，而这里自旋一致。`test_12_source_discrepancy_is_not_hidden` 对应的测试保存差异，不将其设为被忽略的失败测试。通用三带 D10 另外有逐指标循环对照，防止向量化指标错误被当成原文差异。

## 数值证据与未决条件

- phi=pi/12、mu=0、25 K：独立有限差分与 A9 的最大核差约 1.7e-10（41²、步长 2e-5）。
- 161² 的 Néel 反号残差约 1.4e-17；161²→241² 的带隙轨道 x 核变化约 9.4e-11。
- 角度图在 mu=0.08 eV 使用 161²、17 个角度，属于机制图。phi=0 尖峰 161²→401² 轨道核变化约 0.73%，不是全角度高精度冻结基准。
- 示例 a=3.82 Å 为明确输入，上传正文没有给出面内 a。轨道 SI 值正比于 a，自旋值反比于 a，不能无说明地比较幅度。
- 正文/图注对“约两倍”与“接近一个数量级”的描述不一致；未将任何一种措辞作为自动验收阈值。
- finite-T 金属能正常计算，但本次不验收整个低温费米面扫描；零温金属直接拒绝，避免漏掉 delta 型占据导数。

## 相邻模块交叉检验

旧 NMEE 的 CuMnAs Hamiltonian 在显式设置 Jn=0.85 eV、a_nm=1（使旧 dH/dk 的数值等于此处 dH/dq）后必须逐矩阵相同。旧自旋 IGMC `pair_kernels` 提供 Q01=Re(R01 sigma10)；当本测试固定 M=-sigma/2 时，此处 G_S=-旧 sigma_D_per_omega/4。这个代数关系在测试中明确记录，不能将不同论文的相位、磁矩和谐波约定直接混用。
