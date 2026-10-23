# v75 | 原论文符号不可兼容性与守恒碰撞模型的弛豫时间

**2026-10-23。原论文 Xiang et al., PRB 113 L201406 (2026)。整体科学验收仍为 PARTIAL。**

来源分层：原论文印刷公式来自此前基于原始 main PDF 与 16 页 SM 的 v49/v50/v51 审计（主文 SHA256 `b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730`，SM SHA256 `9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4`）。本轮没有宣称重新查看 PDF，也不声称作者承认笔误。弹性散射模型由本轮独立构造，并非原文自动包含的机制。

## A. Eq.(20) → Eq.(31)：明确的代数结论

按原文共同约定，\(\varepsilon_{mn}=\varepsilon_m-\varepsilon_n\)、\(v^a_{nm}=i\varepsilon_{nm}r^a_{nm}\)，以及

\[
\rho^{EB;5}_{mn}=-\tau\,\sigma^c_{mn}\frac{\partial_b f_n-\partial_b f_m}{\varepsilon_{mn}}E^b\bar B^c.
\]

因此对速度直接收缩有

\[
 j_a^{EB;5}/(E^b\bar B^c)
 =-\tau\sum_{m\ne n}\sigma^c_{mn}\,i\frac{\varepsilon_{nm}}{\varepsilon_{mn}}r^a_{nm}(\partial_b f_n-\partial_b f_m)
 =+i\tau\sum_{m\ne n}\sigma^c_{mn}r^a_{nm}(\partial_b f_n-\partial_b f_m).
\]

利用哑指标交换，该式等于 **\(+\tau\sum_{nm}Z^{ac}_{nm}\partial_b f_n\)**，其中 \(Z^{ac}_{nm}=i(r^a_{nm}\sigma^c_{mn}-r^a_{mn}\sigma^c_{nm})\)。印刷的 SM Eq.(31) 却是**负号**。同时反转自旋、磁场、位置算符或全局电荷约定，只会同时改变两条路径，**不会改变相对负号**；若人为定义相反的 \(Z\)，则必须同步重检 Eq.(33) 的自旋导数恒等式，不能只给 Eq.(31) 改符号。带标号互换及带规范 \(U(1)\) 不消除矛盾。本轮随机厄米四带矩阵检查：Eq.(20) 直接收缩 \(-0.203967264526674\)，印刷 Eq.(31) 为 \(+0.203967264526674\)，二者相反且非零。**这定位了原文印刷公式的不可兼容性，但没有证明实际物理磁电系数应选择哪一个最终符号，也不是作者认可的勘误。**

## B. 用明确的弹性碰撞算符，区分量子寿命与输运寿命

选定原文两带 tilted Dirac **导带、单一正能量壳**，并限定 \(|\beta|<1\)：

\[
\epsilon_+(r,\theta)=r(1+\beta\cos\theta),\quad
r_\epsilon(\theta)=\frac{\epsilon}{1+\beta\cos\theta},\quad
w(\theta)\propto\frac{1}{(1+\beta\cos\theta)^2}.
\]

这里 \(w\) 来自二维态密度 \(r\,dr\,d\theta\) 与弹性 \(\delta(\epsilon_k-\epsilon_{k'})\)，并非简单等半径角平均。对质量为零、自旋动量锁定 Dirac 电子，标量杂质的自旋重叠为 \(|\langle u_\theta|u_{\theta'}\rangle|^2=[1+\cos(\theta-\theta')]/2\)。在给定金规则比例系数 \(\gamma>0\) 后定义明确的 **Born-style toy** 跃迁率

\[
 W_{i\to j}=\gamma w_j[1+\chi\cos(\theta_i-\theta_j)],\quad \chi=1\;\text{(Dirac spinor)},\;\chi=0\;\text{(no overlap control)}.
\]

线性碰撞算符 \((Lg)_i=\sum_jW_{i\to j}(g_i-g_j)\) 满足 \(w_iW_{i\to j}=w_jW_{j\to i}\)（细致平衡）、\(L\mathbf1=0\)（恒定占据零模）及 \(w^TL=0\)（粒子数守恒）。\(L\) 在权重内积下半正定。稳态电场驱动求解 \(Lg=v_y-\langle v_y\rangle_w\)，有效的 \(\tau_{\rm tr,eff}\) 从电流投影 \(\langle v_y g\rangle_w/\langle v_y^2\rangle_w\) 定义，而非直接声明某个谱展宽就是 \(1/\tau\)。

**解析特殊情形：** 当 \(\beta=0\) 时 \(w_j=1/N\)，对 \(g\propto\sin\theta\)，有 \(Lg=\gamma(1-\chi/2)g\)。故同一个 \(\gamma\) 在 \(\chi=0\) 下 \(\tau_{\rm tr}=1/\gamma\)，在 Dirac \(\chi=1\) 下 \(\tau_{\rm tr}=2/\gamma\)；此时总散射（quantum）率为 \(\gamma\)，\(\tau_q=1/\gamma\)。这也是 Dirac 体系一般需区分 quantum lifetime 与 momentum relaxation time 的原因。模型与外部散射文献背景一致，例如 Novikov, PRB 76, 245435 (2007) 与 *Spin and Charge Transport in Dirac Materials* (EPFL)。但本轮系数只由这里声明的核给出，**不是**从那些论文借用的数值。

| beta | gamma*tau_tr,eff (chi=1) | 假设 |
|---:|---:|---|
| 0 | 2.000000 | 无倾斜能壳，角向输运模恰为本征模 |
| 0.125 | 1.961508 | 能壳变形，角向本征模混合 |
| 0.5 | 1.492820 | 强倾斜能壳，仅该 toy 内有效 |

特别地，若在 \(\beta=0.125,\epsilon=0.3\) 的同**半径**点间散射，能量范围跨度 \(2\beta\epsilon=0.075\)，所以该捷径根本不是弹性碰撞模型。正确模型将散射限制到变形的**等能**曲线上。另一个重要限定：对于真实 Dirac 标量杂质，\(\gamma\) 仍需从杂质密度、势矩阵元和多带投影中计算；即便角核一致，也不能把其作为从第一性原理得到的实验时间。

## B2. 何时能够严格用单一标量 \(\tau\) 替代碰撞矩阵？

定义电场占据驱动 \(s=v_y-\langle v_y\rangle_w\)，守恒电流子空间中的稳态响应是 \(g=L^+s\)，而不是一般的 \(s/\gamma\)。**严格的单模标量时间**要求 \(Ls=\lambda s\)，且 \(\tau=1/\lambda\)。本轮新增可自动检验的本征模残差

\[
\Delta_{\rm mode}=\frac{\|Ls-\lambda_{\rm R}s\|_w}{\|s\|_w},\quad
\lambda_{\rm R}=\frac{\langle s,Ls\rangle_w}{\langle s,s\rangle_w}.
\]

数值：\(\beta=0\) 时 \(\Delta_{\rm mode}\approx1.25\times10^{-17}\)，模确实闭合、\(\lambda=\gamma/2\)。\(\beta=0.125\) 时 \(\Delta_{\rm mode}=0.00436637917\)、\(\lambda_R=0.03596025876\)，\(\beta=0.5\) 时残差 \(0.01686007827\)。因此有限倾斜下 \(\tau_{\rm tr,eff}=\langle s,L^+s\rangle_w/\langle s,s\rangle_w\) 仅是电流加权投影，不是所有算符共同满足的单一弛豫时间。量化的本征模残差可作为今后完整混合磁电 Liouville 的 **fail-closed** 验收。

## C. 对原论文 \(\tau\) 论断究竟闭合到什么程度

- **VERIFIED：** 原文印刷 Eq.(20)/(31) 的相对符号在共同定义下不可兼容；已有 v50 的 Eq.(20) 极点分支仍是明确的 *formal source route*。本轮额外完成单能壳守恒、细致平衡、实时间响应、带权本征值和独立矩阵求解验证，展示 \(\tau_{\rm tr}\ne\tau_q\) 的可控算例。
- **CONDITIONAL：** 在**声明的 Dirac 单能壳、弹性标量杂质、非相互作用 Born-style 角核**下，可将输运模的 \(\tau\) 解释为碰撞算符的逆本征值；无倾斜情形有解析精确比值。倾斜时需解完整线性碰撞方程，通常不再存在通用标量 \(\tau\)。
- **UNRESOLVED：** 把该单能壳碰撞算符一致地嵌入原文两频率全**矩阵**混合 \(E\bar B\) Liouvillian（含静态磁场准备态、碰撞顶角和自旋依赖）；真实杂质散射数据；作者的原始符号/旋向解释；作者 Fig2(b) 单位约 2π 疑问；轨道 MNHE 全 Wannier 位置磁响应与三带项。因此**没有证明原文的完整 PHE 系数只需取 \(\tau=2/\gamma\)**，更没有完成整篇论文复现。整体 **PARTIAL**。

本轮包含 51 个独立碰撞及符号科学门禁、1 个冻结数据自动重算门禁，另有 4 个需要原仓库 pre-v60 v49/v50 源码的交叉门禁。独立 v60+ ZIP 内 52 PASS、4 SKIP 是真实的分母；完整远端仓库才可执行全部 51 项。

## 下一步的真正理论闭合条件

1. 建立对全两带 Fermi 面、能量积分、真实磁 Zeeman 力矩及 E/B 两种微扰次序一致的粒子数守恒碰撞超算符；检查其 Ward、平衡静态 B 准备、输运 \(\tau\) 与 source SM Eq.(20) 的联系。
2. 对同一碰撞算符，独立解其稳态线性化 Boltzmann 与密度矩阵 Liouville，取可检验的严格极限；比较原文的 Eq.(20) 和印刷 Eq.(31)，不可通过拟合负号强迫重合。
3. 来源勘误仍缺一个外部事实：**作者使用的实际算符旋向及代码，或明确的更正**。若无此信息，严格学术陈述应是“印刷公式不可同时成立”，而不是声称整篇论文错误。
4. orbital MNHE 必须独立引入 Wannier 磁位置修正、轨道矩、跨带三带项和原文物理量纲；这一块目前不能用自旋的碰撞闭合替代。
