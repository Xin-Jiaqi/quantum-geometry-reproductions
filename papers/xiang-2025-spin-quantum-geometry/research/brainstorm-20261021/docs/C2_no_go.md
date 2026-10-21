# C2：交流层间自旋转移的严格审计与 NO-GO 预研

**2026-10-21｜基于 PRL 135, 146303 (2025) 正文、SM 及用户提供的 A–C handoff。** 这是独立科学预研，不是已确认新发现，不代表原论文 V32 复现状态发生变化。不修改 GitHub 仓库。

## 0. 结论先行

1. **可行性（PASS）**：一个保留 $m_zT$、破坏简单表示中的 $PT,C_{2z}$ 的四带绝缘双层，确实能对面内 $E_x$ 产生非零**横向自旋极化的层间自旋转移** $\mathcal T_{-,\mathrm{inter}}^x$，且关闭全部层间跃迁后严格为零。
2. **重要区别（PASS）**：$\mathcal T_{-,\mathrm{inter}}^x$ 本质上是层间链路的自旋通量源项，**不是**直接施加在 Néel 矢量上的磁矩力矩。局域交换力矩必须一起算。
3. **连续性 Ward 关系（PASS）**：$\chi_{\mathrm{inter}}+\chi_{\mathrm{local}}=-i\omega\chi_{S_-^x}$（适用条件见下），在 $m_zT$ 和带隙以下给出 $\chi_{\mathrm{inter}}+\chi_{\mathrm{local}}=O(\omega^2)$。**单独的 $\chi_{\mathrm{inter}}(0)$ 可以非零**。
4. **量子几何限制（PASS）**：对层间 hopping 引入相对自旋 twist 参数 $\theta$，可把静态自旋转移电场响应写成常规 **mixed Berry curvature** $\Omega^{k_x\theta}$ 的 BZ 积分。它并不天然等于 Xiang 的 spin QGT，也不能据此宣布一种新几何。
5. **选题结论（NO-GO for C2 as presently phrased）**：已有 Tang–Zhang–Cheng, *Nat. Commun.* **16**, 7790 (2025) 明确讨论**绝缘体中交流电场驱动 Néel SOT、绝热自旋响应和 AFM 共振**。单独“有限频率、层间跃迁依赖、带几何、驱动横向磁性响应”的组合，目前没有建立足够独立的新意。**不要启动 DFT、材料筛选或宣传新效应。** 如要继续，需要证明现有 NSOT 理论没有的受保护谱关系、可切换绝缘体的真实连接操作或明确实验排他识别。

---

## 1. 原论文的正确起点

Xiang–Jin–Wang PRL 的 $J_a^c=\{v_a,\sigma_c\}/4$ 是面内*自旋流*。其 SM Eq. (14)–(17) 将线性 AC 响应分解为 DC intrinsic spin Hall 背景与附加 displacement 响应：

\[
\delta J_a^c(t)=E_b\left[A_\Omega(\omega)\cos\omega t+A_g(\omega)\sin\omega t\right],
\]

其中附加部分的 spin Berry curvature 通道 T-even，spin quantum metric 通道 T-odd，在无吸收且远低于带间隙时分别通常为 $O(\omega^2)$、$O(\omega)$（可被其他对称性消掉）。SM Eq. (59)–(62) 讨论磁点群分类；其印刷 Eq. (62) 右侧存在 $(0;+)$ 与左侧 $(0;-)$ 的指标不一致，推导须从原始反幺正 Kubo 而非直接机械套用该等式。

**与 C2 区别**：本文自定义 $\mathcal T_{-,\mathrm{inter}}^x$ 是面外层间链路的自旋转移，算符不同，频率核也不同，不能未经推导称为原文的 spin-quantum-metric displacement spin Hall。

## 2. 规范明确的 Kubo 与 $m_zT$ 相位表

用 $\hbar=1$，选择场的微扰 $+E_x(t)X$，Fourier 约定 $e^{-i\omega t}$，且全 BZ 投影绝缘填充。对不显含场的厄米算符 $O$，

\[
\chi_{OX}(\omega)=\int_{\mathrm{BZ}}\sum_{v\in\mathrm{occ},c\in\mathrm{unocc}}
\left(\frac{Q_{vc}}{\omega-\Delta_{cv}+i0^+}-\frac{Q_{vc}^*}{\omega+\Delta_{cv}+i0^+}\right),\qquad
Q_{vc}=O_{vc}X_{cv},\ \Delta_{cv}=\epsilon_c-\epsilon_v>0.
\]

在严格 sub-gap（$\omega<\Delta_{\rm min}$）下：

\[
\chi_{OX}(\omega)=\int_{\mathrm{BZ}}\sum_{vc}
\frac{2\Delta_{cv}\operatorname{Re}Q_{vc}+2i\omega\operatorname{Im}Q_{vc}}{\omega^2-\Delta_{cv}^2}.
\]

取两个物理层 $\tau_z=\pm1$，真实自旋 $\sigma$，面外 AFM 序 $M\tau_z\sigma_z$，定义

\[
\mathcal A=m_zT=\tau_x\sigma_x K,\qquad \mathcal A H(k_x,k_y)\mathcal A^{-1}=H(-k_x,-k_y),\quad \mathcal A^2=+1.
\]

$E_x$ 是 $\mathcal A$-even。若 $O$ 在 $\mathcal A$ 下为 $\eta_O=\pm1$，则 $Q(k)$ 与 $\eta_O Q^*(-k)$ 联系，所以积分后的 $\eta_O=+$ 通道只留实数（同相），$\eta_O=-$ 通道只留纯虚数（正交相）。

| 对象 | $m_zT$ | sub-gap AC 允许相位 | 说明 |
|---|---:|---|---|
| 总 $J_y^x=\{v_y,\sigma_x\}/4$ | 奇 | $\sin\omega t$ | 对应原论文 spin metric 的可能存活通道 |
| 层差 $J_{-,y}^x$ | 偶 | $\cos\omega t$ | 对应 spin BC 的可能存活通道 |
| $S_-^x=\tau_z\sigma_x/2$ | 奇 | $\sin\omega t$ | 横向 staggered spin density |
| $\mathcal T_{-,\mathrm{inter}}^x$ | 偶 | $\cos\omega t$ | 层间链路转移 |
| $S_-^z=\tau_z\sigma_z/2$ | 偶 | $\cos\omega t$ | 平行 AFM 易轴，不是所求横向力矩 |
| $\mathcal T_{-,\mathrm{inter}}^z$ | 奇 | $\sin\omega t$ | 平行自旋转移 |

表中的“允许”并不保证非零；若增加其他磁点群操作，仍需重新检查。

## 3. 微观链路自旋通量与连续性方程

取 $S_-^a=\tau_z\sigma_a/2$，对 spin-independent $H_\perp=t\tau_x$：

\[
\boxed{\mathcal T_{-,\mathrm{inter}}^a=i[H_\perp,S_-^a]=t\tau_y\sigma_a.}
\]

一般 spin-dependent 双层 hopping：

\[
H_\perp=\begin{pmatrix}0&\mathsf T\\\mathsf T^\dagger&0\end{pmatrix},\quad
\mathcal T_{-,\mathrm{inter}}^a=
\begin{pmatrix}0&-i\{\mathsf T,s_a\}\\i\{\mathsf T^\dagger,s_a\}&0\end{pmatrix},\qquad s_a=\sigma_a/2.
\]

这显示层间源项是层间跃迁与自旋的*反对易*组合。定义 $S_-=S_T-S_B$ 时，$\mathcal T_{-,\mathrm{inter}}$ 是层差密度因层间转移的变化率；只有指定方向、界面自旋吸收、局域磁矩交换与边界，才能把它转换为真正的磁性 torque。

对宏观均匀体系：

\[
\dot S_-^x=\mathcal T_{-,\mathrm{inter}}^x+\mathcal T_{-,\mathrm{local}}^x,
\quad
\boxed{\chi_{\rm inter}(\omega)+\chi_{\rm local}(\omega)=-i\omega\chi_{S_-^x}(\omega)}.
\]

当 $[X,S_-^x]=0$，且算符本身没有电场依赖时无附加 equal-time/contact 项。在 $m_zT$ 下 $\chi_{S_-^x}=iC_1\omega+O(\omega^3)$，因此两力矩响应之和 $=C_1\omega^2+O(\omega^4)$。**它并不要求 $\chi_{\rm inter}$ 本身在零频消失。**

## 4. 独立可复验的四带二维周期绝缘模型

\[
\boxed{
H(\mathbf k)=\underbrace{t\tau_x+u\cos k_x\tau_y+\lambda\tau_y\sigma_x}_{H_\perp}
+\underbrace{\tau_z(d\sin k_x+M\sigma_z+\beta\sin k_y\sigma_x)}_{H_{\mathrm{local}}}
}
\]

固定轨道规范，所有跃迁项都是周期函数；没有作随 $k$ 变化的伪 gauge 变换。参数无量纲：$t=1,u=0.21,\lambda=0.33,d=0.49,M=1.3,\beta=0.24$，半填充（2 occupied, 2 empty）。其中 $M\tau_z\sigma_z$ 是上下层反平行的交换场；$\lambda\tau_y\sigma_x$ 是层间自旋依赖 hopping。

直接矩阵计算：

\[
\mathcal A H(\mathbf k)^*\mathcal A^\dagger=H(-\mathbf k)\quad\text{(max error }0\text{)}.
\]

选常见的 $PT=\tau_x(i\sigma_y K)$ 与 $C_{2z}=i\sigma_z$ 表示，在非特殊 $k$ 点，这两种操作的矩阵等式偏差均为 $0.66$（归一化 Hamiltonian），即没有这两种额外保护操作。

对 $S_-^x$：

\[
\begin{aligned}
\mathcal T_{-,\mathrm{inter}}^x
&=t\tau_y\sigma_x-u\cos k_x\tau_x\sigma_x-\lambda\tau_x,\\
\mathcal T_{-,\mathrm{local}}^x&=-M\sigma_y.
\end{aligned}
\]

采用 $E_x$ length gauge、$\omega=0.2$，39×39 BZ k 网格，$\omega<\Delta_{\rm min}\simeq2.4524$：

| 计算量 | Kubo 结果（无量纲） |
|---|---:|
| $\chi_{\mathrm{inter}}(0)$ | $+0.00257758472$ |
| $\chi_{\mathrm{local}}(0)$ | $-0.00257758472$ |
| $\chi_{\mathrm{inter}}(0.2)$ | $+0.00259065168$（虚部约数值零） |
| $\chi_{\mathrm{local}}(0.2)$ | $-0.00257773600$（虚部约数值零） |
| $\chi_{S_-^x}(0.2)$ | $+i\,0.0000645784126$（实部约数值零） |
| $\chi_{\mathrm{inter}}+\chi_{\mathrm{local}}$ | $+0.0000129156825$ |
| $-i\omega\chi_{S_-^x}$ | $+0.0000129156825$ |

**交叉校验**：同一个模型分别按本征能带的 $O_{vc}v_{cv}/\Delta$ 计算，以及在固定轨道基对占据投影 $P(k)$ 数值微分、求解频域 Liouville 并回到算符平均值。25×25 网格、差分 $10^{-5}$ 时，两个 oracle 最大误差约为 $10^{-13}$。场符号是 $+EX$；若改用电子 $-eEX$，响应整体符号跟随变化，不影响选择定则和 Ward identity。

**解耦极限**：把**全部** $H_\perp$ 乘以 $q$，在 $q\to0$ 时四带系统仍有能隙，且 $\mathcal T_{-,\mathrm{inter}}\equiv0$。$q=0$ 的数值 $\chi_{\rm inter}=0$。$q=1$ 时 $\chi_{\rm inter}(0.2)=0.00259065$。这证明不是两个解耦层算符相加产生的转移，但只是定义上的必然条件，并非新效应证明。

## 5. 量子几何的严格对应，以及为什么它降低 C2 的原创性

只对真实层间 hopping 作 relative spin twist：

\[
H(\mathbf k,\theta)=H_{\rm local}(\mathbf k)
+U(\theta)H_\perp(\mathbf k)U^\dagger(\theta),\quad
U(\theta)=e^{-i\theta S_-^x}.
\]

因此 $\partial_\theta H\vert_0=\mathcal T_{-,\mathrm{inter}}^x$。此处 $\theta$ 是施加在*真实层间 hopping* 上的源参数，**不是**单纯把整个 Hamiltonian 作一次基底旋转。

对占据投影 $P(\mathbf k,\theta)$ 定义通常的 mixed Berry curvature：

\[
\Omega_{k_x\theta}(\mathbf k)=i\operatorname{Tr}\left(P[\partial_{k_x}P,\partial_\theta P]\right).
\]

在本文 $+E_xX$、$\hbar=1$ 的约定下，静态响应满足

\[
\boxed{\chi_{\mathcal T_{-,\mathrm{inter}}^x,E_x}(0)
=-\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}\Omega_{k_x\theta}(\mathbf k).}
\]

数值独立验证（21×21 BZ、对 $k,\theta$ 的投影中心差分各为 $10^{-4}$）：

- $\int_{\rm BZ}\Omega_{k_x\theta}\simeq -0.0025775847053$；
- $\chi_{\mathrm{inter}}(0)\simeq+0.0025775847159$；
- 两者符号反向、数值差约 $1.1\times10^{-11}$。

**物理结论**：这个响应确有规范不变的几何表达，但它是关于普通动量 $k_x$ 与层间自旋 twist 参数 $\theta$ 的**常规混合 Berry 曲率**。单凭这个公式并不能区分已有几何磁电极化、绝热 torque 和 generalized-force Kubo，也不是 Xiang PRL 特定的 $g^{ab;c}_{nm}=\mathrm{Re}[(J_a^c)_{nm}v^b_{mn}]/\Delta^2$。

## 6. 文献原创性检查

- **Xiang–Jin–Wang, PRL 135, 146303 (2025)**：已建立 spin QG 与线性 AC displacement spin Hall，但没有直接给出此处的 layer-link torque 作为新可观测量。不能把我们的链路变量说成原文已证明的 spin metric torque。
- **Tang–Zhang–Cheng, Nature Communications 16, 7790 (2025)**：在补偿 AFM 绝缘体系中分析 *staggered Edelstein / intrinsic Néel spin-orbit torque*，**明确用交流电场驱动 AFM resonance**。它与我们试图建立的“绝缘体、场驱动、自旋力矩、无 Ohmic current”直接重叠，是决定性近邻。DOI 10.1038/s41467-025-63171-1。
- **Železný et al., Nature Physics 14, 220 (2018)**；**Spin-transfer and fieldlike torques in antiferromagnets, PRB 104, 174402 (2021)**：AFM 中 torque 与自旋传输已有成熟传统。
- **Fujimoto–Matsuo, PRB 100, 220402(R) (2019)**：AC interfacial STT、交换相关共振已有先例。

当前检查不是穷尽性综述，但已足以否决“首次发现绝缘体 AC interlayer spin transfer / AFM torque”的笼统定位。

## 7. 建议：C2 现阶段不进入真实材料计算

**已验证**：$m_zT$ 相位过滤、四带绝缘模型存在性、层间解耦零基线、两算法 Kubo/Liouville、交换与链路力矩 Ward identity、混合几何等价。

**未证明**：该 link transfer 会对磁矩产生净可测扭矩、实验信号独立于已有 NSOT、存在 spin QGT 专属的新 spectral identity 或 quantization、适合真实可切换磁性双层。

如仍保留 C2，唯一值得投入的下一道门槛是：能否结合**真实磁性堆叠的 connection operator**构造一个在现有 AFM NSOT 下严格不存在的、可被相位读出并由 Ward identity 或边界自旋注入独立鉴别的量。如果不能，应**停止 C2**，优先回到 BSF + 群论筛选能提出什么定理，或重审 B 中真正受层间 hopping 控制的物理响应。

### 可复验脚本

`c2_minimal.py`：仅需 numpy，运行命令 `python c2_minimal.py`。内含模型、$m_zT/PT/C_{2z}$ 检查、Kubo 谱求和、有限差分占据投影 Liouville、关闭层间 hopping、spin-twist mixed curvature 验算。执行时间取决于 k 网格，不修改任何仓库文件。