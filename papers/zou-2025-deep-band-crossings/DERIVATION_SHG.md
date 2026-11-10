# Zou et al. (2025): DBC 对 SHG 的领先发散为何抵消

**来源**：N. Zou *et al.*, *Physical Review Letters* **135**, 226901 (2025), DOI: [10.1103/sqwy-xhdf](https://doi.org/10.1103/sqwy-xhdf)。此处逐项核对正文 Eq. (1)–(4) 及 Supplemental Material (SM) Eq. (S6)、(S7)、(S10)–(S14)、(S19)–(S26)。代码为独立复现，未复制作者程序。

**结论**：当绝缘体内部两条导带 $c,c'$ 逼近简并时，单独的 interband SHG 项 $\chi_e$ 与 mixed interband–intraband 项 $\chi_i$ 可以各自包含 $1/|\mathbf k-\mathbf k_0|$ 型发散。对相同的物理张量、相同的占据和频率约定，SM (S19)–(S26) 保证所有这些**领先发散项**在 $\chi_e+\chi_i$ 中相消。此结论不排除有限但较大的 SHG 响应；它也不代表 shift current 的虚跃迁发散会消失。

## 1. 三带空间、符号与发散的来源

考虑 $v$ 已占据，$c,c'$ 未占据，导带简并差 $\delta=E_{c'}-E_c\to0$，光学间隙 $\Delta=E_c-E_v>0$ 保持有限。使用 SM 的约定：

$$
\omega_{nm}=(E_n-E_m)/\hbar,\qquad
r_{nm}^{a}=-i\,v_{nm}^{a}/\omega_{nm},\qquad
f_{nm}=f_n-f_m.
$$

$D^a r_{nm}^b=\partial_{k_a}r_{nm}^b-i(A_n^a-A_m^a)r_{nm}^b$ 是 U(1) 协变导数。SM 定义的 Hermitian connection 为

$$
C_{mn}^{abc}=r_{mn}^{b}D^a r_{nm}^{c}.
\tag{1}
$$

**索引关键**：$C_{mn}$ 的前因子是 $r_{mn}$，但导数作用于 $r_{nm}$。实现时若将两者同时写成 $nm$，会得到错误的 SHG 抵消结果。

从 SM (S10) 的 sum rule 选取包含 $1/\omega_{nn'}$ 或 $1/\omega_{m'm}$ 的虚跃迁项，再用 $v_{pq}^a=i\omega_{pq}r_{pq}^a$，在 $\delta\to0$ 的领先阶得到 SM (S12)：

$$
\boxed{\;
\widetilde{D^a r_{nm}^{b}}
=i\sum_{n'}r_{nn'}^a r_{n'm}^b
-i\sum_{m'}r_{nm'}^b r_{m'm}^a.\;}
\tag{2}
$$

其中 $n'$、$m'$ 是与对应带准简并的搭档；在三带模型中 $c\leftrightarrow c'$，价带 $v$ 没有搭档。波浪号只保留 $1/\delta$ 型可能发散的领先项；有限余项不计入。

例如，根据 (1)–(2)：

$$
\widetilde C_{vc}^{abc}=i r_{vc}^b r_{cc'}^a r_{c'v}^c,
\qquad
\widetilde C_{cv}^{abc}=-i r_{cv}^b r_{c'c}^a r_{vc'}^c,
\tag{3}
$$

与 SM (S13) 一致。因为 $r_{cc'}^a\sim1/\delta$，两项分别可能发散。对于 $c'$ 也有交换带标记得到的对应关系。**导带相互接近并不令 $v\leftrightarrow c$ 的真实光学跃迁能量趋零**；这个奇异性来自中间的 $c\leftrightarrow c'$ 虚跃迁。

## 2. 将 SHG 拆成两个贡献

SM (S5)–(S7) 采用

$$
P^\alpha(2\omega)=\left[\chi_e^{\alpha\beta\gamma}(-2\omega;\omega,\omega)
+\chi_i^{\alpha\beta\gamma}(-2\omega;\omega,\omega)\right]E^\beta E^\gamma.
$$

$\chi_e$ 为显式三带 interband 项，$\chi_i$ 为 interband–intraband mixed 项。省去两者共有的 $e^3/(2\hbar^2)\int d\mathbf k/(2\pi)^d$，完整 interband 被积表达式（SM S6）为

$$
\chi_e\propto\sum_{nml}
\frac{r_{nm}^{\alpha}(r_{ml}^{\beta}r_{ln}^{\gamma}+r_{ml}^{\gamma}r_{ln}^{\beta})}{\omega_{ln}-\omega_{ml}}
\left[\frac{2f_{nm}}{\omega_{mn}-2\omega}
+\frac{f_{ln}}{\omega_{ln}-\omega}
+\frac{f_{ml}}{\omega_{ml}-\omega}\right].
\tag{4}
$$

SM (S7) 中可能携带 $1/\delta$ 发散的部分，仅来自 Hermitian connection；其余含 $g_{nm}$、$\Omega_{nm}$、$\Delta_{mn}$ 的项在当前**价带–导带光学对的间隙保持有限**的条件下不含这种中间导带分母。保留奇异部分得到 SM (S19)：

$$
\begin{aligned}
\widetilde\chi_i\propto i\sum_{nm}f_{nm}\Bigg[&
\frac{2(\widetilde C_{nm}^{\gamma\alpha\beta}+\widetilde C_{nm}^{\beta\alpha\gamma})}
{\omega_{mn}(\omega_{mn}-2\omega)}\\
&+\frac{\widetilde C_{mn}^{\gamma\beta\alpha}+\widetilde C_{mn}^{\beta\gamma\alpha}}
{\omega_{mn}(\omega_{mn}-\omega)}\\
&-\frac{\widetilde C_{mn}^{\alpha\gamma\beta}+\widetilde C_{mn}^{\alpha\beta\gamma}}
{2\omega_{mn}(\omega_{mn}-\omega)}\Bigg].
\end{aligned}
\tag{5}
$$

必须同时保留 (4) 和 (5)。若只计算 $\chi_e$，单独出现的 $r_{cc'}\sim1/\delta$ 不能当作最终 SHG 的物理奇异性。

## 3. $2\omega$ 分母：SM (S20)–(S21)

把 (4) 中的 $2f_{nm}/(\omega_{mn}-2\omega)$ 和 (5) 的第一项按相同的 $n,m$ 整理。利用 $\omega_{n'n}\to0$ 时 $\omega_{n'n}-\omega_{mn'}\simeq-\omega_{mn}$（SM S20），得到可能发散的共同括号：

$$
\frac{2f_{nm}}{\omega_{mn}(\omega_{mn}-2\omega)}
\Bigg\{i\widetilde C_{nm}^{\gamma\alpha\beta}
+r_{nm}^{\alpha}\left[
-\sum_{n'}r_{mn'}^{\beta}r_{n'n}^{\gamma}
+\sum_{m'}r_{mm'}^{\gamma}r_{m'n}^{\beta}
\right]+(\beta\leftrightarrow\gamma)\Bigg\}.
\tag{6}
$$

令方括号中的双和式为 $B_{mn}^{\gamma\beta}$。由 (1)–(2) 可直接写出

$$
\widetilde C_{nm}^{\gamma\alpha\beta}
=i r_{nm}^{\alpha}B_{mn}^{\gamma\beta}
\quad\Longrightarrow\quad
 i\widetilde C_{nm}^{\gamma\alpha\beta}=-r_{nm}^{\alpha}B_{mn}^{\gamma\beta}.
\tag{7}
$$

将 (7) 代回 (6)，$\beta\gamma$ 两种顺序都分别归零：

$$
\boxed{\widetilde\chi_{e,2\omega}+\widetilde\chi_{i,2\omega}=0.}
\tag{8}
$$

这里的等号仅针对 $\delta^{-1}$ 领先项；在舍弃 $O(\delta/\Delta)$ 分母修正时保留的有限部分不包含在内。

## 4. $\omega$ 分母：SM (S22)–(S26)

从 (4) 移除已经处理的 $2\omega$ 部分后，仅剩 $f_{ln}$ 与 $f_{ml}$ 两组。由于可能发散的内连接只出现在 $c,c'$ 之间，三带求和有两类配对关系：

**第一类：** $f_{ln}$ 取 $m=l'$，$f_{ml}$ 取 $n=l'$（SM S23）。用准简并分母近似，调换哑指标并应用 (2)，SM (S24) 给出

$$
\widetilde\chi_{e,\omega}^{(A)}
=\sum_{nm}\frac{i f_{mn}}
{\omega_{mn}(\omega_{mn}-\omega)}
\left(\widetilde C_{mn}^{\gamma\beta\alpha}
+\widetilde C_{mn}^{\beta\gamma\alpha}\right).
\tag{9}
$$

而 (5) 第二项与 (9) 仅相差 $f_{mn}=-f_{nm}$，故

$$
\widetilde\chi_{e,\omega}^{(A)}+
\widetilde\chi_{i,\omega}^{(A)}=0.
\tag{10}
$$

**第二类：** $f_{ln}$ 取 $m=n'$，$f_{ml}$ 取 $n=m'$（SM S25）。同理调换指标并应用 (2)，SM (S26) 得到

$$
\widetilde\chi_{e,\omega}^{(B)}
=-\sum_{nm}\frac{i f_{mn}}
{2\omega_{mn}(\omega_{mn}-\omega)}
\left(\widetilde C_{mn}^{\alpha\gamma\beta}
+\widetilde C_{mn}^{\alpha\beta\gamma}\right).
\tag{11}
$$

它与 (5) 的第三项再次因 $f_{mn}=-f_{nm}$ 相消：

$$
\widetilde\chi_{e,\omega}^{(B)}+
\widetilde\chi_{i,\omega}^{(B)}=0.
\tag{12}
$$

(8)、(10)、(12) 覆盖 SM (S19) 中全部可能含 $1/\delta$ 的分母结构，因此

$$
\boxed{\operatorname{dvg}(\chi_e^{\alpha\beta\gamma}+\chi_i^{\alpha\beta\gamma})=0.}
\tag{13}
$$

证明的条件包括：导带双重准简并、有限 valence–conduction 间隙、相同带标与占据、对 $\beta\gamma$ 的对称化、将完整 SHG 两个通道相加、频率处于所讨论的普通有限 $\omega$ 区域。准确简并点应改用非 Abelian 表达，不能在那里单独计算 $r_{cc'}$。

## 5. 用 SM 紧束缚模型进行独立验证

SM Sec. IV.5 的三轨道模型是

$$
H(\mathbf k)=\begin{pmatrix}
 A&-it\sin k_x&-id_1\sin k_y\\
 it\sin k_x&-A&-id_2\sin k_x\\
 id_1\sin k_y&id_2\sin k_x&E_g+\epsilon_1(C_0+C_1\cos k_y)
\end{pmatrix},\quad
 A=C_0+C_1\cos k_y+C_2\cos k_x,
\tag{14}
$$

参数 $C_0=-0.75$、$C_1=0.5$、$C_2=1$、$t=1$、$d_1=d_2=0.3$、$E_g=-6$，通过 $\epsilon_1$ 调整导带/价带的光学速度比 $v_v/v_c$。$M_x=\mathrm{diag}(1,-1,1)$；自旋忽略时 $H(-\mathbf k)=H(\mathbf k)^*$。

| $\epsilon_1$ | SM Table S1 $v_v/v_c$ | 独立计算 |
|---:|---:|---:|
| 0.75 | 0.701 | 0.700498 |
| 0.90 | 0.849 | 0.849232 |
| 1.03 | 0.978 | 0.978166 |
| 1.50 | 1.445 | 1.444550 |
| 2.00 | 1.941 | 1.941068 |
| 3.00 | 2.935 | 2.935042 |

在 $\epsilon_1=1.03$、$k_x=q$、$k_y=k_{y0}+0.37q$、$\hbar\omega=0.7\epsilon_0$ 处逐点计算 SM (S6) 的 $\chi_e$、SM (S19) 的 $\widetilde\chi_i$。公共物理前因子省略，数值仅用于检查奇异阶次：

| $q$ | $\lvert\chi_e\rvert$ | $\lvert\widetilde\chi_i\rvert$ | $\lvert\chi_e+\widetilde\chi_i\rvert/\lvert\chi_e\rvert$ |
|---:|---:|---:|---:|
| $10^{-3}$ | $6.26\times10^{-3}$ | $6.25\times10^{-3}$ | $5.14\times10^{-3}$ |
| $10^{-4}$ | $6.25\times10^{-2}$ | $6.25\times10^{-2}$ | $5.16\times10^{-4}$ |
| $10^{-5}$ | $6.25\times10^{-1}$ | $6.25\times10^{-1}$ | $5.16\times10^{-5}$ |

**相同的检查单独应用于 $2\omega$ 和 $\omega$ 两组，均得到有限余项。** 数值验证还包括：完整 SM (S10) 与平行输运规范下的 $k$ 差分一致；$\widetilde C$ 与完整 $C$ 的相对差随 $q\to0$ 线性下降；U(1) 随机相位、$c/c'$ 互换、$\beta/\gamma$ 互换不改变结果。

## 6. 与 shift current 的关键区别和严格边界

SM (S15)–(S17) 的 shift current 含有 $\operatorname{Im}\widetilde C_{vc}$ 等共振项，两个导带贡献在特定偏振分量下**可同号累加**；SHG 则必须将 (4) 和 (5) 一起求和，其三组领先奇异项成对消失。因此，DBC 对 shift current 的增强机制不能直接推广为 SHG 必然发散。

本复现证明了 **$1/\delta$ 的消除**。它没有计算 SM (S7) 中全部有限项，也没有给出完整 SHG 频谱；有限峰值、精确共振、SOC、HSE、GeTe 和实验可观测量需要另行检验。对于磁 shift current，正文/SM 另有结论，不应从这里的 SHG 抵消推出其消失。

## 7. 运行与检查

在独立论文目录下：

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/reproduce_core.py --outdir results
```

单独的返回 JSON 和 CSV 留下各 $q$ 的数值，图文件在 `results/shg_cancellation.png`。论文原文及 SM PDF 不纳入 GitHub 仓库；所有图为自算结果。

