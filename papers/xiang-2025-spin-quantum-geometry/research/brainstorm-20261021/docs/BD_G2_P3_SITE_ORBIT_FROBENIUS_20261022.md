# B+D G2 | p3 层群的 Wyckoff 位点、隐藏自旋响应与 Frobenius reciprocity

**记录日期：2026-10-22**。**科学验收：源对齐的幺正群代数 PASS；实际铁电切换/磁性模动力学/新效应均未证明。** 群论实验在 [group-theory-operations-toolkit Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32)，受保护 main 尚未合并。严格保留 MnPS3 R1 已关闭/未收敛，且不提交任何新计算。

## 1. 我们到底计算了什么？

参考 [AFLOW 2D plane group p3/#13 原始位置表](https://aflowlib.org/prototype-encyclopedia/Tutorials/Two_Dimensional_Periodic_Systems_Part4/)，对应 [International Tables 层群 p3/#65](https://onlinelibrary.wiley.com/iucr/itc/Ea/ch4o1v0001/sgtable4o1o065/)。约定常规六方基 \(\mathbf a_1=(1,0),\ \mathbf a_2=(-1/2,\sqrt3/2)\)；\(C_3\) 的分数表示
\[
R_3(x,y)=(-y,x-y).
\]
对比 **同一个 \(C_3\) 点群、相同磁性位点数目** 的两种 **明确不同的结构占位**：
- **A**：每层三个位点来自一个一般 \(3d\) 轨道（\((x,y),(-y,x-y),(-x+y,-x)\)），在 \(C_3\) 下三者循环交换；
- **B**：每层三个位点分占 \(1a=(0,0),1b=(1/3,2/3),1c=(2/3,1/3)\)，各自在 \(C_3\) 下按晶格平移等价。

构造双层六个位点的*合成*层 A 型反铁磁序（上三处磁矩 +z，下三处 −z），\(C_3\) 仍是磁性幺正对称。保证这是 **只按 C3 约束的受控 representation toy**，没有证明其它可能的镜面、额外平移、反幺正关系在真实样品中不存在。特殊三个中心若同元素且全部局域装饰相同，还可能出现更小原胞，必须用真实 ligand/位点环境确认独立性。它们**绝不是** MnPS3 的真实 L2 磁序和已认证 FE 双态。

输出是每个物理位点的局域电子 \(S_i^z\) 或横向 \((S_i^x,S_i^y)\)，输入是面内均匀**极性电场** \(E_x,E_y\)。使用真实位点 Seitz 置换矩阵 \(P(g)\)、自旋的轴矢量旋转和完全显式的约束
\[
D_{\rm spin/site}(g)K=K D_E(g).
\]
注意只是**允许空间**，不是物理 Kubo 响应、磁振子激发或实际电学探测强度。

## 2. 独立的 Hom 零空间和特征标结果

| 允许系数维数 | 轨道 3d | 轨道 1a+1b+1c |
| --- | ---: | ---: |
| 双层六个局域 \(S_z\) 观测量 | **4** | **0** |
| 单层三个局域 \(S_z\) 观测量 | **2** | **0** |
| 双层六个局域横向 \(S_{x,y}\) 观测量 | **8** | **12** |
| 单层三个局域横向 \(S_{x,y}\) 观测量 | **4** | **6** |
| 总和 \(S_z\)、每层求和的 \(S_z\) | **0** | **0** |

例如 generic \(3d\) 给 \(P_3\) 三循环的迹 \(\chi=(3,0,0)\)，special \(1a+1b+1c\) 为 \(\chi=(3,3,3)\)，均匀 \((E_x,E_y)\) 为 \(\chi_E=(2,-1,-1)\)。使用
\[
\dim \mathrm{Hom}_{C_3}(E,Q)=\frac13\sum_{g\in C_3}\chi_E(g)\chi_Q(g)
\]
与奇异值分解的**独立**零空间计算完全一致。 \(3d\) 的 \(S_z\) 可表现为 \((-1,\ 0.44868,\ 0.55132)\) 这种纯局域响应图样，**每层和都严格为零**；数值仅是对称性基归一化，不是材料的预测磁矩。

对真正具有磁动力学意义的横向 \(S_{x,y}\)，计算时必须把自旋 \(C_3\) 轴矢量旋转 \(R_3\) 与位点置换 \(P_3\) 做 Kronecker 积，**不能**将纵向的允许系数数目误说成磁振子模式数。

## 3. 更高收益的计算简化：Frobenius reciprocity

由于每个磁性 Wyckoff 轨道表示可写作局部稳定子的诱导表示，经典
\[
\boxed{
\mathrm{Hom}_G(V_E,\mathrm{Ind}_{H_r}^G V_{\rm spin,r})
\simeq \mathrm{Hom}_{H_r}(\mathrm{Res}_{H_r}^{G}V_E,V_{\rm spin,r})
}
\]
使得 **全层群局域允许响应筛选不一定要先构造所有高维站点排列矩阵**。可先读取 site stabilizer 和局域物理算符，计算小维数允许空间，再对候选构造全群响应张量。这个效率性质是**已知 Frobenius 互易**，不是新定理。我们用 [独立的稳定子计算脚本](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/site_frobenius_gate.py) 对 \(3d\) 与 \(1a/b/c\) 的 **纵向 (2,0)、横向 (4,6)** 逐一复核了完整位点置换结果。

**严格范围：** 目前只验证了 \(\Gamma\)、幺正 \(C_3\) 例子，未建立反幺正 magnetic layer-group corep 或有限 \(q\) 的 nonsymmorphic 诱导表示程序。

## 4. 本轮第二个重要 NO-GO：允许系数 ≠ 可激发模，也 ≠ 线偏振各向异性

任何 \(K\) 由两个面内电场分量输入，所以 \(\mathrm{rank}K\le2\)；纵然横向 Hom 维数从 8 变为 12，也未增加电场同时可激发的独立向量数。数值测试对两种横向轨道都产生 \(\mathrm{rank}K=2\)。

对真实归一化的**实**局域响应矩阵，\(C_3\) 等变约束还严格要求
\[
K^T K=\lambda I_2.
\]
在对称不变的等权局域观测范数下，\(\|K E_{\parallel}\|^2\) **与线偏振方向无关**；随机允许矩阵的残差约 \(10^{-15}\)。这不否定某个单独位点的方向性。完整**复数交流**响应则仅约束 \(K^\dagger K=aI+b\,iJ\)，可能允许圆偏振手性差异：纯表示基例子给两个 helicity 特征值 \(0,2\)。后者不证明真实吸收、材料光磁性或磁振子圆二色性，因为物理 Kubo、反幺正与观测权重均尚缺失；[Yuan et al. PRL 137, 106901 (2026)](https://doi.org/10.1103/12cl-b9jj) 等已有相近的磁振子光学理论。

## 5. 为什么这还不是 FE 切换故事？

**对称等价铁电双态不可能被用来任意改变 Wyckoff 轨道的稳定子阶数。** 真正父群连接操作必将稳定子共轭、保持轨道重数。因此严格磁性等价的 ±P 态，不可能仅凭这一操作使一个 \(3d\) 轨道变成三个独立 \(1a+1b+1c\) 特殊位点。因此“FE 切换使四个局域系数一键归零”**在此示例中是无效故事**。若是磁序被固定、SOC 破坏连接或发生真实非等价结构重构，应重新检查实际磁群与能量；这些问题这里尚未计算。

## 6. 2020–2026 特别近的同行研究

- [Watanabe & Yanase, *PRB* 98, 245129 (2018)](https://doi.org/10.1103/PhysRevB.98.245129)：磁性多极矩与电磁响应已有一般群论。
- [*Hidden physical effects in noncentrosymmetric crystals*, *PRB* 102, 165143 (2020)](https://doi.org/10.1103/PhysRevB.102.165143)：全局禁戒但局域允许的自旋、Berry 等 hidden effects 已有一般讨论。
- [Urru & Spaldin, *Annals of Physics* 447, 168964 (2022)](https://doi.org/10.1016/j.aop.2022.168964)：磁性八极多极矩与补偿的局域磁电响应。
- [Li et al., *Nano Letters* 25, 12762 (2025)](https://doi.org/10.1021/acs.nanolett.5c03407)：电场驱动三聚体局域磁矩的生成/湮灭，直接抢占“控制局域三站点磁序”的大方向。
- [Tan et al., *PRB* 113, 085115 (2026)](https://doi.org/10.1103/26wm-zhv5)：隐藏径向自旋极化与子晶格对称性，已有真实材料例子。
 
**同行层面结论：** 此 G2 发现属于**正确且高效的已知表示论实现**；单独以局域自旋 hidden effect 命名论文难以越过新颖性门槛。必须加上有物理证据的堆叠切换关系、新的可测响应/互易差别，才值得继续扩展全层群。

## 7. 下一轮低成本科学搜索，分支按可证伪性推进

1. **物理更严格的 FE 连接图：** 从 BSF 现有 80 LG 的结构连接数据中筛选真正的 ±极化关系，匹配实际磁性位点稳定子和 full SOC / antiunitary 映射，首先检查是否只是对应的局域 Kubo/Onsager 重标号。无 FE 来源的 \`3d\to1a+1b+1c\` 一律不得纳入“可切换器件”样本。
2. **真正的模式或 Kubo 响应：** 由一个真实来源的 magnetic site orbit + 物理可归一化的横向磁力坐标，建立 \(\Gamma\) 模表示/力常数；从实空间原子/交换构造一个少参数模型，不给自由谱参数追新奇；若对称性允许但全部微观矩阵元为零则 NO-GO。
3. **推广磁群的正确工具：** 稳定子程序可先扫 **unitary** 的 Wyckoff 轨道，随后才处理 antiunitary corep 与有限 \(q\) Bloch 平移相位。清晰区分 **80 ordinary LG、528 magnetic LG、448 collinear spin LG**，不把同一个群表的不同计数当发现。
4. **独立原创性门槛：** 只有能给出严格 FE 对应结构、无法被普通点群多极矩/隐藏自旋群论/Onsager 自动推出的测量差异，才进入全层群与实材料阶段。

**当前决定：** G2 群论数学 PASS；作为普适的“Wyckoff 能允许隐藏磁响应”属于已知物理；**暂不把全层群视作可投稿新理论**。保留代码和全部负结果作为 G3 决策基准。R1 不重启，R2-P0 只读。
