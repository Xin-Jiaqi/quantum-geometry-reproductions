# Round 5: SM Tables I–II 全表逐分量审计（2026-10-20）

## 结论

从用户提供的原论文补充材料 `supple-2(1).pdf` **第 8–10 页直接核对** Tables I/II 的磁点群分组和所印张量表达式，再与第四轮使用 toolkit v0.18.0（commit `551a0daf8d4ba5a5f24b2811694bfe438616de48`）生成元独立计算的张量子空间对照。**不能把机器预测反填成论文原始表格**。

- 来源表 23 行张量模板：Table I 11 行、Table II 12 行；一共覆盖 C/D 各 **58 个允许 MPG**。其他各 64 个禁戒 MPG 也显式对照零空间。
- 总共 122 MPG × 2 = **244** 组张量空间逐项检验：**243 完全一致，1 项源文张量参数的独立性不足**（详情如下）。这里“完全一致”指二维子空间的 **9×9 正交投影矩阵相等**，检查全部 9 个分量及相等/相反关系，不只是允许/禁戒数量或 tensor rank。
- Monoclinic `2`, `m`, 以及相应含撇磁群，SM 采用唯一轴 **y**，toolkit 的 snapshot 采用唯一轴 **z**。已明确在源矩阵上施行 `P[y↔z]`，保持轴约定一致后再做投影子空间比较；此转换不是为了拟合任意分量。
- 完整 Table I: **122/122**（包括 58 允许+64 禁戒）对照通过。Table II: **121/122**，具体异常 `-42'm'` (No. 52)。SM Table III 的纵/横向分类不受该异常影响。**Table III 排版逐符号抄录没有作为独立的第二组来源测试**，不应把 Table I/II 与 III 的重述算作两个独立源证据。

## 原文 Table II 的一个遗漏关系（p.9, 第 5 个分组）

SM Table II 将 `2'2'2, mm2, -42'm'` 放到同一行，展示模板

\[
\sigma^D_{\mathrm{SM,row\ 5}}=
\begin{pmatrix}0&A&0\\ B&0&0\\0&0&0\end{pmatrix}
\]

其中 A = σ_xy、B = σ_yx 是表内 **两个不同的参数符号**。对于 `2'2'2` 和 `mm2`，这确实是 2 维允许空间。**对于同一行的 `\bar4 2' m'`**，四次反轴操作还要求 A = B。因此它是上面 2 维模板的 **1 维子空间**，具体为

\[
\boxed{\sigma^D_{\bar4 2' m'}=a
\begin{pmatrix}0&1&0\\1&0&0\\0&0&0\end{pmatrix}}
\]

解析推导可独立于数值 SVD：对该群的生成元选取

\[
S_4=\begin{pmatrix}0&1&0\\-1&0&0\\0&0&-1\end{pmatrix},\qquad
(C_{2x}T): R=\operatorname{diag}(1,-1,-1),\quad\theta=1.
\]

对 D 的实 quadrature 约定采用 \(\sigma^D=\eta_D^\theta\det(R)R\sigma^D R^T\)、\(\eta_D=-1\)。设面内子块 \(\bigl(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\bigr)\)。

- 无撇 \(S_4\) 将其送到 \(\bigl(\begin{smallmatrix}-d&c\\b&-a\end{smallmatrix}\bigr)\)，因此 \(a=-d\)、\(b=c\)。
- \(C_{2x}T\) 将面内块送到 \(\bigl(\begin{smallmatrix}-a&b\\c&-d\end{smallmatrix}\bigr)\)，因此 \(a=d=0\)，\(b,c\) 尚不受约束。
- 两个条件合并后，恰有一个自由参数 \(b=c\)。另对 z 耦合分量施加生成元可知它们全部消失。数值结果的 9D projector 与这个解析唯一基完全吻合。

**界限：** 原文的非零分量位置、允许 Hall/transverse 分类没有错；缺的是属于这个特定 MPG 的 **\(\sigma_{xy}=\sigma_{yx}\)** 约束。不能把此处说成整张 Table II 错，也不能改变 Table III 允许群数。可能是作者用一个共同容许形状归组，但 Table II 标题宣称各群的 tensor shape，故这里至少需要独立参数说明。

## 本轮没有解决的科学问题

1. 正文 Eq.(14) 统一的 \(\eta_T=-1\) 反幺正因子，与来源 SM Table I 的纯 T 群及第四轮独立 C 响应协议仍有歧义：**源公式不等价于直接解释出来的源分类**，需分析复频率反幺正响应的定义才能定性为笔误。
2. SM 第 6 页 Eq.(34) 的 \(Q\) 前因子与正文 Eq.(5) 的显式定义不同；第三轮的 Wilson-loop、直接 Liouville/Kubo 证据支持正文 Eq.(12) 分支，**不能因为 Fig.2(d) 更接近被乘二的诊断值就修复图**。
3. 目前仍不能在该离线容器直接运行第三方原版 `group_theory_operations` Python API。第四轮的源 toolkit 数据是可信的固定生成元 snapshot，方法级独立核验存在，但 **native live toolkit integration = BLOCKED/optional skip**。

## 可复现输入与验收

源文件：`src/igmc_lab/paper_table_source.py`（由 SM 原件逐行人工辨认的张量模板；含页码与 23 个分组，不从计算结果导出）。

```bash
cd papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
PYTHONPATH=src python scripts/audit_round5.py
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/verify_release.py
```

结果：`results/raw/round5_source_table_244.csv`（244 行，包括禁戒）、`round5_source_table_audit.json`、新冻结 `results/frozen/round5_source_table_reference_v1.json`。旧 Round1–4 冻结结果的 SHA-256 保留且不重算替换。此处报告的是 **symbolic source-vs-model audit**，不意味着作者的 Fig.2(d) 数值验收。
