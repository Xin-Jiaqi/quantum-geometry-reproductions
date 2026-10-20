# 第四轮｜SM Tables I–III 与 Group Theory Operations Toolkit 全 122 MPG 审计

2026-10-19，目标：Xiang et al., PRL 134, 116301 (2025) 正文 Eqs. (11)–(14)、SM §V pp.8–10。仅复现论文原有的磁点群对称性，不拓展层结构或新材料。

## 1. 真正接入了哪个 GitHub 来源

用户维护的公开仓库 `Xin-Jiaqi/group-theory-operations-toolkit`，0.18.0，固定 Git commit `551a0daf8d4ba5a5f24b2811694bfe438616de48`。通过 GitHub connector 读取实际主分支（不是搜索引擎推测）的 `data/magnetic_point_groups.json` 与 `data/group_operations.json`，源 Git blob SHA 分别为 `406309c2686cbdaf7634dbea10ad8b059ee40e4d`、`1c43c66646e214c8f5939c06a17a8e69fb4c5122`。原仓库已具有 122 种 MPG 的反幺正标签和通用 `magnetic_tensor_basis()` API。

因为当前执行容器不能联网 `git clone` / `pip install`，只将相关的 **122 个磁点群生成元、三套笛卡尔矩阵及预期群阶数**制作成明确版本固定的轻量 snapshot（`data/toolkit_mpg_generators_v018.txt`）；没有假装已在容器里安装并运行原仓库。`src/igmc_lab/toolkit_bridge.py` 在这些真实生成元上**重新实现独立的 Neumann 约束矩阵和 SVD 零空间**，不是复制群论 toolkit 的张量代码。它额外提供可选的 live 接口，以用户安装 toolkit 后用原 API 逐个比较 244 个 `dimension`。

源 group 数据含 32/32/58 种 I/II/III 型 MPG，共 122；本次逐一重建所有群生成元闭包与阶数，没有发生错误。小数据 snapshot 与溯源 SHA 均作为本论文的固定输入，**不触碰原用户的群论仓库**。

## 2. 物理表象与时间反演约定

设两种实际测量 quadrature 为

\[
J_\alpha(t)=(-g\mu_B)\bigl[\sigma^C_{\alpha\beta} B_\beta\cos\omega t
+\sigma^D_{\alpha\beta} B_\beta\sin\omega t\bigr]\quad\text{(sign convention of paper recorded separately)}.
\]

这里仅用于对称性判定：`J` 是时间奇的极矢量，`B` 是时间奇的轴矢量，`dB/dt` 是时间偶的轴矢量。conduction C 对时间反演为偶系数，displacement D 对时间反演为奇系数。令空间矩阵 R，时间反演标签 θ=0/1，二维秩张量空间变换

\[
T'=\eta_X^\theta\det(R)R T R^{\mathsf T},
\quad \eta_C=+1,\quad\eta_D=-1.
\]

这是实 in-phase / quadrature **系数**的静态磁群 Neumann 分类；并不是对任意复频率 \(\sigma(\omega)\) 忽略复共轭、频率反号的完整反幺正关系。正文 Eq.(14) 用同一 `-1` 作用于两类系数，另保留 `paper_printed_eq14` 分支，不暗中修改原文。

用户 toolkit 原生接口与这里的对应是：

```python
from group_theory_operations import magnetic_tensor_basis
C = magnetic_tensor_basis('4mm1\'', 'polar_vector', 'axial_vector',
                          output_time_parity='odd', input_time_parity='odd')
D = magnetic_tensor_basis('4mm1\'', 'polar_vector', 'axial_vector',
                          output_time_parity='odd', input_time_parity='even')
# C.dimension == 1 ; D.dimension == 0
```

安装独立工具后调用 `python scripts/audit_round4.py --require-live-toolkit`，应逐项比较 122×2 = 244 个结果；当前仅跑过离线生成元独立求解，此项 live 双实现比对 **待执行**。

## 3. 全群分类结果

| 由群对称性允许的类型 | Conduction C | Displacement D |
|---|---:|---:|
| 仅纵向 \(\alpha=\beta\) | 24 | 23 |
| 仅横向 \(\alpha\ne\beta\) | 18 | 19 |
| 纵向与横向均有 | 16 | 16 |
| 全部禁止 | 64 | 64 |
| **合计允许** | **58** | **58** |

这与 SM Table III 的 3 类逻辑和代表性分组相符。完整 122 行请见 `results/raw/round4_mpg_122.csv`，其中每一行包含磁群编号、符号、群阶数、允许张量维数及纵/横向类别；完整可归一化的三维数值张量基在 `round4_mpg_tensor_shapes.json`，**不是只给几个手挑例子**。

全 122 群的 244 次 SVD 得到的每个基向量都做了逐生成元不变性残差检验（阈值 \(10^{-6}\)）。群闭包也按包括时间反演的 `(R,θ)` 群律独立验证。

## 4. 与 SM 的具体张量结构一致性

- `4mm1′`：C 仅反对称 `xy/yx`，1 个参数；D 完全被时间反演禁止。
- `4m′m′`：C 仅反对称 `xy/yx`，1 个参数；D 对角 `xx=yy` 与 `zz` 两个独立参数。
- `3m1′`、`6mm1′`：C Hall 型 1 个，D 为 0。
- `1′`：C 全 9 个，D 0；`-1`：C、D 均 0；`PT=-1′`：C 0，D 全 9 个。
- SM p.8 Table I 的单斜 `2` 写成唯一轴 y，允许 `xx,xz,yy,zx,zz`；toolkit 的标准生成元采用 `2_z`。显式交换 **y↔z** 后完全一致。`m` 的镜面张量结构也在相同坐标置换下完全一致。**不能将不同晶轴 setting 误判成原文错误。**

## 5. Eq. (14) 的扩展证据与限制

将正文 Eq.(14) 的统一 `-1` 当作 C 的反幺正变换因子时，在全 122 个标准 MPG 中有 **58 个群**的 C 张量空间**维数**与分通道规则不同。例如 `4mm1′` 会被错误判为 C=0，但 SM Table I 明确允许其非零 Hall 分量。这个结果扩大了前一轮单独 Rashba 反例的范围；它依赖本文已明确的**实系数 / quadrature 响应协议**。

不能仅凭此断定作者对一般复频率 antiunitary response 的 Eq.(14) 没有别的未声明约定。因此记录为 `SOURCE / PROTOCOL DISCREPANCY`，不是擅自修正论文。

## 6. 局限与可重复性

- `PASS`：基于固定 toolkit 生成元完成 **122 个磁点群、244 个张量空间**的独立对称性计算，群阶数闭包检查和若干 SM I–III 的逐结构测试。
- `PARTIAL`：SM Tables I/II 的每一个符号化参数关系（尤其非晶轴 setting）**尚未完成逐单元人工排版核对**，虽然全部 122 个群都有数值空间和参数维数表。
- `BLOCKED`：直接在当前离线容器里运行 `group_theory_operations.magnetic_tensor_basis` 需要从已安装的原仓库环境提供；本轮测试该项 `skipped`，未冒充通过。
- `OPEN`：Fig. 2(d) 原作者绝对幅度，与第三轮的因子 2 来源问题无关，本轮未强行关闭。

```bash
# 当前包本地，完全离线可跑
PYTHONPATH=src python scripts/audit_round4.py
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/verify_release.py

# 另在有网络的开发环境中（可选，把 toolchain pin 在单独环境）
python -m pip install 'git+https://github.com/Xin-Jiaqi/group-theory-operations-toolkit.git@551a0daf8d4ba5a5f24b2811694bfe438616de48'
python scripts/audit_round4.py --require-live-toolkit
```

保存标准：`data/toolkit_mpg_generators_v018.txt` / `results/raw/round4_mpg_122.csv` / `results/raw/round4_mpg_audit.json`, `results/raw/round4_mpg_tensor_shapes.json` / `results/frozen/round4_mpg_reference_v2.json`。源仓库版本一旦改变，应新建 snapshot，而不是悄悄覆盖这轮固定数据。

## 7. 全 244 个张量结构也可重建，而不只是维数

为避免 SVD 返回的随机正交基妨碍核对，`canonical_component_basis` 从每个 SVD 零空间中选择一组线性独立的张量分量作为自由参数，并通过矩阵求逆把这些自由分量标准化为 1，其余分量由群对称性给出（通常为 0、+1、−1 或代数比例）。所有标准化后的 244 个张量子空间再次通过逐生成元不变性检验。

`results/raw/round4_mpg_tensor_shapes.json` 包含**122 种 MPG × C/D 两类的完整数值张量基**，便于逐项核对 SM Table I/II 的符号参数；其 SHA-256 也在 `results/frozen/round4_mpg_reference_v2.json` 冻结。这是已完成的机器版完整结构表，而**逐项人工核对排版、磁群标准 setting 和参数符号**仍为下一轮。无需逐个猜测作者表格中的哪一项可以非零。

## 8. 完全不同的第二种群论维数核验

除 9×9 实矩阵约束的 SVD 零空间以外，还用磁点群**完整闭包**计算 Reynolds 平均中的特征标：

\[ d_X=\frac{1}{|M|}\sum_{(R,\theta)\in M} \eta_X^\theta\det(R) [\operatorname{Tr}R]^2. \]

独立地在全部 122 个 MPG、C/D 两个通道检查 `d_character == d_SVD`，**244/244 一致**；只把生成元列表换一条算法求秩不算这种验证。另对 source Eq.14 alternative 的 C 分支也作特征标检查。

Round 4 在加入这个新独立核验之前已经冻结一个早期审计版本 `results/frozen/round4_mpg_reference_v1.json`，加入 character check 和数值 tensor 结构以后新建 `round4_mpg_reference_v2.json`；不回写先前冻结版本。

## 8. 二次独立确认：磁群特征标方法（v2 checkpoint）

本轮对每个磁点群，不再仅使用 Neumann/SVD 的相同九维矩阵，而是独立遍历生成元闭包内的全部 `(R, θ)`，利用特征标平均公式

\[
 d_X=\frac{1}{|M|}\sum_{(R,\theta)\in M}\eta_X^\theta\det(R)\,[\mathrm{tr}(R)]^2,\quad
\eta_C=+1,\quad\eta_D=-1.
\]

其推导来自输出极矢量 `J` 的特征标 `tr R` 与输入轴矢量 `B` 的特征标 `det(R) tr R` 相乘，因此不需要复制 toolkit 的张量算法，也不调用 `tensor_matrix` 或 SVD。122×2 个结果**全部为整数、且与 SVD 张量空间维数一致**；相应第二组特征标维数也写入 CSV。额外将 `paper_printed_eq14` 的 conduction 分支独立进行 122 个特征标维数测试，和其 SVD 维数一致。

这比同一矩阵实现的两次调用更独立，但仍是同一组物理表象的两个数学实现，不可用作 Zeeman 量子几何公式或 Fig.2(d) 物理因子的独立验证。

原第四轮参考结构发生**可追溯的扩充**：新增 `C_char_dim / D_char_dim` 两列和 `character_crosschecks=244`，使用 `results/frozen/round4_mpg_reference_v2.json` 冻结。前几轮冻结模型、SI 与 Wilson/Liouville 数值 **没有改动**。
