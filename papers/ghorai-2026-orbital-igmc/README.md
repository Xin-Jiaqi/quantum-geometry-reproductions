# Orbital IGMC: 核心物理接入

Ghorai, Sarkar & Agarwal, *Intrinsic gyrotropic magnetic current of orbital origin*, **PRB 114, 245420 (2026)**. DOI: [10.1103/xr5v-gkv5](https://doi.org/10.1103/xr5v-gkv5).

**MODEL_READY 的范围：指定模型、明确公式协议的可运行复现。** 已接通交变磁场引起的轨道极化及其位移电流、CuMnAs 的 PT 双重简并、带隙内响应和 Néel 反号。原文 D10 与 A9 的轨道幅度存在未解决差异；两入口分别保留，调用时必须指定，不能据此宣称唯一物理幅度已核准。

## 运行与 API

仓库根目录：

```bash
python scripts/core.py test orbital-igmc
python scripts/core.py run orbital-igmc --output outputs/orbital-igmc.json
```

本目录独立运行（NumPy；跨论文测试另需原仓库及 SciPy）：

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/reproduce_core.py
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/plot_core.py
```

```python
from orbital_igmc import CuMnAs, geometry, integrate, mesh, response_si

model = CuMnAs(a_angstrom=3.82)  # 显式示例输入；主文未给面内晶格常数。
g = geometry(*model.operators(mesh(161)), group_sizes=(2, 2))
r = integrate(g, mu_eV=0.0, temperature_K=25.0, protocol="printed_D10")
chi = response_si(r, model, drive_energy_eV=0.0001, convention="harmonic")
# 改用 protocol="density_A9" 可直接求解原文密度矩阵驱动。
```

`geometry` 输入固定正交轨道基下的 H、两个 dH/dq、自旋磁矩及明确的能量块。`integrate` 返回 [xz, yz] 两分量；`printed_D10` 另外返回 metric、moment、multiband 三项。metric 梯度使用完整周期 BZ 上的分部积分，不能将其输出当作逐 k 的位置修正。

## 核心机制与证据

磁场改变电子位置，形成 P(t)；磁位移电流满足 j=dP/dt，对应原文 Eqs. (30)–(32)、(D9)–(D11)。它在绝缘体中仍可非零，静态极限电流为零。本文模型取 t0=0.08、t1=1、alpha=0.8、Jn=0.85 eV；原仓库 NMEE 的默认 Jn=0.6 eV 不能直接沿用。

原始数组和验证结果见 [results/core.json](results/core.json)，脚本每次重新计算。下表取 phi=pi/12、mu=0、T=25 K、hbar omega=0.1 meV、显式 a=3.82 Å、161×161 BZ：

| +omega 谐波 Im chi_xz，mA/(m T) | printed_D10 | density_A9 |
|---|---:|---:|
| 轨道 | 1.03411 | 0.517055 |
| 自旋 | -0.462693 | -0.462693 |
| 总和 | 0.571417 | 0.0543619 |

这些数值是独立计算结果，未用作者曲线拟合。论文未公开原始数据，面内 a 未明确给出，因此不标记为原 Fig. 3 的定量验收。

测试覆盖：Eq. (38) 能谱、解析速度/有限差分、BZ sewing、任意 U(2) 双重简并基旋转、轨道基幺正变换、Néel 反号、镜面对称性、恢复 T 后消失、绝缘体平台、空/满占据、DC 为零、线性频率与 SI 长度标度、网格收敛、密度矩阵/有限差分投影对照、原 D10 逐指标求和，以及现有 NMEE/IGMC 模块的接口交叉检验。共 23 个测试，其中 2 个需要完整仓库。

## 两个必须显式的约定

1. **公式协议。** `printed_D10` 忠实保留三项；`density_A9` 直接计算磁驱动密度矩阵及极化，不通过对 D10 乘系数得到。CuMnAs 带隙内，D10 的轨道值为 A9 的两倍，自旋一致；一般三带模型不适用统一倍数修正。详见 [SOURCE_AUDIT.md](SOURCE_AUDIT.md)。
2. **谐波与单位。** B=B0 cos(omega t)，`harmonic` 定义 j=chi B0 exp(+i omega t)+c.c.；`phasor` 定义 j=Re[chi B0 exp(+i omega t)]，后者为前者两倍。q=k a；H 与 dH/dq 用 eV；b=B/Bunit，Bunit=hbar/(e a²)。输出为二维面电流/磁场 A/(m T)，不自动除层厚或推断电压。

## 边界

仅低频本征磁极化通道；需 omega tau≪1、hbar/tau≪带间隙及 hbar omega≪带间隙。代码不提供真实碰撞积分、共振/无序输运、完整电磁电路或 0.2 mV 的实验论证。零温金属需要独立费米线积分，目前明确拒绝；有限温金属需另做网格收敛。

CuMnAs 使用持续简并的双重态子空间，并对近带交叉或块内非标量速度拒绝计算。一般固定轨道基三带测试只验证声明的算符协议，不能自动扩展到包含未知位置连接/自旋算符的 Wannier 材料。没有修改其他论文的冻结结果，也没有把本目录当作一般磁输运位置修正的替代品。

原论文 PDF/SM 不随仓库分发；本目录是独立实现，非作者官方代码。许可证不继承其他论文子目录的声明。
