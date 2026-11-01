# IGMC 当前核心调用（2026-11-01）

**正文主图：** [图片、数据与重算说明](figures/README.md)

**已接通论文两个模型的核心物理。** 传导电流C、位移电流D、完整非共振频率响应及SI换算均可调用；无需先完成材料、无序或作者全部图像复现。

| 入口 | 输入/输出与边界 |
|---|---|
| `core.rashba_analytic` / `rashba_integrate` | 连续Rashba模型传导电流；解析式与极坐标积分 |
| `core.dirac_response` | 晶格Dirac模型；返回C、D/Ω、D、完整复Kubo |
| `core.pair_kernels` + `kernel_response` | 可换H和速度的两带内核；固定物理Pauli自旋基底，非简并、非共振 |
| `fermi_contour.integrate_contour` | 论文Dirac金属精确费米边界积分；避免粗网格假收敛 |
| `integrate_conduction_fermi_line` | 独立Eq13费米面核对 |
| `energy_from_frequency` + `physical_response_si` | 明确Hz/rad·s⁻¹、g和Å单位；输出论文绘图约定的电流/磁场，A/(m·T) |

最短重算入口（仅需NumPy）：

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/reproduce_core.py
```

本轮[结果](results/core_access.json)采用g=2、f=1 THz、T=0。Rashba平台3.88457 A/(m·T)；Dirac μ=0时C≈0、D=0.0177538932 A/(m·T)；μ=1 eV时C=-33.4862352、D=-0.0641784473。金属C的费米海/费米面差8.37e-10 Å⁻¹；格林函数/波包逐k差1.24e-13；轨道Liouville/谱求解积分差4.36e-17 Å⁻¹。

这里D是低频一阶系数。完整1 THz的Im K给出0.0177560969 A/(m·T)。**旧Round8把两者差异归因于101²粗网格并不准确：**101²与201²的低频D只差约9e-14 A/(m·T)，而完整频率与低频值相差2.20e-6。这是有限频率修正，不能混称积分误差。旧冻结数值保留。

## 正文与SM独立回查

本轮重新逐页查看SM p6 Eq34、p7 Eq46/47：Eq34的Q=2Re(rσ)与正文Q=Re(rσ)不一致；Eq46的cosine变换每支应为π，而Eq47恢复了正确的1/2谐波。我们的计算按两支B/2和正文Q归一化，多方法相互核验。**这些局部不一致不能证明作者Fig2(d)绘图代码用了两倍系数。** 作者原始曲线缺失，幅度原因仍未知。

Eq14反幺正约定保留为局部未决；实际调用区分C与D的频率相位/时间反演性质，已有122群闭合、244个张量空间及独立特征标检验。不要用一个共同符号覆盖两通道。

## 运行环境

系统Python里的`group_theory_operations`只剩空命名空间，缺`magnetic_tensor_basis`。这是可选外部群论交叉检查的安装问题，不是IGMC求解器缺失。核心例子不依赖它；干净Python3.11环境完成76项测试（75通过、可选外部检查1项跳过）。完整本地toolkit源码可通过`PYTHONPATH=src:/path/to/group-theory-operations-toolkit/src`显式载入，避免改动全局环境。本轮Python3.12显式载入本地toolkit 0.16.0后，76项全部通过，包括244个外部API维数比较；Python3.11的全部历史冻结核验也通过。原122群离线检查始终运行。

历史README、验收文件和冻结清单保持原样；本页及根目录RESEARCH_HANDOFF是当前状态入口。本轮不加入新的散射或材料任务。
