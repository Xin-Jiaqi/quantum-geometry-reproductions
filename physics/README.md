# 当前物理能力目录

状态依据：独立求解或解析对照、模型参数检查、数值收敛和运行入口。`MODEL_READY`只适用于所列协议；新模型仍需检查自己的适用条件。[清单](capabilities.json)是根运行器和CI共用的数据源。

| ID | 直接API（在相应论文目录导入） | 核心证据 | 输入/输出 |
|---|---|---|---|
| `nmee` | `nmee_lab.evaluate_dynamic`；`evaluate_clean`；`evaluate_insulating_static` | 轨道N²求解、有限k微分、PT投影、六带、规范检查 | H/eV、k/nm⁻¹、同规范位置与自旋；AC为s/μB每E·B |
| `spin` | `methods.response.SpinModel`；`response_at_k`；`integrate_response` | SM四源项、独立轨道求解、非零四带及规范检查 | H、位置连接、自旋；一/二阶spin current或polarization |
| `magnetotransport` | `qgmap.orbital_phe_v4.orbital_phe_eq56`；`orbital_mnhe_v5.orbital_mnhe_3d`；`orbital_geometry_v54.band_geometry` | 有限B微分、四带spin MNHE、几何投影对照 | 选定spin/orbital通道；自然单位和SI桥分开 |
| `curvature-flux` | `flux_lab.core.qgt_eigen`；`causal_response.interband_kubo`；`source_kme_expansion` | 本征态/投影子、Kubo/Liouville、互易性、物理反幺正检查 | 同规范H、速度与物理Pauli自旋；几何/带间GHE/KME |
| `igmc` | `igmc_lab.core.pair_kernels`；`dirac_response`；`fermi_contour.integrate_contour` | 波包/格林函数/Liouville、费米海/费米面 | 两带物理Pauli基底、非共振频率；C、D及A/(m·T)换算 |

详细模型设置、打印公式分支与本轮例子见[根目录索引](../README.md)。`python scripts/core.py run <ID>`从根目录执行；它的作用是选择运行环境和记录来源，不改变物理公式。

## 已有方法扩展，按需进入

这些底层工具可复用，但不替代上述论文响应的完整验收。

| 方法 | 所在位置 | 使用边界 |
|---|---|---|
| 任意非简并N带metric/Berry/orbital moment | [QGMap v54](../papers/xiang-2026-quantum-geometric-map/docs/INDEPENDENT_MULTIBAND_ORBITAL_GEOMETRY_V54.md) | 同基底H、导数；moment×curvature子项不等于完整MNHE |
| 简并簇带间几何与响应 | [Cui nonabelian_response](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/nonabelian_response.py) | 声明簇展平和带间排除；不是完整金属响应 |
| 122磁点群及C/D张量空间 | [IGMC toolkit_bridge](../papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/src/igmc_lab/toolkit_bridge.py) | 明确通道时间反演性质与坐标setting |
| lead/NEGF、Ward、静态和有限AC模型工具 | [QGMap历史方法目录](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/d65fc43/physics/README.md) | 各自有限模型与边界；不是统一bulk/device输运理论 |
| Wannier H/position/spin导入检查 | [Cui模块](../papers/cui-2026-curvature-flux-zeeman-qgt/README.md) | 文件存在不等于物理算符来源已确认 |

## 未接入不应被误标为已实现

噪声关联函数、非线性光学涨落—耗散关系、交替磁三阶Hall、超导BdG Hall等需要各自的观测量、模型与验收。已有几何、密度矩阵或NEGF基础只是可复用部件。见[文献覆盖清单](../docs/LITERATURE_COVERAGE.md)。

旧`integration/physics_interface/contracts.evaluate`仍是设计占位，不调用它获取物理数值。
