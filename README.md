# Quantum Geometry · 可复用物理能力库

[![Core regression](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/workflows/core-capabilities.yml/badge.svg)](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/workflows/core-capabilities.yml)

从论文中提取、独立核验并复用量子几何与输运计算。当前 **7 个核心模块可运行**，每个都有明确的输入、频率/单位约定、独立对照和可重算例子。

`MODEL_READY` 表示指定模型与响应协议可用；原图完全复现、通用材料能力和源文争议分别记录。本仓库不是作者官方实现。

**看图判断能力：** [五篇正文主图图册](docs/MAIN_TEXT_FIGURES.md)——能带、响应扫描、几何分布；附重算脚本与原始数组。以物理机制可调用为标准，逐图标明公式分支和未决差异。

## 从物理问题选择入口

| 我想计算 | 模块 | 当前能力 | 主要边界 |
|---|---|---|---|
| 电场×磁场诱导自旋磁化 | [NMEE](papers/jia-2026-nonlin-magnetoelectric-edelstein/README.md) | 多带/PT简并、完整混合AC、静态/印刷公式分支 | AC结果不自动等于本征DC |
| 电场诱导自旋极化/自旋流 | [Spin](papers/xiang-2025-spin-quantum-geometry/README.md) | 一阶/二阶、物理自旋算符、完整有限频率响应 | 有限η不自动等于τ⁰ |
| 平面Hall与磁非线性Hall | [Magnetotransport](papers/xiang-2026-quantum-geometric-map/README.md) | 指定模型的spin/orbital PHE、MNHE及多带几何 | 未完成一般多带轨道位置修正 |
| Zeeman几何、Dirac通量及互易响应 | [Curvature–Flux](papers/cui-2026-curvature-flux-zeeman-qgt/README.md) | 四类几何、绕数/通量、带间GHE/KME | Table I全表、完整金属耗散未核准 |
| 交变磁场诱导传导/位移电流 | [IGMC](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/CORE_ACCESS.md) | Rashba/Dirac、非共振响应、SI换算 | 原Fig2(d)幅度来源仍未决 |
| 交错磁体磁性光伏：线偏振 MIC | [Ezawa BPVE](papers/ezawa-2025-altermagnet-bpve/README.md) | quantum metric × 带间速度差、磁序反号、Eq. (17) 数值基准 | 仅两带模型；Eq. (19) 存待核查差异；无 Wannier 生产链路 |
| 层电偶极×轨道磁矩的混合几何与电流 | [IMHE](papers/das-2026-layer-orbital-imhe/README.md) | 多带混合导数、A16完整乘积展开、有限外场独立核验 | 固定参考磁矩的有效模型；非R5G/SI/全磁输运；原文差异单列 |

**能力详情与API：** [能力目录](physics/README.md) · [机器可读清单](physics/capabilities.json)
**论文覆盖：** [王健组近期文献与缺口](docs/LITERATURE_COVERAGE.md)——五篇核心接入不代表近期PRL/PRB全部覆盖。

## 运行

在仓库根目录，使用Python 3.11或3.12：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-core.txt
python scripts/core.py list
python scripts/core.py run igmc --output outputs/igmc.json
python scripts/core.py test igmc
```

模块ID依次为 `nmee`、`spin`、`magnetotransport`、`curvature-flux`、`igmc`、`altermagnet-bpve`、`imhe`。

```bash
python scripts/core.py run all --output outputs/all.json
python scripts/core.py test all
```

`run` 真正执行各篇已核验例子，输出计算值、单位、边界与源代码提交；不会仅返回已有快照。全部运行需要数分钟。`check`只检查目录连接，物理验证使用`test`和`run`。输出不会覆盖原冻结数据。

要换自己的Hamiltonian、算符或参数，从[复用指南](docs/REUSE.md)进入各模块Python API。这个运行器不把不同物理协议强行变成同一个响应公式。

## 仓库布局

```text
physics/                  # 当前能力目录与机器清单
scripts/core.py            # 运行与测试入口
papers/                   # 七篇独立实现、正文/SM核查、冻结基准
studies/                  # 历史跨论文探索，不作为核心能力验收
integration/              # 旧接口设计与材料输入契约
outputs/                  # 用户本地运行结果（不入库）
docs/REUSE.md             # 新任务如何调用
docs/LITERATURE_COVERAGE.md# 已覆盖与未覆盖论文
```

[验证与未决问题](RELEASE_STATUS.md) · [逐篇交接](RESEARCH_HANDOFF.md) · [新增能力规范](CONTRIBUTING.md) · [GitHub Actions](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions)

## 维护原则

- 新任务先选观测量、驱动、频率/占据和单位，再选模型；不能从HR文件猜物理自旋或位置算符。
- 源文不一致、我们的误读、协议差异和缺失作者数据分别记录，不通过调系数追平图片。
- 旧数值基准和逐篇证据保留。首页仅维护当前入口；此前逐轮首页见[整理前版本](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/d65fc43/README.md)。
- Publisher PDF/SM不在此镜像。许可证按子目录现有声明处理，不将IGMC的BSD许可证自动扩展到全仓库。
