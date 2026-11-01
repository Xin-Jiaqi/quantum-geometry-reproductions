# 核心能力与未决问题

更新：2026-11-01。五个模块在声明的模型/协议下为`MODEL_READY`；不能从CI通过推导出全文、真实材料或所有输运极限已复现。

**新增可视化验收：** [五篇正文主图图册](docs/MAIN_TEXT_FIGURES.md)包含 8 张正文对应图和 2 张额外检验图，附 PNG/SVG/NPZ、逐篇参数与重算命令。
图册调用已有内核；NMEE Fig. 1(c) 保留为静态对照，不阻塞核心能力接入，也不声称已修正原文 DC 曲线。

## 已验证的核心

| 模块 | 回归规模 | 独立证据与例子 |
|---|---:|---|
| NMEE | 38 | 轨道Liouville、有限k导数、PT/规范/六带；481²→721²相对变化5.83e-7 |
| Spin | 76 | SM四源项、轨道求解、物理自旋流；非零四带网格与步长检查 |
| Magnetotransport | 1214 | spin四带、orbital有限B对照；Eq56链式法则更正；保留历史数值 |
| Curvature–Flux | 375 | 本征态/投影子与Kubo/Liouville；四带64²→96²相对变化6.86e-8 |
| IGMC | 76 | 波包/格林函数/Liouville、金属费米海/费米面、SI与频率桥 |

IGMC在无外部群论工具的干净环境中跳过1项可选实时交叉检查；122群离线检查仍运行。本地显式载入完整toolkit时76项均通过。测试数量是回归规模，不是独立物理证明的数量。

当前提交的远端结果以[Reusable physics cores](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/workflows/core-capabilities.yml)及各篇旧冻结工作流为准。根运行器和CI实际执行例子并上传JSON，不以文件存在替代数值运行。

## 争议与边界

| 论文 | 已厘清 | 尚未闭合 |
|---|---|---|
| NMEE | 完整有限频率E/B通道；SM9/11局部对角项符号冲突单列 | 作者Fig1c本征DC、物理散射极限 |
| Spin | SM26/45明确省略项；完整动态内核保留多带过程 | 印刷积分/图像幅度；通用金属τ⁰提取 |
| Magnetotransport | **撤销旧Eq56错误判断**：∂ab f含f''链式项 | SM20/31符号、49/52因子、Fig2 SI来源；一般多带轨道位置项 |
| Curvature–Flux | **Eq16与112一致**；Eq19与126的一阶项冲突，独立动态支持126；无效m′/C3反例已排除 | Table I全表、完整耗散金属KME |
| IGMC | SM34/46局部系数问题；完整频率/低频差异不应归因于粗网格 | 作者Fig2(d)幅度来源、Eq14整体反幺正约定 |

不通过经验倍率追平作者图片。待核问题不妨碍已验证的指定模型使用，也不能被隐藏为“全篇PASS”。

[当前调用目录](physics/README.md) · [逐篇详细证据](RESEARCH_HANDOFF.md) · [文献覆盖](docs/LITERATURE_COVERAGE.md) · [整理前完整历史状态](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/d65fc43/RELEASE_STATUS.md)
