# IGMC 第四轮复现报告｜磁点群全覆盖（2026-10-19）

研究来源：Xiang et al., PRL 134, 116301 (2025) 主文 Eq.(14)、SM §V Tables I–III。外部群论参照：`Xin-Jiaqi/group-theory-operations-toolkit` 0.18.0，固定 commit `551a0daf...`。

## 核心成果

1. 读取用户 GitHub **真实的 122 个磁点群生成元和时间反演标签**，将最小输入版本化；重新用 9 维张量零空间求解所有群，逐群验证磁群闭包、矩阵正交性、张量基不变性，并对 244 个张量空间分别使用完整磁群特征标平均做第二条独立维数验证。
2. 结果：conduction 允许 58 个 MPG，分成**纵向 24、横向 18、两者均有 16**；displacement 同样允许 58 个 MPG，分成**纵向 23、横向 19、两者均有 16**；各有 64 个被完全禁止。SM 代表群的具体维数/结构与计算一致。
3. 源文 Eq.(14) 对两个通道一律使用时间反演因子 −1；按本工作已明确的实响应分通道协议，若机械用于 C，则 122 个 MPG 中 **58 个群的允许张量维数发生变化**，包括 `4mm1′` 从应为允许 1 变成 0。仍然保留源文定义为命名分支，不宣称自动纠错。
4. 澄清源码与 SM 的 **单斜轴选择不同**：toolkit 的 `2_z` / `m_z` 与论文唯一轴 y 通过坐标 y↔z 变换互相匹配，不应当误报为科学矛盾。

| 关卡 | 结论 |
|---|---|
| 122 MPG 群生成元闭包 | PASS |
| 244 个 C/D 张量空间及特征标独立维数检查 | PASS（Neumann/SVD 及群平均 244/244 一致） |
| SM Tables I–III 代表群及分类数量 | PASS（代表性 / 计数），完整符号表逐项核对 PARTIAL |
| 源 Eq.(14) 反幺正因子 | SOURCE/PROTOCOL DISCREPANCY；58 个 MPG 维数差异 |
| 外部 GitHub 原生 toolkit Python API 实机比对 | BLOCKED（离线未安装，显式跳过） |
| Fig. 2(d) 原图幅度 | 继承上一轮 PARTIAL，未修改物理因子 |

## 交付

- 源码：`src/igmc_lab/toolkit_bridge.py` 与 `data/toolkit_mpg_generators_v018.txt`
- 全群结果：`results/raw/round4_mpg_122.csv`，`round4_mpg_audit.json`、`round4_mpg_tensor_shapes.json`
- 可回归冻结：`results/frozen/round4_mpg_reference_v2.json`
- 详细原理与接口：`docs/ROUND4_TOOLKIT_SYMMETRY_AUDIT.md`

整体论文依旧 **PARTIAL**：群论张量计算已覆盖 122 磁点群，但图像绝对幅度、SM 个别前因子和原作者数据来源问题未消失。仅本地执行，未提交 GitHub；旧 NMEE 与前三轮数值基准保持不变。

### 全张量形状补充

本轮最后还新增 `results/raw/round4_mpg_tensor_shapes.json`，提供全部 244 个 C/D tensor space 的一组标准化独立张量基（不是仅计算维数），并独立对这些基再次执行磁群 Neumann invariance check。连同旧测试共计 **40 个发现，39 个通过，1 个期望跳过**（外部 toolkit live API）。

**额外独立验证：** 使用群平均特征标公式 `det(R)[Tr(R)]²`，重算所有 244 个允许维数，与 SVD 结果 **244/244 一致**；当前统计 **40 tests, 39 pass, 1 optional skip**。
