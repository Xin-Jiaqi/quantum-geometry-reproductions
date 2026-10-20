# 第七轮｜金属费米面收敛与整篇论文验收矩阵

## 核心完成内容

- 对 Eq.(17) 两带 Dirac 模型，严格用 `E(k)>|µ|` 确定零温占据域。把 `E²=µ²` 在 `cos(akx)` 中写为二次方程，精确确定每个 ky 切片的全部费米面交点，并按所有 Lifshitz/端点事件切分 ky 积分区间。避免原先的规则 k 网格通过偶然命中或遗漏费米面造成的尖峰漂移。
- **独立两路：** 直接积分主文 Eq.(11) Zeeman 曲率的占据域；由主文 Eq.(13) 对费米线做 delta-function 积分。两条数值积分路径不共用 C 的被积核。
- 对 `µ=0, ±0.8, ±0.95, ±1, ±1.05, ±1.1, ±1.2 eV` 收敛，包含 1 THz 完整有限频率分母（非共振）。对 `160→224` 外部 Gauss 精度的最大残差：`C=1.0623e−9 Å−1`，`D/ω=3.4350e−10 Å−1 eV−1`；占据域 C 与独立费米线 Eq.(13) 最大误差 `5.3496e−10 Å−1`。

## 关键可复现数值

模型参数：`β=.25 eV`、`v_F=.825 eV Å`、`M=.1 eV Å²`、`a=1 Å`、`g=2`、`T=0`、`f=1 THz`（循环频率，`ℏω=.004135667694 eV`）。纵轴使用原文约定 `−gµ_B σ`。

| 化学势 | C_xy, A/(m T) | D_xx, A/(m T) |
|---|---:|---:|
| 0 eV（绝缘体） | 约 0 | +0.0177539 |
| ±1 eV（金属） | **−33.48623534** | **−0.06417844743** |

在离散 μ 扫描点中，conduction 峰位约 `|µ|=1.05 eV`，displacement 峰位约 `|µ|=1.01 eV`；并非精确优化峰位。费米面模型对 `µ↔−µ` 有精确对称性。

**仍未解决**：Fig.2(d) 图中纵轴视觉幅度约大一倍；不能用 SM Eq.(34) 的局部 ×2 中间式取代多种独立方法验证的正文 Eq.(12)。先前发现的 Eq.(14) 时间反演因子与 SM Table II `−42′m′` 参数关系疑点继续保留。

## 完整研究状态

详见 `docs/FINAL_ACCEPTANCE_MATRIX.md`，按每张图、每个 SM 小节、群论和软件发布的独立验收要求标明 `PASS` / `PARTIAL` / `BLOCKED`。本轮完成金属费米面计算方法的可信收敛，并未获得原始作图数据，因此**全篇仍为 `PARTIAL`**。

## 本轮交付

- `src/igmc_lab/fermi_contour.py`、`scripts/audit_round7.py`、`tests/test_round7_metal_fermi.py`。
- `results/raw/round7_metal_convergence.csv`、`round7_metal_curve.csv`、`round7_metal_audit.json`。
- `results/figures/round7_fig2b_metal_Fermi_boundary.png`、`round7_fig2d_metal_Fermi_boundary.png`、`round7_metal_convergence.png`。
- `docs/ROUND7_METALLIC_FERMI_CONTOUR.md` 和 `docs/FINAL_ACCEPTANCE_MATRIX.md`；新增冻结验收数据，往轮数据不改。

GitHub 尚未上传，未发送作者邮件，也未发表/镜像论文或 SM 原 PDF。
