# Stage8l 零作业核验与伪自旋矩阵元门控 —— 中文小结

**日期**：2026-10-22　**状态**：PASS（全部四项 SHA 冻结值匹配；未新增任何 SCF 作业）

## 1. 任务书对应

`STAGE8M_DEEPSEEK_ROUGH_DFT_RESPONSE_PILOT.md` B 段：先完成 Stage8l 的零作业证据纠错（带清单 off-by-one）与占据带→空带伪自旋矩阵元筛选；通过后才允许 Stage8m 后处理与响应估计。

## 2. 带清单纠错（零作业）

- Stage8k 证据 JSON 中 `E_band35_eV` / `S_grouped` 键实际存放的是 **1-based 第 36 带** 的数据（`spinor_rows[2]` 是 0-based 列表第 3 行 = 带 36）。能量与自旋数值有效，仅标签差一。
- 7 条非简并记录的能量与自旋全部匹配 legacy 带 36（`matches36=True`）；glide 伪 K′ 点简并无法消歧（预期行为，不判错）。
- 产出：`STAGE8L_BAND_LEDGER_CORRECTED_20261022.json`（原始 Stage8k JSON 未改动）。

## 3. 占据/空带窗口核验

- 占据带 33–36（1-based）占据数 = 1.000；空带 37–44 占据数 = 0.000。
- 各采样 k 点直接带隙 1.15–2.99 eV（AB 与 glide 均在绝缘窗口内）。

## 4. 伪自旋跨隙门控（四个 k 点，AB 与 glide 各一）

方法：WAVECAR 自旋行 Euclidean 归一化的伪自旋（非 PAW 物理自旋、非速度、非 NMEE）；比较量 = 占据(33–36)→空(37–44) 块 Frobenius 范数（简并子空间内 pair 相位依赖，块范数为规范稳健量）。

| k 点 | AB 块 Frobenius (S_x, S_y) | glide 块 Frobenius | 直接带隙 (AB) |
|---:|---:|---:|---:|
| Γ–M 2/3（非 K′） | 0.206 | 0.203 | 2.93 eV |
| 近 M（非精确 M） | 0.147 | 0.149 | 2.95 eV |
| Γ–K 内点 | 0.137 | 0.141 | 2.20 eV |
| K（真 K） | 0.029 | 0.027 | 1.15 eV |

典型最强 pair：33→42、34→41、35→38（|S| 0.03–0.05，dE 2.9–3.5 eV）。

## 5. 科学判读与边界

1. **跨隙伪自旋翻转通道存在**：四个 k 点、两种堆垛的 Frobenius 范数均非零。
2. **AB 与 glide 顶点强度同量级**（0.206/0.203、0.147/0.149、0.137/0.141、0.029/0.027）：符合"带间 SOC 混合由单层自旋轨道强度主导、堆垛主要改变速度顶点与 BZ 求和"的预期。
3. **谷点 K 处弱**（~0.03）：Ising 型面内自旋压制的定性体现。
4. **必要非充分**：跨隙自旋翻转通道存在 ≠ Γ^stat 非零；还需速度顶点与全 BZ 求和（Stage8m D 段）。本小结不宣称新颖性或物理效应成立。

## 6. 产物

- `STAGE8L_SANITIZED_GATE_SUMMARY_20261022.json`（净化 JSON：SHA、窗口、块范数、判读；无原始路径/大型文件）
- `STAGE8L_BAND_LEDGER_CORRECTED_20261022.json`
- 门控脚本：`stage8l_pseudo_interband_gate.py`（SHA 校验 fail-closed）
