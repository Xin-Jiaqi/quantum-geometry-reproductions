# 第八轮 | SM 半经典与 Green 函数路径闭合

## 新增关键结果

1. 按照 SM Eqs. (21),(24) 直接计算两频率谐波的 wavepacket 位置修正，与从 SM Eq. (44) 原始矩阵 resolvent 构造的 lesser Green's current 在两带 Dirac 模型的完整复频率下逐点一致。覆盖 6 个张量分量、101² 网格；最大绝对差小于 `1.4e−13`（模型核单位）。
2. 用新的 GF 方法得到 `µ=0, f=1 THz` 的 `Im K_xx = 6.30109046490e−5 Å⁻¹`，`−gµB σ_D_xx=0.01775610 A/(m·T)`（101² 网格）。与第六轮收敛数值 `0.01775389` 吻合，残差是网格精度而非 ×2。
3. 来源公式审计：SM Eq. (34) 把 `Q=Re(rσ)` 写成 `2Re(rσ)`；本轮证实是 **wavepacket 展开中的局部重复倍乘**。主文 Eq. (10) 的 `F` 正弦极化项与 SM Eqs. (21),(24) 相比，在同一 `Z` 定义下还有符号差异，但对应的电流效应是 `O(ω²)`，不改变本文保留的 `O(ω)` displacement。
4. SM Eq. (46) 的 `Bcos(ω1t)` 傅里叶变换系数印成每个谐波 `2πB`，按其变换约定应为 `πB`；Eq. (47) 又使用了正确的 `B/2`。独立 FFT 校核支持正确谐波振幅 `0.5B`。

**结论**：SM §III 和 §IV 在本文两带、无共振的模型响应层面已构成自洽的、可以测试的独立路线。Eq. (34)、Eq. (46) 的局部错误尚无证据证明会影响作者原始绘图；Fig. 2(d) 原始绝对幅度继续保持 `BLOCKED`。不得将本次定位的倍数差当作匹配原图的许可。

## 代码与复现

- `src/igmc_lab/green_wavepacket_audit.py`、`scripts/audit_round8.py`、`tests/test_round8_green_wavepacket.py`、`parameters_round8.json`；
- `docs/ROUND8_SEMICLASSICAL_GREEN_CLOSURE.md`（含逐式推导）；
- `results/raw/round8_GF_wavepacket_audit.json`、`results/figures/round8_GF_wavepacket_audit.png`。

后续重点不是重复画曲线，而是固定此 validated-method subset，完善最后的 Eq. (14) 复频率响应定义并取得 Fig. 2(d) 原始参数或数据。没有联系作者或修改远端仓库。
