# 2026-11-08 正文材料图更新

本轮仅重组既有结果，没有新科学计算。正文图7替换早期有限准备场的投影控制，使用真正零磁场自洽参考的96/112态完整材料响应；图8补能谱—电子态与跨层混合干预。

运行 `python plot_material_figures.py` 重建两图，输出到本目录的 figures。PNG阅读副本位于正文附件；PDF/SVG及逐点CSV在此存档。需要 numpy、matplotlib。

数据来源（相对 studies/bsf-nmee-connection）：
- docs/revision-20261106c/figure_data.json：图7b/c历史独立有限场材料对照。
- results/mechanism/zero-field-20261107/refined/controls96-analysis.json、full-controls-analysis.json：图7a及图8a。
- results/mechanism/mixing-channels-20261107/material-channels.json：图8b。
- 正文图9/10沿用 results/materials/mos2-interface-20261108 与 mos2-compression-20261108 的既有图源。

六种设置展示方法敏感性，不是误差棒。各通道是修改Hamiltonian后的完整响应，不可相加。α=0保留界面其余作用；不等于自由单层或纯静电极限。两种材料支持电子态重组的重要性，未建立唯一通道、统一阈值或实验读出方案。
