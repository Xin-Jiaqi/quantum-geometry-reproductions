# 2026-11-06：恢复中间 KM 堆垛，核清补偿的矢量含义

## 用户授权与本轮范围

先解决为什么主图用 KM AB–BA 路径四分之一构型，以及当前主文是否物理自洽；之后才在对话讨论同行工作与项目整体架构。本轮没有提交 DFT、没有登录日本资源，也没有启动后续电子机制研究。全部新计算为本地既有八带模型的有限核查。用户原始批注稿保持不变。

术语：KM AB/BA 是模型端点；MoSSe S/S 堆垛 A/B 是材料镜像伙伴。KM λ=1/4 的 m001 伙伴是 λ=3/4。

## 决定及关键新证据

1. 恢复中间 KM 构型为主图；同一固定结构呈现局域反向与层间补偿，比当前 KM AB 更适合直接对照。原 KM AB/BA 数据作为图 S4，不删掉门控未重现纯补偿这一结果。
2. λ=0.24、0.26，及固定 λ=0.25 时 m=44.9、45.1 meV，均支持两类过程。它是有限参数邻域中的存在性结果，不是所有滑移堆垛的普遍规律。λ=0.30 或 m=46 meV 的同一杂化范围未得到下层反向，负面结果也保留。
3. 发现并纠正旧称谓：原生 X 电场下旧约 10.25 meV 是 C_XX 总投影零点，C_YX 非零。不能称整个自旋相消。旧数值不改写，历史 worklog 加醒目更正。
4. 保留的竖直镜面 m_v 的面内法线为 native XY 的 (1/2, sqrt(3)/2)。主图沿法线施加 E，B 沿 z；镜面逐层禁止横向与面外自旋。新方向重新计算，未重标旧曲线。主文 x、y 是 native XY 旋转 +60 度后的右手系；MoSSe 原材料坐标不变。
5. 新主图：m=45、λ=.25，t=10→20 meV 时下层 −3.146→+1.077（单位 1e−5/(nm V T)）；t=10 固定、U=0→5 meV 时两层保持下负上正，总系数在 U=4.5–4.75 meV 间过零，线性插值 4.7266，仅报约4.7。更大 U 时单层也会反向，完整扫描在 S6b。
6. U=4.75 meV 独立读出三个自旋分量，横向与面外层/总系数约 1e−16，确认主图几何下的选择定则。完整原生面内张量转动后的禁止分量最大约1.1e−15。H 镜像协变残差约机器精度。

## 正文叙述边界

模型能证明两种过程在同一堆垛中可共存，也给出不依赖精确1/4位置的相邻证据。模型不含机械总能量，不证明这一固定中间构型是自由双层的基态，也不证明可滑移路径或能垒。材料中的稳定候选承担实际材料算例，但不等同于这个模型分数坐标。当前主线是“界面可反转指定层自身的带间磁电响应；总自旋变号也可能只来自补偿”。其可预测的能级/矩阵元条件仍待机制研究，不宣称 PRL 首创性已建立。

## 数据及计算源

所有 C 的原始单位为 (nm V T)^−1，质量、t、U 为 meV。完整保留以下探索过程：

| 数据 | 物理用途 | 计算源 |
|---|---|---|
| normal_main.csv、gate_extra_main.csv | 新主图：法线探测，较细积分 | km_normal_check.py main；km_gate_extra.py main |
| neighbourhood_robust.csv、gate_extra_robust.csv | 相邻 λ 和 m，统一检查精度 | km_neighbourhood.py robust；km_gate_extra.py robust |
| normal_robust.csv | 更远参数的适用边界 | km_normal_check.py robust |
| tensor_check.json | 旧原生XY完整面内张量，独立静态核 | km_tensor_check.py |
| normal_out_of_plane_robust.csv | 法线电场下的面外自旋 | km_normal_z.py robust |
| principal_main.csv、principal_robust.csv、out_of_plane_robust.csv | 沿镜面迹线的探索，面外自旋允许；未用作完整补偿图 | km_principal_check.py main/robust；km_out_of_plane.py robust |
| registry_robustness.csv | 旧X投影的相邻路径及少量离路径检查 | km_registry_check.py |
| vertical_mirror.json | 原始零场H镜像核查记录 | 本轮独立矩阵检查；由下列脚本再核验 |
| forbidden_at_compensation.json | 新过零邻域的纵向、横向、面外自旋及H镜像复核 | check_mirror_and_forbidden.py |

计算源保存在 calculation_sources，物理内核复用 study 中的 run_fixed_stack_layer_cxx.response，不改写基础模块。目录位于本仓库时自动定位依赖；在仓库外设 NMEE_REPO_ROOT 为仓库根。NMEE_AUDIT_OUTPUT 可指定输出目录，避免覆盖归档值。除张量/registry/选择定则脚本外，其余带 main 或 robust 位置参数，按上表运行。main 使用81×80×120谷区积分参数、0.03125 T导数；robust 使用41×45×72、0.0625 T。源脚本保留本轮具体点集，不包装成通用材料工具。

图由自足的 plot_revision.py + figure_data.json 生成 PDF/SVG/PNG；build_revision.py 说明怎样从上一版素材和新CSV装配本版。选择规则、缩放和展示区间另存 FIGURE_CONTRACT.md 与 figure-selection.json。MoSSe 图2、3及S1–3、S5的数据沿用第二轮。

## 核查与交付

- check_revision.py 通过：主图21点、邻域37点，所有点共用绝缘占据；主数据最小共同能隙约5.20 meV，层和与总值偏差小于1.6e−15。
- 新主图、KM AB/BA 对照与方向修正图完成视觉核查；全部图保留独立caption，PDF/SVG文字可编辑，PNG供讨论预览。
- 图源预检12项通过、0失败。两个提示为PNG仅320dpi且未导TIFF；本阶段是讨论稿，以矢量PDF/SVG为高质量输出，因此接受，不声称已完成正式投稿文件检查。
- 对照旧图、校验旋转顶角时始终保持全局Zeeman顶角为总Sigma_z，未把测量层投影到施加磁场一端。
- Markdown链接、美元公式边界和表格列数检查通过。原始批注稿及第二轮均保留。新读者入口仍使用稳定文件名“正文讨论稿-修订版.md”。

## 给后续接手者

不要把旧±10.25 meV投影零点与新约±4.7 meV法线探测零点混为一谈；新负根由已验证的水平镜像关系给出，本轮没有重跑完整器件环境扫描。以前关于公共标量偏置平移、镜像奇偶性的数学关系仍适用于相应标量响应，但旧数字不能直接作为完整自旋相消的证据。

下一步先与用户讨论：如何解释反向修正、两类过程共存的电子条件，以及相对 LEE/带边混合等工作的物理增量。用户未在本轮要求立刻展开新一轮机制计算或扩材料。当前已完成的是选点与观察量的澄清及阶段稿件更新。
