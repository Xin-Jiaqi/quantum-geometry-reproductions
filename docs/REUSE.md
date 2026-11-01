# 在新任务中复用

先用根目录例子证明环境和目标通道能运行，再替换模型。每个研究任务固定代码提交，并保留参数、物理约定和收敛结果；根运行器的JSON已包含提交与工作区是否改动。

## 选择接口

1. 磁化对E×B：NMEE的`CovariantBloch2D`与`DynamicSpec`，参见[实例](../papers/jia-2026-nonlin-magnetoelectric-edelstein/examples/use_dynamic.py)。精确简并可由完整动态方法处理。
2. 自旋对E或E²：[SpinModel例子](../papers/xiang-2025-spin-quantum-geometry/scripts/reproduce_core.py)。Hamiltonian、连接和自旋均是可替换函数；current_axis选择自旋流，未指定时为自旋极化。
3. 自旋/轨道磁输运：[QGMap例子](../papers/xiang-2026-quantum-geometric-map/scripts/reproduce_core.py)。先选通道，再复用其中几何与响应内核；不要把三带moment子项称为完整轨道响应。
4. 几何/拓扑/带间互易响应：[Curvature–Flux例子](../papers/cui-2026-curvature-flux-zeeman-qgt/scripts/reproduce_core.py)。输入H、两方向D_kH和三个同规范物理Pauli算符；源频率与驱动频率相反。
5. 磁场驱动电流：[IGMC例子](../papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/scripts/reproduce_core.py)。当前可替换两带Hamiltonian和速度，但物理自旋固定为Pauli；默认接口不是任意Wannier自旋模型。

## 每个新算例保留六项信息

- H(k)及其坐标/能量单位、BZ积分测度。
- 测量算符与外场耦合算符，使用同一轨道规范；说明位置连接。
- 复频率、cosine/Fourier幅度、输出频率以及静态或动态协议。
- μ、温度、简并策略；区分绝热展宽与物理散射。
- 至少一个适用的独立对照，以及目标参数附近的网格/微分步长收敛。
- 原始响应结果、SI换算因子与局限，分开保存。

改模型后不必重做作者全部图像；需要重新验证的是该模型改变的假设。不要从同一个公式的两种写法获得“独立验证”。

## 结果与历史

本地结果写入`outputs/`。要成为可复用基准时，按[新增规范](../CONTRIBUTING.md)加入相应模块。`papers/*/results`中的历史冻结值是回归参照，不能重算覆盖后再宣称通过。

当前根命令执行完整示例而非读取快照；每次只导入对应论文目录。脚本输出中的`MODEL_READY`是目录状态，运行失败仍会报错退出。统一数据对象的设计稿保留在`integration/`，未实现的适配器不视为物理能力。
