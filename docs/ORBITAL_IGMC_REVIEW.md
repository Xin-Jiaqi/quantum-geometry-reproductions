# Quantum geometry 仓库 review 与轨道 IGMC 接入

日期：2026-11-09。审查基线：`9ffa75d3d48d8ff9c2152a43eebfb3cbf705ff5a`；集成时同时保留 main `c088f4202423bfb20a1077c127236cc93fe6f36d` 新增的 IMHE，最终登记 8 个模块。

本次核查能力登记、贡献规范、统一运行器/CI，以及与新论文直接相关的 CuMnAs Hamiltonian、自旋 IGMC 内核和多带轨道几何实现。未对所有历史论文重新逐式审稿；并行新增的 IMHE 仅合并保留，未纳入本次代码审计。

| 既有模块 | 与本次关系 | 本次处理 |
|---|---|---|
| NMEE / Jia 2026 | 已有 CuMnAs、PT 简并处理；默认 Jn=0.6 eV | 加跨模块 Hamiltonian/速度对照；新论文用 Jn=0.85 eV |
| Spin / Xiang 2025 | E/E² 驱动自旋响应 | 保持原接口，不把 E 响应当作 B 响应 |
| Magnetotransport / Xiang 2026 | 已有非简并多带几何、轨道 PHE/MNHE；一般位置修正仍有边界 | 借鉴固定轨道基/失效即拒绝规范，不更改原能力宣称 |
| Curvature–Flux / Cui 2026 | Zeeman 几何、带间磁电响应 | 保留原频率/算符协议，不合并公式 |
| Spin IGMC / Xiang 2025 | 最接近的新论文前序：两带自旋磁位移电流 | 增加 Q01 与新自旋核归一化交叉测试；不覆盖原幅度争议 |
| Ezawa BPVE 2025 | 已有量子度规磁性光伏 | 与本次 B 线性响应分开，旧冻结数据保持原样 |
| IMHE / Das 2026 | 主分支并行新增的层电偶极/轨道磁矩混合响应 | 原代码、数值与能力声明原样保留；未重新验收 |
| 新增 Ghorai 2026 orbital IGMC | 轨道磁位移电流及 CuMnAs 反铁磁读出 | 独立目录、两公式协议、统一运行入口和双 Python 版本 CI |

## Review 发现

1. 审查基线的根目录测试仍断言仅 5 个 ID，而当时能力登记已经有 6 个。最终与 IMHE 的并行更新合并，入口测试覆盖全部 8 个 ID，不通过跳过测试规避。
2. 旧 NMEE 与本文 Hamiltonian 形式一致，但参数和自旋磁矩约定不同，不能直接复用默认参数。跨论文测试显式映射参数。
3. 本文正式版 D10 与直接 A9 求解的轨道响应不等价。该差异与简并随机旋转、cosine/phasor 整体因子无关；保留两入口和可执行差异记录，详见 [公式审计](../papers/ghorai-2026-orbital-igmc/SOURCE_AUDIT.md)。

## 新能力覆盖

`orbital-igmc` 覆盖指定固定轨道基模型的轨道/自旋磁极化及低频位移流；带隙平台、Néel 反号、镜面选择规则和频率/单位均有对照。D10 保留三项，A9 有矩阵函数/有限差分投影独立验证。

不扩成材料高通量、完整共振/无序输运、电磁电路或全图复现；不上传出版社 PDF、不改其他论文冻结结果。`MODEL_READY` 只适用于明确选定的协议，唯一物理幅度保持 `OPEN_SOURCE_DISCREPANCY`。

运行：`python scripts/core.py test orbital-igmc`；`python scripts/core.py run orbital-igmc --output outputs/orbital-igmc.json`。
