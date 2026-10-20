# IGMC 论文复现第五轮（2026-10-20）

## 实际完成

对照原论文 SM pp.8–10 的 Tables I/II，逐行人工建立 23 个原表张量模板（C 11 行、D 12 行），覆盖 C、D 各 58 个允许磁点群。再结合已有 122 磁点群的完整零空间结果，逐一比较 **9 个笛卡尔张量分量及其约束**，另外显式验证各 64 个禁戒群。

| 验收目标 | 结果 | 状态 |
|---|---:|---|
| Table I conduction 张量：全部 122 群（含禁戒） | 122/122 | PASS |
| Table II displacement 张量：全部 122 群（含禁戒） | 121/122 | PARTIAL，定位 1 项表格约束缺失 |
| 两张表合计 | 243/244 | 一处源表异常，不是数值误差 |
| Table III 允许/禁戒数量与纵/横向类型 | 沿用前轮 122 群与分组 | PASS 分类；未独立逐文字抄录 |
| 原 toolkit 运行时本地 Python API | 未安装 | BLOCKED / optional skip |
| 原文 Fig. 2(d) 幅度 | 未改变 | PARTIAL |

## 发现新的、可解析复核的问题

SM Table II p.9 把 `-42'm'` 和 `2'2'2`、`mm2` 列为同一一般张量形式，其中 \(\sigma_{xy}\)、\(\sigma_{yx}\) 两者独立；但对于 `-42'm'`，它们必须满足 \(\sigma_{xy}=\sigma_{yx}\)。使用群生成元 `S4` 与 `C2xT` 可逐项手算证明，数值零空间则给出唯一基

\[
\sigma^D_{\bar 4 2'm'} = a \begin{pmatrix}0&1&0\\1&0&0\\0&0&0\end{pmatrix}.
\]

该问题 **仅影响独立参数数目，不改变非零分量位置或 Table III 分类**；它与 Eq.(14) 的符号协议以及 SM Eq.(34) 的因子二问题是三个不同的审计项目。

## 复现证据与运行

- [详细审计](docs/ROUND5_SM_TABLE_SOURCE_AUDIT.md)
- [来源表转录](src/igmc_lab/paper_table_source.py)
- [逐群逐通道对照表](results/raw/round5_source_table_244.csv)
- [机器审计摘要](results/raw/round5_source_table_audit.json)
- [冻结审计](results/frozen/round5_source_table_reference_v1.json)

从论文目录运行：`PYTHONPATH=src python scripts/verify_release.py`。原 PDF 不进入可公开发布的 ZIP 或 GitHub；未修改任何旧冻结计算结果，未推送 GitHub。

**结论：** 本轮群论 source-level 审计从“代表群 + 总数”推进至完整 244 组张量形式；整篇仍为 `PARTIAL`，Fig.2(d) 不得声称完全复现。
