# Quantum Geometry Reproductions

A paper-by-paper collection of **independent, reproducible quantum-geometry calculations**.

本仓库面向量子几何及相关自旋、磁电输运响应的理论与数值复现。核心原则是：**一篇论文一个独立模块；一项科学结论对应明确的计算协议、证据和验收状态。** 新论文可以复用经过验证的工具，但不覆盖已有论文的冻结结果。

> **Scientific status ≠ test status.** 自动化测试通过，表示指定实现与冻结基准一致；不代表论文所有公式、图像或真实材料预测均已严格复现。本仓库不是论文作者的官方实现。

## Paper index

| Paper | Physics / reproduced scope | Scientific status | Validation |
| --- | --- | --- | --- |
| [Jia et al., *PRL* **136**, 226302 (2026)](papers/jia-2026-nonlin-magnetoelectric-edelstein/) | Nonlinear magnetoelectric Edelstein effect (NMEE); Kane–Mele / CuMnAs models; explicit published-vs-static response protocols | **PARTIAL** · driven-DC / static-protocol discrepancy remains open | [Paper README](papers/jia-2026-nonlin-magnetoelectric-edelstein/README.md) · [CI](.github/workflows/nmee.yml) |
| [Xiang et al., *PRL* **134**, 116301 (2025)](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/) | Intrinsic gyrotropic magnetic current (IGMC); Rashba / Dirac models, Zeeman quantum geometry, numerical and symmetry audits | **PARTIAL** · published Fig. 2(d) magnitude and some source conventions remain open | [Paper README](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/README.md) · [CI](.github/workflows/igmc.yml) |

每一行的状态指向**整篇论文的科学复现**，而不是局部模型计算的测试结果。详细的公式约定、数值误差、已解决与未解决问题均以各论文目录的报告为准。

## Repository architecture

```text
quantum-geometry-reproductions/
├── README.md
├── papers/                           # 独立、可回归的论文复现
│   ├── jia-2026-nonlin-magnetoelectric-edelstein/
│   └── xiang-2025-intrinsic-gyrotropic-magnetic-current/
├── studies/                          # 基于已有模块的跨论文探索
│   └── bsf-nmee-connection/
└── .github/
    └── workflows/
        ├── nmee.yml
        ├── igmc.yml
        └── bsf-nmee-connection.yml
```

**`papers/`：独立复现。** 每篇论文拥有自己的依赖、计算程序、测试、冻结参考值、物理协议与验收报告。不同论文的内部 Python 包布局可以不同，不强制抽象出尚未经验证的“通用响应核”。

**`studies/`：后续研究。** 跨论文连接、新模型或物理假设放在这里，通过明确的接口调用已有论文模块，不反向修改原论文的冻结物理基准。例如 [BSF × NMEE connection](studies/bsf-nmee-connection/) 研究的是合成模型中的对称性与层滑移响应关系，并非已验证的真实材料铁电开关预测。

**`.github/workflows/`：独立验收。** 每个模块使用自己的 GitHub Actions 工作流与路径过滤。新增论文原则上不改变现有论文的测试门槛或 CI 定义。

## Quick start

需要 Python **3.10+**。以下命令在每个论文目录内单独运行；建议分别使用虚拟环境，避免不同论文的依赖相互干扰。

**NMEE · Jia et al. (2026)**

```bash
cd papers/jia-2026-nonlin-magnetoelectric-edelstein
python -m pip install -e .
python scripts/verify_release.py
```

**IGMC · Xiang et al. (2025)**

```bash
cd papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
python -m pip install -e .
python scripts/verify_round9.py
```

上述命令验证的是各自**明确界定的模型计算与历史冻结数据**。绘图、参数扫描和额外收敛测试参见对应论文的 [NMEE 使用说明](papers/jia-2026-nonlin-magnetoelectric-edelstein/README.md) 与 [IGMC 使用说明](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/README.md)。

## Adding a new paper

每篇新论文创建一个并列目录：

```text
papers/<first-author>-<year>-<short-title>/
├── README.md                   # 文献来源、复现范围、运行方式、科学状态
├── pyproject.toml              # 可独立安装的依赖
├── src/ or <package>/          # 论文专用计算实现
├── scripts/                    # 运行与验收入口
├── tests/                      # 单元测试及物理回归
├── results/                    # 允许公开的独立生成结果、冻结基准
└── docs/                       # 公式协议、来源、收敛性与开放问题
```

这是**推荐的最小组织约定**，不是要求将已有论文重构成统一内部包结构。接入新论文时：

1. **冻结复现范围**：记录 DOI / SM、模型参数、符号与单位、频率及温度约定、论文图表与公式的对应关系。
2. **建立科学验收表**：区分解析验证、数值收敛、原图定量一致性、真实材料推广及尚待作者确认的问题。
3. **独立实现与交叉检查**：保留输入和来源哈希、关键数值的冻结基准，禁止通过任意比例因子强行拟合原图。
4. **加入模块级测试与 CI**：提供单条离线验收命令，并新建按该论文目录触发的 `.github/workflows/<paper>.yml`。
5. **更新本页索引**：添加论文链接、复现范围、`PASS / PARTIAL / BLOCKED` 科学状态和独立 CI；不修改其他论文的冻结基准。

新论文的具体问题与运行细节写在**自己的 README**，根 README 只提供导航和共同规则。跨论文研究则新增到 `studies/<topic>/`。

## Scientific reproducibility rules

- **明确区分协议。** 平衡态静态响应、外场驱动的 DC 响应以及 AC 有限频率响应不能仅因符号相似就混为一个系数。
- **冻结且可审计。** 记录模型、占据、展宽、(k) 网格、物理单位、积分方法、输入哈希及数值误差。修正理论问题时增加新基准或说明迁移原因，不静默覆盖历史结果。
- **科学状态独立于 CI。** 建议使用 `PASS`（声明范围内的科学验收已完成）、`PARTIAL`（有可靠子集但仍存在开放问题）、`BLOCKED`（缺少关键输入或尚无法开展必要验证），并在文档里注明状态作用的范围。
- **模型复现不等于第一性原理预测。** 仅有 `wannier90_hr.dat` 与 Wannier 中心不足以保证定量正确的 NMEE / IGMC；相关真实材料计算还需要**同一规范下**经过验证的自旋及位置 / Berry 联络算符。
- **来源与版权清晰。** 公开独立生成的代码、图和适合发布的基准；不直接提交未经授权的出版商 PDF、原始图像、私人输入、受限制软件或内部研究归档。

## Connection to other workflows

已有的 [Wannier workflow](https://github.com/Xin-Jiaqi/wannier-workflow) 与 [NLO response workflow](https://github.com/Xin-Jiaqi/nlo-response-workflow) 采用冻结 manifest、SHA256、占据带数与 top-K 模型交接。本仓库可以复用其**只读输入审计和工作流约定**，但目前不声称已经具备通用真实材料 NMEE / IGMC 计算能力。实际算符完备性与规范一致性必须分别通过物理验证。

磁点群计算可参考 [Group Theory Operations Toolkit](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit)；借用的数据与版本须在对应论文模块中注明来源、固定版本并保留许可证归属。

## License and attribution

这是独立研究复现仓库，**不是原作者发布的软件或数据**。目前仓库根目录**没有统一的开源许可证**；各模块的许可范围应分别核对：

- [IGMC](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/LICENSE)：独立实现采用 **BSD-3-Clause**，第三方群论数据保留其自身署名与许可声明。
- [NMEE](papers/jia-2026-nonlin-magnetoelectric-edelstein/)：当前没有单独声明的开源许可证，不应从仓库公开可访问推断出默认再使用许可。

具体引用、原始文献 DOI、代码来源及未解决科学问题，请优先查阅相应论文目录。
