# 第九轮：validated-method subset 发布候选版

本轮不调整任何物理响应的归一化或数值基准，工作集中在发布边界、可重复执行与 fail-closed 验收。

- 新增 `release/release_contract.json`，机器可读标注 `PARTIAL`、Fig.2(d) `BLOCKED`、Eq.(14) `PARTIAL`、无经验放大因子。
- 新增 `scripts/verify_round9.py`，先审核发布内容和物理状态，再调用 Round1–8 冻结验证与全部单元测试。
- 新增 `tests/test_round9_release.py`，执行元数据、外部许可归属、证据文件存在性等 5 项回归。
- 准备 GitHub workflow `.github/workflows/igmc.yml`，只针对新 Xiang 目录，既有 NMEE 的工作流和冻结数据不会被修改。
- 真正 GitHub push 与远端 Actions 尚未执行；不构成远端 CI PASS。

见 `docs/RELEASE_CANDIDATE_ROUND9.md` 和本轮后附验证日志。

## 验证记录

- 本地冻结验证：`PASS_LOCAL_DRAFT`；单元测试 `Ran 76 tests`，`OK (skipped=1)`，即 75 通过、1 因 native toolkit 未安装跳过。
- 分发 ZIP 重新解压后，运行相同冻结验证仍通过。
- 离线全新 venv 初次安装由于环境没有 `setuptools.build_meta` 失败，并非 IGMC 代码或测试失败；随后明确使用本机预装的 `/opt/pyvenv/lib/python3.13/site-packages` 作为**外部依赖来源**执行 `pip install --no-deps --no-build-isolation -e .`，安装成功、包可导入、冻结验证通过。该测试是干净源码检验但不是“从完全空白依赖环境在线安装”的测试；CI 用 Python 3.11 常规安装流程，尚未在 GitHub 上运行。
- 保护上一轮 26 个此前带 SHA 的文件未发生改动；新增 R9 receipt 对新发布文件与第八轮核心冻结参考追加保护。
