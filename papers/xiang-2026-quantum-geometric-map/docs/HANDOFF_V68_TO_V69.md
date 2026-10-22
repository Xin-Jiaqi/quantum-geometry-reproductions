# v68 → v69 完整科学交接｜2026-10-22

仓库 `Xin-Jiaqi/quantum-geometry-reproductions`。Xiang PRB 原论文状态：**PARTIAL**。

## v68 新成果

首度在 v66 **完全相同投影器件极化电流算符** `J(E,B)=i[H(B)-E Pb(B),Pa(B)]` 下，换用 v67 的**真实半无限两电极**，计算完整有限频率、非平衡线性响应的 lesser Keldysh 双能量积分，不使用标量 gamma。均匀标量长度规范电场 **仅作用于有限器件**，两电极维持平衡；此协议既不同于 v66 现象学浴，也不同于 v67 的两端 reservoir-voltage Landauer 电导。不把两种浴的同算符结果数值不等视为误差，也不宣称源论文体相响应已复现。

Keldysh:
`delta G^<(eps+omega,eps)=-Gr_+ Pb G^<_- - G^<_+ Pb Ga_-`; `rhoE=-i Integral delta G^</(2pi)`; `sigma_device = Tr(i[H,Pa]rhoE)+Tr(-i[Pb,Pa]rho0)`；磁场导数包含 H、Pa、Pb、ΣL/R、ΓL/R 和二能量谱函数全部磁接触。

**Two-energy microscopic Ward**：`omega rhoE-[H,rhoE]+[Pb,rho0] = C_leads`，右侧由 reservoir Dyson lesser 自能**独立构建**，不是定义成左侧残差。ω=.29 三种规范的矩阵残差 ~9e-11；不连电极 site3 的局域 Ward ~1e-11，全局 charge accumulation 等于显式 lead charge injection。精确引入非对角投影位置导致电流电场接触 `-i[Pb,Pa]`。

模型默认 `(mu,T,B,q/hbar,omega,tlead)=(.13,.37,.23,-.84,.29,1.1)`，四轨道 Peierls 器件，三种规范。收敛能量积分 **640 Gauss 节点/子区间**（在 eps,eps+omega 的带边分段），64 节点失败为明确负控制。冻结数值 `sigma_device=-0.0052926496159-0.0011942573991i`，`dB_sigma_device=+0.0070700539001-0.0041462908484i`，独立 5 点全 B 差分验证。谱归一约 2e-9。静态 omega=0：`sigma≈+0.002996376699`，`dB_sigma≈+0.001901350531`，**重新计算非线性 finite E 器件 Green** 后 E/B 四角差分三种规范通过。

**VERIFIED**：所定义有限非相互作用双端 Keldysh AC、完整 2-energy Dyson Ward、数值算符及磁接触完备。**CONDITIONAL**：电场仅设备、外部给定轨道偶极、半无限单通道理想弹性电极、有限积分与范围。**UNRESOLVED**：完整体相均匀场、真正可观测宏观轨道流、实材料磁 Wannier 全嵌入及原 Xiang 2026 PRB 三带/sign/factor 辩论；原论文 **PARTIAL**。

代码：`src/qgmap/lead_keldysh_eb_v68.py`，40 tests `tests/test_lead_keldysh_eb_v68.py`，冻结 `scripts/reproduce_v68.py` + `results/frozen/v68_same_device_keldysh_eb.json`，独立 1 test `tests/test_lead_keldysh_eb_v68_frozen.py`，科学报告 `docs/FINITE_FREQUENCY_SAME_DEVICE_KELDYSH_EB_V68.md`。

v69 应补充一个**真正不同的、空间有限完整 Hamiltonian 的 lead-oracle**（例如足够长 lead + 有限接触能量的谱 Lehmann，并控制有限长度、能级展宽与时间分辨），及对 bulk uniform electric field/physical orbital current 的可操作同观测量匹配。不要用重新积分同一个 Keldysh 公式冒充微观独立方法。项目未证明原论文已完全复现。

## 发布核验

PR、最终 HEAD 的 GitHub Actions、合并后 main Actions 链接、SHA 与完整测试数必须从真实 GitHub 日志记录；不要使用本地通过代替远端通过。独立交接 ZIP 重新解压后运行总测试及 SHA256 校验，并注明离线文件和远端 Git blob 是否一致。远端冻结 JSON 按数值验收、非擅自宣称字节完全一致。
