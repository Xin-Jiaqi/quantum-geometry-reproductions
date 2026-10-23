# Stage8m 粗粒度响应估计 —— 运行文档

## 管线

```
VASP 6.4.3 ALGO=None 后处理（LWANNIER90 + LWRITE_SPN，ISTART=1, NELM=1, LDIAG=F）
   -> wannier90.mmn / wannier90.spn / wannier90.eig / wannier90.win
stage8m_response.py   # 响应估计（DeterminantPlaquette 引擎 + 441 点 BZ 求和）
stage8m_analyze.py    # 收敛/零值/对称性后分析
```

## 引擎（stage8m_plaq.py）

- 物理量：`Gamma^{stat}_{iz;j} = d/dB_z [ i Tr P[ dP/dlam_i, dP/dk_j ] 的 BZ 积分 ]`，
  frozen-KS 哈密顿量 `H(k,lam,B) = diag E(k) + lam*S_i/2 - muB*B*S_z`（g=2、spin-only Zeeman）。
- 中心对称参数环（与 Stage7i 独立引擎一致，rel < 5e-8 回归）。
- mmn 模式：k 边用两步链 `k1 -> k_c -> k2`（wannier90.mmn 只提供单步最近邻重叠），
  O(dk^2) 精度；链接全带重叠按窗口截断。
- SPN 存 Pauli 矩阵元（对角 ~0.817>0.5 验证为 sigma，非 sigma/2），源乘 1/2。

## 运行（远端）

```bash
python3 stage8m_response.py --dir <run_dir> \
  --windows 36,48,64,80,96 --Blist 0.125,0.25 --dllist 1e-4,2e-4,4e-4 \
  --source x,y --occ 36 --kmesh 21,21 --jobs 8 --out stage8m_AB_response.json
```

`--jobs N` 用进程池并行（fork 上下文、共享数组经 initializer 注入一次）。

## 输出

- 每配置 `axis_x/axis_y`：`Gamma_mesh_mean`（网格均值，1/(eV·T·Ang^-1)）、
  `Gamma_BZ_Ang2`（全 BZ 估计 = weight·sum，1/(eV·T·Ang^2)）、
  `min_gap_eV`（fail-closed）、`min_sval`（polar 分解最小奇异值）。
- `sym_residual_abs = |Gamma_x + Gamma_y|`。
- `_meta`：单位约定、kmesh、SPN Hermiticity 偏差、S_z 对角上界。

## 判读

- GO 需要：非零且 B/dl/window 收敛；glide total null 成立；AB/glide 差异可解释；
  引擎回归（Stage7i）与 fail-closed 测试全绿。
- SI 换算：`Gamma_BZ_Ang2 × (1e-9 m/Ang)^2 / (1.602176634e-19 J/eV) = m^2/(J·T)`。

## 测试

```bash
python3 -m unittest discover -s tests -p 'test_stage8m*'
# 全部 173 测试（skipped=7）通过
```
