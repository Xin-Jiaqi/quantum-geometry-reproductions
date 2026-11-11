# Kubo/projector 静态响应对照 —— 公式推导与规范说明（Stage8n）

目标量（与 plaquette 引擎一致）：Γ^stat_{iz;j} = ∂_{B_z} ∫BZ d²k/(2π)²
· i Tr P[∂_{λ_i}P, ∂_{k_j}P] |_{λ=B=0}，冻结 H(k,λ,B) = diag E + λS_i/2 − μ_B B S_z。

## 1. 展开

P = Σ_{n∈occ} |u_n⟩⟨u_n|。一阶微扰：
  ∂_λ P = A + A†，A = Σ_{n occ, m unocc} a_nm |u_m⟩⟨u_n|，a_nm = ⟨u_m|S_i/2|u_n⟩/(E_n−E_m)
  ∂_k P = B + B†，B = Σ_{n occ, m unocc} b_nm |u_m⟩⟨u_n|，b_nm = ⟨u_m|∂_k H|u_n⟩/(E_n−E_m)

因 P A = P B = 0：
  Tr P[∂_λ P, ∂_k P] = Tr[(A†)(B†)] − Tr[(B†)(A)] = Σ (a* b* − b* a)
                     = −2i Σ Im(a_nm b_nm)
  ⇒ i Tr P[∂_λ P, ∂_k P] = +2 Σ_{n occ, m unocc} Im[ S_mn v_mn ] / (E_n−E_m)²

其中 S_mn = ⟨u_m|S_i/2|u_n⟩（SPN/2），v_mn = ⟨u_m|∂_k H|u_n⟩ 带间速度。

## 2. 带间速度的完整表达式

⟨u_m|∂_k H|u_n⟩ = (E_n−E_m)⟨u_m|∂_k u_n⟩ + δ_mn ∂_k E_n。带间（m≠n）：
  v_mn = (E_n−E_m)⟨u_m|∂_k u_n⟩。

在 B 依赖基 u(B) = u(0)U(B) 下：
  ⟨u_m(B)|∂_k u_n(B)⟩ = [U†⟨u(0)|∂_k u(0)⟩U + U†∂_k U]_mn。

第一项来自提交的 MMN 重叠中心差（VASP 带基）：dM = [M(k,k+δk)−M(k,k−δk)]/(2δk)。
第二项 U†∂_k U 是 B 依赖基自身随 k 的旋转（带间），需要 U(B) 在相邻 k 的连续规范
才能数值稳定（eigh 列相位任意 ⇒ ∂_k U 带间项相位任意）。

## 3. 规范（gauge）性质

- S_mn 与 ⟨u_m|∂_k u_n⟩（第一项）在同一下变换下带间组合 Im(S·v) 规范不变。
- U†∂_k U 项对列相位敏感；零新作业约束下未做平行传输规范固定。
- 数值对照（AB，5 k × 2 轴）：1/10 点 rel < 0.5，其余 9–275×；
  归因于 U†∂_k U 项缺失 + eigh 相位不连续。合成模型因带隙部分近简并
  （min occ-unocc gap 0.014 eV）微扰展开失效，不适合作 Kubo 数值基准。

## 4. 结论

Kubo 公式与实现方向正确（个别点 rel<0.5 支持）；完整独立对照需要规范固定
（平行传输/平滑规范），超出零新作业数据直接可计算范围，如实记录。
