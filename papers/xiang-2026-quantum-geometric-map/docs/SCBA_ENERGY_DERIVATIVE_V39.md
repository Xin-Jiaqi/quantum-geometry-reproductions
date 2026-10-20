# v39 — Energy-dependent SCBA, density Ward derivative, and an exact-disorder sea benchmark

**Scientific acceptance:** PASS only for a precisely defined *two-band, scalar-white-noise, noncrossing SCBA energy-derivative method*. **The complete disordered Kubo–Bastin/Středa conductivity and the published nonlinear mixed-field EB response remain OPEN.** This is an independent extension of the v38 clean Fermi-sea test rather than a result extracted from Xiang *et al.*, PRB **113**, L201406 (2026).

## 1. Physical model and strict protocol

The deliberately separate test system is the QWZ Bloch Hamiltonian

\[
H(\mathbf k)=\sin k_x\,\sigma_x+\sin k_y\,\sigma_y+(m+\cos k_x+\cos k_y)\,\sigma_z,
\]

with scalar uncorrelated onsite disorder of variance \(W\). Retarded SCBA on a complete periodic momentum mesh is defined by

\[
G^R(\epsilon,\mathbf k)=[(\epsilon+i\eta)I-H(\mathbf k)-\Sigma^R(\epsilon)]^{-1},\qquad
\Sigma^R(\epsilon)=W\langle G^R(\epsilon,\mathbf k)\rangle_{\mathbf k}.
\]

Here \(\eta>0\) is an *external Green-function regulator*, not an independently verified physical collision rate. The model uses \(e=\hbar=1\), momentum in lattice radians, and the v38 dimensionless convention for the normalized sea response \(\sigma/(e^2/h)\).

## 2. New Ward-consistent energy derivative

Differentiating the **self-consistent** equation, rather than pretending that \(\Sigma\) is energy independent, gives

\[
\partial_\epsilon G^R=-G^R\,\Gamma^{RR}_0\,G^R,\qquad
\Gamma^{RR}_0=I-\partial_\epsilon\Sigma^R,
\]

\[
\boxed{\Gamma^{RR}_0=I+W\langle G^R\Gamma^{RR}_0G^R\rangle_{\mathbf k}.}
\]

The program solves this final linear equation in the full four-dimensional Pauli basis, independently converges \(\Sigma(E+h)\) and \(\Sigma(E-h)\), and compares the resulting finite-difference \(\partial_E\Sigma\) and \(\partial_E G\) with the Ward derivative. The **RR density vertex \(\Gamma_0^{RR}\)** is not a replacement for any transverse electromagnetic current vertex; this distinction matters when claiming a conserving Hall conductivity.

## 3. Numerical benchmark (fixed model, not paper's Fig. 2)

For \(m=-1\), \(21^2\) mesh, \(W=0.12\), \(\eta=0.12\), energy derivative tested at \(E=0.7\), and real-energy quadrature with 256 nodes over \([-7,\mu]\):

| Quantity | Numerical value |
| --- | ---: |
| \(\max |\Gamma^{RR}_0-I|\) at \(E=0.7\) | 0.32958193 |
| finite-difference self-energy derivative discrepancy | \(9.37\times10^{-8}\) |
| finite-difference propagator derivative discrepancy | \(1.63\times10^{-6}\) |
| insulating \(\mu=0\), dressed *bare-current* sea kernel | −0.84151865 |
| insulating \(\mu=0\), naive \(\partial_E G=-G^2\) kernel | −0.81484870 |
| insulating derivative correction | **−0.02666995** |
| metallic \(\mu=1.1\), dressed *bare-current* sea kernel | +0.97431959 |
| metallic \(\mu=1.1\), naive kernel | +1.09769378 |
| metallic derivative correction | **−0.12337420** |

These are **two versions of the same defined mathematical sea kernel**, not two physically complete Hall conductivities. The larger metallic difference is sensitive to finite broadening and the Fermi-level region. Our direct real-energy integration is explicitly converged against a denser quadrature, and the \(W=0\) limit reproduces v38's analytic two-pole primitive.

## 4. Truly independent exact finite-disorder check

For a finite \(L\times L\) torus, we Fourier-rotate the clean Bloch Hamiltonian and both exact velocity matrices into a \(2L^2\)-dimensional Hilbert space. We then add **explicit random scalar onsite potentials** with standard deviation \(\sqrt W\), exactly diagonalize each disordered Hamiltonian, and evaluate the same *bare-current sea trace* from all unordered eigenstate pairs using the independently derived two-pole primitive.

The \(W=0\) exact finite-torus result agrees with the same finite \(k\)-mesh clean v38 primitive to floating-point precision, providing an independent check of Fourier normalization and pairwise spectral algebra. For \(L=4\), \(W=0.12\), six fixed seeds, the exact finite-system average is **−1.05681902**, with sample standard deviation **0.02397452**. It differs from the thermodynamic-mesh SCBA sea result \(-0.84151865\); this is **not an SCBA failure diagnosis**, because mesh sizes, finite-size boundary sampling, disorder ensemble and missing transverse ladder contributions are not matched. We do not fit or assert equality.

## 5. Explicitly unfulfilled gates

- Transverse current vertex and relevant gauge/contact or seagull diagrams consistent with the disorder self-energy.
- Full Kubo–Bastin/Středa I (Fermi-surface) and II (sea) response for *metallic, disordered* Hall physics, with correctly specified response limits.
- Crossed scattering processes and physical \(q,\omega\to0\) order of limits.
- General orbital three-band pieces and magnetotransport in real Wannier materials.
- Original paper's **second-order \(EB\) DC PHE** and the unresolved Eq. (11)/Fig. 2(b) ~\(2\pi\) normalization discrepancy.

To reproduce, in the PRB module:

```bash
python -m unittest discover -s tests -p 'test_scba_streda_v39.py' -v
python scripts/reproduce_scba_streda_v39.py --mesh 21 --W .12 --eta .12 --quadrature 256
```

The companion JSON file `results/frozen/v39_benchmark.json` stores parameter labels, residuals and model snapshots; no publisher figures or original data are included.
