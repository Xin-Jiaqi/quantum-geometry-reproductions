# BSF × NMEE: Connection-operator interface (stage 1)

**Scope:** link existing symmetry and NMEE packages. This is **not** a new material prediction or a validated ferroelectric switching path. The first example is the horizontal-mirror connection between two *symmetry-related* candidate interface translations, and an **independently evaluated synthetic eight-band Hamiltonian** that realizes the same unitary transformation.

## Reuse, do not duplicate

- Symmetry: [`Xin-Jiaqi/group-theory-operations-toolkit`](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit), pinned at `551a0daf8d4ba5a5f24b2811694bfe438616de48`, v0.18.0 (`equivalent_interface_translation`, `point_group_operations`).
- NMEE: [`papers/jia-2026-nonlin-magnetoelectric-edelstein/`](../../papers/jia-2026-nonlin-magnetoelectric-edelstein), **unchanged** (both `printed` and `full_static` branches remain distinct).
- `bsf_nmee/connection.py` is only a *thin adapter*: explicitly carries fractional and Cartesian rotations; it does **not** reproduce a point group database, copy a stacking classification engine, or redefine NMEE physics.

## Connection and predicted tensor relation

For `m_001` of the hexagonal parent `6/mmm`, the toolkit's interface map gives

\[
\tau_1=(1/3,2/3)\xrightarrow{m_z}\tau_2=(2/3,1/3),
\quad \mathrm{top}\leftrightarrow\mathrm{bottom}.
\]

**Precondition:** the mirror must act on the *complete monolayer atomic motif* in the specified setting, not just on its Bravais lattice. The orbit map alone does not demonstrate energetic equivalence or switchability for an arbitrary material.

For unitary spatial operation `R_cart` and tensor ordering `Gamma[spin, B, E]`, spin and B are both **axial**, E is **polar**. Hence their two `det(R)` factors cancel:

\[
\Gamma^{\mathrm{target},L'}_{\alpha\gamma;\beta}
=R_{\alpha a}R_{\gamma c}R_{\beta b}
  \Gamma^{\mathrm{source},L}_{ac;b},
\quad L'=p_R(L).
\]

With `m_z=diag(1,1,-1)`, **for the component `Gamma_{xz;x}` only**:

- `Gamma_total` is **odd** between mirror-related stackings;
- `Gamma_top - Gamma_bottom` is **even** between mirror-related stackings;
- `Gamma_top(tau2) = -Gamma_bottom(tau1)` and vice versa.

Do not confuse this with pure layer relabeling `tau_x`, which does *not* act on spin and therefore has different parity.

## Run

From the repository root:

```sh
python -m pip install -e papers/jia-2026-nonlin-magnetoelectric-edelstein
python -m pip install 'git+https://github.com/Xin-Jiaqi/group-theory-operations-toolkit.git@551a0daf8d4ba5a5f24b2811694bfe438616de48'
cd studies/bsf-nmee-connection
python -m unittest discover -s tests -v
python scripts/mirror_covariance_nmee.py --mesh 11
```

CI runs against **the real upstream toolkit**, not mock results. `tests/test_tensor_contract.py` additionally checks the pure physical-tensor transformation contract without an external dependency.

## First Hamiltonian check

`mirror_covariance_nmee.py` takes a synthetic eight-band model based on two coupled Kane–Mele copies and **generates its mirror partner by the physical unitary** `U = tau_x (layer) ⊗ I_sublattice ⊗ sigma_z (spin)` (global spinor phase irrelevant). Its shifts are *labels of an abstract connection*, not a mapping from atomic coordinates to hopping functions `T(k,tau)`. No sliding path, polarization energy or barrier is computed.

The script calculates `Gamma_xz;x` independently in both configurations, for top, bottom and total spin, on a uniform 11×11 BZ grid, using the **existing PRL code unchanged** and **both** response branches. It verifies the absolute mirror/particle-layer relation; these small-grid values are **not converged response amplitudes**.

The result file `results/mirror_11.json` records the mesh and parameter conventions, model energies, two response branches and numerical residuals. To convert this stage into a genuine BSF sliding model, next provide `T(k,tau)` from structure or an explicitly symmetric orbital model and verify covariance of H(k,tau) and all observables, **before** claiming polarization switching.

## Important limits

- Only **unitary**, layer-normal-preserving operations. Antiunitary magnetic operations and nonsymmorphic atomic mappings are not yet covered by this adapter.
- NMEE's existing finite-T driven DC vs complete static Eq. (9) ambiguity remains unresolved; these response branches must never be combined.
- The original paper reproduction and its frozen tests are **not modified**.

## Stage 2: microscopic relative-sliding hopping (independent preliminary model)

[`bsf_nmee/sliding_tight_binding.py`](bsf_nmee/sliding_tight_binding.py) now provides a *different* eight-band model in which interlayer hopping is independently generated from real-space orbital displacements for each stacking `shift`. It is **not** the stage-1 abstract unitary-rotated toy model. The hopping uses five-shell lattice sums and the *full* interorbital Bloch phase, including intra-cell orbital positions and relative displacement. Both spin-orbit-reversed monolayers are related by a physical `m_z` mirror with top/bottom exchange.

[`bsf_nmee/valley_quadrature.py`](bsf_nmee/valley_quadrature.py) integrates small valley disks with nonlinear radial resolution and the rest of the BZ on a uniform grid. A static zero-temperature finite-Zeeman-field derivative of the **interband spin-electric susceptibility** is evaluated using the already-reviewed NMEE routines. This is *not* the paper's nonequilibrium DC branch. See [`docs/STAGE2_PHYSICS.md`](docs/STAGE2_PHYSICS.md) and `results/sliding_stage2_validated.json` / `results/sliding_convergence_extended.json` for numbers and physical caveats.

A key finding in this particular synthetic model is a nonzero **net** static `Gamma_xz;x` induced by interlayer hopping between otherwise mirror-compensated layers. At relative shift `(0.23,0.41)`, `t_perp=0.028 eV`, zero layer potential difference, the well-refined reference is approximately `+2.03110e-4 (nm V T)^-1`; `t_perp=0` yields numerical zero. With only same-sublattice tunneling, it shrinks to `+2.83e-8` in this model. Neither result demonstrates ferroelectric stability or a universal enhancement, and a prior informal `5.80e-5` estimate should **not** be treated as a frozen result for this microscopically specified model.

Reproduction with the paper package installed, from this study directory:

```bash
python -m unittest discover -s tests -v
python scripts/run_stage2_benchmark.py   # ~30 s in a single-thread BLAS environment
```

No original PRL reproduction modules are modified. True static-response amplitude depends on the microscopic hopping, SOC, and position convention; the existing final-DC-formula disagreement remains unresolved.

## Stage 3: sign change *within* a gapped continuous sliding path

At fixed `t_perp=20 meV`, zero layer potential difference and a fixed chemical
potential, the independently calculated **static insulating** `Gamma_xz;x`
changes sign as `tau(u) = u(0.23, 0.41)` is varied from `u=0` to `u=1`.
The sampled bandgap does not close along this path. The zero is near
`u ≈ 0.1974`, **not** a symmetry-forced or topological transition. Each
of K and K′ changes sign independently; their contributions are equal, not
mutually canceling. BZ-averaged Frobenius-norm matching confirms that omitting
cross-sublattice tunneling suppresses the response by orders of magnitude in
this model. These are synthetic static susceptibility results, not ferroelectric
switchability or the author's driven-DC NMEE prediction.

See [stage-3 physics report](docs/STAGE3_PHYSICS.md),
[complete numeric scan](results/sliding_stage3_sign_reversal.json), and
[`scripts/run_stage3_sign_reversal.py`](scripts/run_stage3_sign_reversal.py).
The frozen Jia et al. PRL reproduction folder has **not been modified**.

Run from this directory after installation of the parent paper module:

```bash
OPENBLAS_NUM_THREADS=1 python scripts/run_stage3_sign_reversal.py
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
```


## Stage 8a: literature-only material candidate prescreen

The first candidate shortlist is documented in [material-screening/STAGE8_CANDIDATE_PRESCREEN.md](material-screening/STAGE8_CANDIDATE_PRESCREEN.md), with a machine-readable evidence/unknowns matrix in [material-screening/STAGE8_CANDIDATE_MATRIX.json](material-screening/STAGE8_CANDIDATE_MATRIX.json).

- **Primary:** 3R-WSe2 homobilayer (strong SOC evidence, semiconducting band structure, reported AB/BA ferroelectric domains).
- **Control:** 3R-MoS2 homobilayer (well-studied switching and the same clean 1H-TMD stacking-symmetry route).
- **Secondary:** 3R-WS2.
- **Exploratory/separate classes:** Cd2Se2, 1T-PtS2, and MoS2/WS2 heterobilayer.

This is a **literature and symmetry pre-screen only**. No DFT, Wannierization, or real-material NMEE has been calculated. The next gate is to retrieve exact structures and test the complete atomic motif with the existing group-theory toolkit before planning any electronic-structure run.
