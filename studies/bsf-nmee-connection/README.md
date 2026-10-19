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
