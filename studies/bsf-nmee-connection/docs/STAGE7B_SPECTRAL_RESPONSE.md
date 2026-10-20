# Stage 7b: exact k-adapted four-band subspace is NOT enough for static NMEE

**Scope:** synthetic sliding 8-band Kane–Mele, $T=0$, insulating equilibrium/static finite-Zeeman susceptibility. NOT the unsettled driven-DC intrinsic NMEE, DFT, or stable sliding ferroelectricity.

## Motivation and exact decomposition

Stage 7's full 8-band quadratic K/K' expansion accurately reproduces the signed valley NMEE. The Stage 7 SW gate showed that a **fixed** decoupled K-center four-band basis need not be perturbative across valley disks.

To avoid unfairly blaming a badly chosen fixed basis, Stage 7b constructs the **exact, k-dependent middle-four eigenspace** (energy-ordered bands 2,3,4,5). It then partitions the *original eight-band static response* into four **disjoint, exact** interband sectors at both finite $\pm B_z$:
1. middle occupied (2,3) $\to$ middle empty (4,5), called LL;
2. remote occupied (0,1) $\to$ middle empty;
3. middle occupied $\to$ remote empty (6,7);
4. remote occupied $\to$ remote empty.

The sum is algebraically exactly the parent occupied-to-empty response. This band-*sector* diagnostic is **not** the same as the fixed-basis same-/cross-sublattice hopping decomposition, nor is an energy-resolved partial term an independent gauge-unrestricted observable.

## Independent numbers

Two nonoverlapping radius **0.5 nm^-1** disks at K and K'; signed integral of $\Gamma_{xz;x}$ using native $H(k)$ and $\partial_{k_x}H$ with $\pm0.25$ T finite-Zeeman differentiation (26 radial Gauss-Legendre nodes, 72 angular nodes per valley):

| u | $t_\perp$ | full eight-band | LL only | involving remote states | LL has wrong sign? |
|---:|---:|---:|---:|---:|---|
| .1 | 5 meV | -6.916904e-7 | -5.664028e-7 | -1.252875e-7 | No |
| .1 | 20 meV | -1.704537e-5 | +1.099542e-5 | -2.804079e-5 | **Yes** |
| .4 | 5 meV | +2.164836e-6 | -1.330624e-5 | +1.547108e-5 | **Yes** |
| .4 | 20 meV | +9.843337e-5 | +1.077986e-4 | -9.365245e-6 | No |

At $(u=.4,t_\perp=5\ \mathrm{meV})$ the four components are:
LL $-1.330624\times10^{-5}$, remote-occ/low-empty $+3.74739\times10^{-6}$, low-occ/remote-empty $+1.60859\times10^{-5}$, remote-occ/remote-empty $-4.36218\times10^{-6}$. All four channels can matter.

Independent all-band numbers agree with the *previous* Stage 7 disk integral for $(u=.1,t=5\mathrm{meV},R=.5)$ (old $-6.916784\times10^{-7}$ vs present $-6.916904\times10^{-7}$), and for $(u=.1,t=20\mathrm{meV},R=.5)$ (old $-1.704504\times10^{-5}$ vs present $-1.704537\times10^{-5}$). These are **disk-local** checks, not new BZ convergence claims.

Numerical convergence: $(u=.1,t=20\mathrm{meV},R=.5)$ full -1.704537449e-5 at 26×72 versus -1.704537456e-5 at 40×96; LL +1.099541750e-5 versus +1.099541742e-5. Changing $\pm B$ step from .25 to .125 or .5 T does not change signs. Likewise at $(u=.4,t=5\mathrm{meV},R=.75)$ full remains positive and LL negative at all three field steps. The independent implementation checks $\partial_kH$ by finite differences.

## Spectral projector continuity

On 30 radial × 60 angular k-points in each radius 0.75 nm^-1 disk, evaluate $P(k)$ and the decoupled four-band projector $P_0(k)$ **at the same k**. The smallest sampled adjacent middle-four-to-remote gap and smallest singular value squared of $P_0^\dagger P$ (principal-angle cosine squared) are:

| u | $t_\perp$ | min outer gap | min cos² |
|---:|---:|---:|---:|
| .1 | 5 meV | 11.54 meV | .9754 |
| .4 | 5 meV | 12.47 meV | .9821 |
| .1 | 20 meV | 8.48 meV | 1.97e-7 |
| .4 | 20 meV | 2.61 meV | 2.66e-5 |

**Especially important:** the weak-5-meV cases have a very well aligned k-dependent projector but the LL-only response can still have the WRONG sign. Downfolding error is not only a poor wavefunction-subspace match; virtual remote-band matrix elements are essential. A finite positive sampled gap does not establish a rigorous whole-disk spectral bound. The stronger-20-meV cases require revisiting the reference basis or keeping more states.

## Covariant velocity: necessary but not sufficient

For any locally smooth exact invariant middle-four frame $W^\dagger W=I$, $h=W^\dagger HW$ and $\mathcal A_j=W^\dagger\partial_{k_j}W$:

\[
\boxed{W^\dagger(\partial_{k_j}H)W=\partial_{k_j}h+[\mathcal A_j,h].}
\]

Thus, diagonalizing $h$ at every $k$ and using only derivatives of its eigenvalues discards interband velocity matrix elements. With a polar-aligned local frame at $q=(.05,.02)$ nm^-1, central difference step .001 nm^-1, the checked covariant-velocity residual is $\leq1.2\times10^{-5}$ eV nm, while exact middle-band velocity offdiagonal matrix elements are approximately 0.096–0.145 eV nm. This is a local **identity and finite-difference check**, not an integrated four-band response result.

**Crucial caveat:** even $W^\dagger vW$ together with $W^\dagger SW$ reproduces *only* the LL interband terms. Remote occupied/empty transitions remain in the full Kubo curvature. A genuinely equivalent four-band response needs response-level **dressed vertices/field-dependent downfolding and associated geometric/contact terms**, not just the Hamiltonian or a covariant-velocity correction.

## Reproduction and next single task

- New module: bsf_nmee/spectral_sector.py (uses ORIGINAL 8-band source and velocities).
- Tests: tests/test_stage7b_spectral_sector.py (derivative check, spectrum/occupation, LL sign failure and remote contribution; discovered by existing study CI).
- Independent no-project-import implementation, full JSON/CSV evidence, numerical convergence, local projector/covariant-velocity checks and tests are archived separately as BSF_NMEE_STAGE7B_20261020.zip.

**Next task:** work out the response-preserving reduced resolvent/vertices (or keep validated eight bands), with quantitative control of the omitted-band contribution. Do **not** claim a successful four-band NMEE theory. Preserve prior frozen results; no DFT.
