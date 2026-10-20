# Stage 7: K/K' valley k·p reproduction of the sliding NMEE mechanism

**Scope:** Eight-band synthetic sliding Kane–Mele bilayer; zero-temperature, insulating **equilibrium static** mixed spin–momentum response. This is **not** the PRL's unresolved driven-DC response and is **not** a Wannier fit or demonstrated ferroelectric.

## Scientific finding

An independently evaluated second-order Taylor Hamiltonian about the two valleys,

\[
H_\nu(\mathbf q;u)=H_\nu(0;u)+q_i V_{\nu i}(u)+\tfrac12q_iq_j W_{\nu ij}(u),\quad V_i=\partial_{k_i}H,\quad W_{ij}=\partial_{k_i}\partial_{k_j}H,
\]

reproduces the **valley-local static NMEE**, its weak-coupling AA/BB–AB/BA interference coefficient, the sliding sign-reversal location, and the previously derived C3/mirror sliding-order constraints. Unlike a fitted response curve, **all k·p coefficients are obtained by differentiating the microscopic Hamiltonian** and velocities are consistently differentiated from the same Taylor Hamiltonian.

The local model still has **eight bands**. It is a controlled low-momentum expansion, **not** yet an analytic 2/4-band downfolding or a universal closed-form coefficient.

## Reproducing the checks

From this directory, with the existing paper `nmee_lab` installed and `numpy` available:

```bash
python -m unittest discover -s tests -v
python scripts/run_stage7_kp.py --output results/sliding_stage7_valley_kp.json
python scripts/run_stage7_capture.py
python scripts/run_stage7_channels.py
python scripts/run_stage7_robust.py
python scripts/run_stage7_root_and_trace.py
python scripts/plot_stage7.py
```

The full scan may be relatively expensive. The frozen numerical results and physical constraints are detailed in [`docs/STAGE7_PHYSICS.md`](docs/STAGE7_PHYSICS.md), and provenance/checksums in [`results/STAGE7_MANIFEST.json`](results/STAGE7_MANIFEST.json).

## Link to material Wannier calculations (concept only)

Given consistent Wannier orbitals and orbital positions, obtain `H(k)`, `∂k H(k)` and `∂k∂k H(k)` around K/K' (or any relevant direct-gap extrema) from the **same Bloch gauge**. A physical response additionally requires projected **spin and position operators/Berry connections**. `hr.dat` alone is not enough; no first-principles calculation has been performed here.
