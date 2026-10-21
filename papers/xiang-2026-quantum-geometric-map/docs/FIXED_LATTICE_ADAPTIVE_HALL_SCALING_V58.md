# v58 | Physics-first planar Hall bars: fixed lattice-site density, explicit multi-site contacts

**Scientific objective:** replace v57's sparse diamond-ring size family with square-lattice Nx×Ny *planar rectangular samples*, preserving local nearest-neighbour hopping t=1, lattice spacing a=1, on-site energy 0, signed magnetic flux per plaquette q B a²/hbar=-0.1044, chemical potential mu=.65 and T=.12. Fixed *lattice-site density* a^-2 is **not** the same as self-consistently fixed carrier density. The present model is spinless ordinary coherent orbital **charge Hall**, not the source-paper intrinsic quantum-geometric orbital MNHE/PHE.

## Four-terminal physics

Terminals L,T,R,B couple via distinct centered, mirror-symmetric, non-overlapping interior-edge contact site sets. Two explicitly different physically legitimate contact-scaling protocols are available:
1. **per_site**: each contacted site contributes the same constant local wideband Γ_ai=gamma; adding contacts increases each reservoir's total coupling.
2. **fixed_reservoir_sum**: each reservoir retains total sum_i Γ_ai=gamma; expanding contacts reduces the width per site.

For an explicitly positioned Peierls finite lattice and local diagonal energy/B-independent wideband Γ:
\[
G^r=[E-H(B)+i\sum_\alpha \Gamma_\alpha/2]^{-1},\qquad
T_{\alpha\beta}=\mathrm{Tr}(\Gamma_\alpha G^r\Gamma_\beta G^a)
=\sum_{ij}\Gamma_{\alpha i}\Gamma_{\beta j}|G_{ij}^r|^2.
\]
Take d_B G^r=G^r (d_BH) G^r, calculate the analytic derivative of *each* transmission, then self-consistently solve floating voltage-probe conditions I_T=I_B=0, V_L=+.5,V_R=-.5. In addition to Hall ratio V_H/ΔV, report the actual current-normalized Hall resistance R_xy=V_H/I_L in h/q² units and its derivative:
\[
\partial_BR_{xy}=[(\partial_BV_H)I_L-V_H(\partial_BI_L)]/I_L^2.
\]
**Neither device Hall voltage nor R_xy is automatically equal to intrinsic bulk Hall resistivity or conductivity.**

## A crucial independent numerical correction: 80-point Gaussian integration fails

Sharp mesoscopic transmission resonances can be narrower than the fixed Gauss node separation. For 5×5 with 3 contacted sites per edge and gamma=.6 per site, the 80-node Gauss result was V_H/ΔV=+0.01452986 and d_B(V_H/ΔV)=−0.06136070, while *independently converged adaptive* integration gives **+0.01484164** and **−0.06470123**. The earlier Gauss result must be **withdrawn as unconverged**.

The new adaptive spectral integrator evaluates the whole Fermi-window weighted T and d_BT vector, splits its energy domain at the finite-graph H(B) eigenvalues and uses adaptive scipy quad_vec, with positive epsrel/epsabs and fail-closed convergence. Independently tightened tolerance, both magnetic gauges and separate five-point full B differences verify the result. Fixed Gauss is explicitly labeled **uncertified**, never the default CLI output.

## Converged fixed-spacing size family

At mu=.65,T=.12,B=.12,q/hbar=−.87,t=a=1, full symmetric interior-edge contacts, per-site gamma=.6:

| Nx×Ny | Contacts L/T/R/B | Hall voltage / bias | d_B Hall voltage / bias | Hall resistance /(h/q²) |
| --- | --- | ---: | ---: | ---: |
| 5×5 | 3/3/3/3 | +0.01484164 | −0.06470123 | +0.01590859 |
| 7×5 | 3/5/3/5 | −0.01383948 | −0.17857225 | −0.01797448 |
| 7×7 | 5/5/5/5 | −0.03055589 | −0.46960112 | −0.02243102 |
| 9×5 | 3/7/3/7 | −0.01155765 | −0.14481043 | −0.01474357 |
| 9×9 | 7/7/7/7 | −0.05233421 | −0.37761076 | −0.02651467 |

At **the same 5×5** graph, merely changing full-contact per-site gamma=.6 to **fixed-reservoir total gamma=.6** changes V_H/ΔV from +0.01484164 to +0.02042262. Keeping per-site gamma=.6 but selecting **one central contact site per edge** instead produces +0.13034854 (adaptive reference). This is a contact effect at fixed local material model, not bulk-material intrinsic variation.

The size rows are **not monotonically convergent**; although the geometry and local model now use a consistent family, a finite fully coherent, contact-sensitive device series does **not** prove a thermodynamic bulk limit. Fixed mu does not imply identical finite open-system electronic carrier density.

## Independent acceptance

**25 new locally passing tests** cover exact fixed local hoppings/site density and plaquette flux; symmetric disjoint electrode sets; exact independent Caroli matrix trace vs site-pair formula; single-contact reduction to v57; full symmetric/Landau magnetic gauge invariance of *individual T_ab and d_BT_ab*; Onsager–Casimir, optical theorem and zero probe/total lead current; odd magnetic Hall and B0 zero for mirror-symmetric samples; proper Rxy normalization and quotient derivative; independent *five-point complete B voltage/current reruns*; two different lead-width normalization modes; rectangle rotation; correct a²B scaling; explicit numerical underconvergence negative gate; tighter adaptive quadrature, adaptive two-gauge and complete FD tests; invalid parameters rejected.

- **VERIFIED:** specified multi-site contact, coherent Peierls finite-device Landauer implementation, analytic field derivatives and numerical adaptive convergence.
- **CONDITIONAL:** finite-coherent sample Hall voltage and Rxy, which depend on size, contact coverage and self-energy.
- **UNRESOLVED:** thermodynamic bulk orbital MNHE/PHE, fixed carrier density, realistic disorder/dephasing, conserving microscopic scattering, nonlocal/energy dependent leads and full Wannier geometry. The original Xiang source PRB remains PARTIAL.

## Reproduce

    python -m unittest tests.test_planar_hall_scaling_v58 -v
    python scripts/reproduce_planar_hall_scaling_v58.py --output results/frozen/v58_fixed_lattice_full_contact.json
    python scripts/reproduce_planar_hall_scaling_v58.py --sizes 5x5,7x7,9x9 --coverage .4
    python scripts/reproduce_planar_hall_scaling_v58.py --sizes 5x5,7x7,9x9 --gamma-mode fixed_reservoir_sum

Three frozen JSON reports preserve full model inputs, contact site counts, reservoir Γ totals, Hall voltages, Rxy, B derivatives, adaptive quadrature errors and evaluation counts. A standalone v58 ZIP handoff contains the *complete Chinese scientific report*, the exact code/tests and all three full-resolution JSON benchmarks.

**v59:** energy-dependent / B-dependent realistic leads, disorder and phase-breaking, a separate periodic-bulk Kubo oracle and controlled physical limits. This v58 finite coherent square-lattice conventional Hall is a verified numerical device baseline, not yet a claim of the original paper's quantum-geometric magnetotransport.
