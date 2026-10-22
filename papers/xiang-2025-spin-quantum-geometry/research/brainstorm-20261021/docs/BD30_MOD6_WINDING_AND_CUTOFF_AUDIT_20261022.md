# BD30 | Operator-parity winding modulo six and Pd₂Te₂ numerical-continuity correction

**2026-10-22 | THEOREM/ALGEBRA_PASS; TOY_SINGLE_POINT_NONZERO_REPRODUCED; MICROSCOPIC_VORTEX_TRANSITION_NOT_YET_VALIDATED; NOVELTY_UNRESOLVED.** No new DFT/HPC work.

## 1. Scientific question: how can an individual result grow into a general theory?

Rather than reject a general layer-resolved response-field theory because individual ferroelectric/Edelstein/optical-vortex examples exist, identify an **exact operator-specific proposition**. The 2023 BSF theory and 2026 LEE PRL are *direct foundational theory* for stacking symmetry and layer-resolved current-induced spin; existing specific materials or shift-current vortices are not automatic priority conflicts for a different S2/S3 general theorem. The user requests this proposition-relative review standard.

Fix physical source to **uniform in-plane DC charge current** and output to layer-even/odd **out-of-plane axial spin density**:
\[
\delta S_\pm^z(\mathbf u)=\mathbf F_\pm(\mathbf u)\cdot\mathbf j_\parallel,\qquad S_\pm^z=S_T^z\pm S_B^z.
\]
For an actual identical-layer family having \(C_{3z}\), a vertical mirror and a horizontal **layer-swapping** \(m_z\), let \(z=u_x+i u_y\) and \(f_\pm(z)=F_{\pm,x}+iF_{\pm,y}\). The crystal family's operations map the entire registry family, not necessarily stabilize a given displaced structure:
\[
f_\pm(e^{2\pi i/3}z)=e^{2\pi i/3} f_\pm(z),\quad f_-(-z)=-f_-(z),\quad f_+(-z)=+f_+(z).
\]
The vertical mirror adds \(f_\pm(-\bar z)=\overline{f_\pm(z)}\) in the selected basis.

### Stronger, loop-level selection law

For **continuous, real, nowhere-vanishing** \(\mathbf F_\pm\) on a circular path centered at the high-symmetry registry, the integer winding \(w=(2\pi)^{-1}\oint d\arg(f)\) satisfies:
- \(C_3\) covariance: \(w\equiv1 \pmod3\).
- Layer-swap oddness/evenness: \(w_-\) odd, \(w_+\) even.
- Therefore
\[
\boxed{w_-\in1+6\mathbb Z,\qquad w_+\in4+6\mathbb Z.}
\]

This constrains **all permitted winding values**, not just the lowest-order \(+1\) and \(-2\) obtained from \(f_-\sim i\gamma z\) and \(f_+\sim\beta\bar z^2\). It is the application of **standard equivariant-map topology**, **NOT** a newly invented universal topology theorem. The effective sixfold restriction is an action **on registry parameter space**, not a claim that the physical crystal has \(C_6\). Finite-frequency complex tensor responses and magnetic antiunitary symmetries require distinct treatment.

### A mathematical six-satellite mechanism, not a material claim

The symmetry-allowed real-coefficient example
\[
f_+(z)=\beta\bar z^2+c z^4
\]
has a \(-2\) local index at \(z=0\) (if \(\beta\ne0\)) and **six index-\(+1\) satellite zeros** at radius \(\sqrt{|\beta/c|}\), hence an outer winding \(+4\). The six satellites are an explicit demonstration of **winding reallocation in multiples of six**, not evidence that real Pd₂Te₂ has them. The [standalone exact test](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/stacking_z6_winding_gate.py) checks \(w=-5,+1,-2,+4\), all assumed symmetries and the local \(+1\) Jacobians.

## 2. Independent re-evaluation of the real-BiDB-geometry Te-p SOC toy

Using the same BiDB id 980 atomic coordinates and 4 Te p sites × spin (24 bands), independent local reconstruction reproduces the historical **hard-cutoff** fixed-current Boltzmann numbers at \(u_x=.10\) Å:
- \(8^2\): \(\chi_{S_-^z,j_y}=-0.0065826711\) and \(\chi_{S_+^z,j_x}=-0.0000623071\);
- \(12^2\): \(-0.0065247867\) and \(-0.0000857695\).
These are **arbitrary unfitted model units**, with Pd d orbitals OMITTED, **not** a true Pd₂Te₂ material prediction. At 10² over illustrative chemical potentials \(\mu=1.0,1.25,1.5,1.75,2.0\) eV the layer-odd signal remains nonzero and same sign (about −0.0024 through −0.0067); SOC=0 returns numerical zero. This strengthens *model existence* but cannot establish a new physical mechanism.

### **Critical numerical bug: the hard cutoff invalidates global smooth registry topology**

The original script deletes all hopping at distance \(r>7.2\) Å; at \(u_x=0.105\to0.110\) Å, the pair list suddenly drops **74→70** (four Te-p bonds with pre-deletion amplitude ~0.0055 in model units). This makes \(H(\mathbf u)\), and potentially the derived \(\mathbf F(\mathbf u)\), **discontinuous** in the very registry range used for earlier contour plots. Thus the *hard-cutoff global winding data cannot, by itself, certify a continuous topological field*. This does **not** retract verified single-point nonzero values or symmetry-forced zeros.

The [GitHub model script](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/pd2te2_te_p_soc_boltzmann_gate.py) now supports optional \`--smooth-cutoff-start 6.5\`: a \(C^2\) quintic taper from 6.5 to 7.2 Å. Default **hard** cutoff is preserved *only* for archived numerical reproduction and prints a topology warning; future registry topology MUST use smooth cutoff. Original BiDB is read-only and not committed.

### Sensitive weak-even coefficient: do not report a six-vortex transition yet

At toy \(\mu=1.31\) eV, smooth cutoff, circle radius \(0.125\) Å:
| \(k\)-mesh / contour angles | computed layer-even winding | status |
|---|---:|---|
| 10×10 / 12 | **+4** | min norm \(5.8\times10^{-7}\) (near numerical cancellation) |
| 12×12 / 12 | **−2** | changed by \(k\)-quadrature |
| 12×12 / 24 | **−2** | denser angular sampling gives same result; max adjacent phase increment ~0.581 rad |
| 12×12 / 12 at \(\mu=1.50\) eV | **−2** | much stronger signal, min norm ~\(1.0\times10^{-4}\) |

On the \(30^\circ\) symmetry ray, the model does have near-zero roots \(r=0.11263\) Å (10²) and \(r=0.11498\) Å (12²), and an additional \(0^\circ\) ray root \(r\approx0.10413\) Å appears at 12²; local Jacobians of investigated roots are positive. Their complete index accounting, possible additional roots, and gauge/near-degeneracy stability have **not** been established; coarse angular sampling can alias higher windings. Therefore neither the position/count of physical satellite vortices nor a chemical-potential-driven topological transition is verified. Do not promote these preliminary branches into the Abstract.

A \(C^2\) taper-start comparison 6.3/6.5/6.7 Å changes weak-even amplitudes; different shifted \(k\) meshes also change near-cancellation values. This is a genuine **evidence blocker**, not a reason to reject the S2 general-theory direction.

## 3. Closest prior art: separate exact theory competition from related examples

**LEE PRL 137, 076401 (2026), DOI 10.1103/qhvs-pfkq:** the [public full text](https://arxiv.org/html/2607.17209) already contains the homogeneous second-rank axial response tensor, Boltzmann fixed-current formula (Eqs. 1–2), minimal model Eq. 3, bilayer intra-/interlayer group mapping and Type I/II feasibility for all 80 LGs. The user-supplied full research notes \`Layer Edelstein Effect.md\` additionally summarize its SM Sec. I (ordinary axial tensor), Sec. III (80 LG stacking table), Sec. IV (C2DB) and V–VI (MoSSe models). Those notes do **not** record an explicit **continuous-registry real-spin response winding congruence modulo six or six-satellite theorem**, but the **original LEE SM PDF was not independently obtained**; absence from private notes cannot establish absence from source. **Exact SM novelty status = UNRESOLVED.**

**Ji et al., PRL 130, 146801 (2023), BSF:** establishes configuration-space G/A/B/C/M/N and stacking symmetry. It is direct foundations for the geometry but is not automatically a prior calculation of this exact nonequilibrium operator-dependent spin-response topology. Same rule for optically driven shift-current vortices and field-driven Edelstein theories: assess exact operator, source, dynamical regime, abstraction level and proposition before calling a direct competitor.

## 4. Correct next scientific gate

1. **Numeric:** use smooth (at least \(C^2\)) hopping; diagnose degeneracies and the band-diagonal Boltzmann gauge issue; full 2D root search + refined phase contours and local Jacobians on two or more genuinely converged \(k\) grids. A 12-angle contour alone is not sufficient when \(w\) could be \(+10\).
2. **Microscopic/material:** find a source-parametrized Pd-d+Te-p SOC electronic structure, or replace with a better-validated bilayer; current Te-p model is a useful existence example, not a material prediction.
3. **Priority/source:** obtain original 2026 LEE SM (rather than infer from annotated reading notes) and compare the exact operator-dependent topology with Sec. I/III. Do not run an 80/528 census or magnetic domain-wall DFT until the specific effect survives.
4. **Computing permissions:** MnPS3 R1 remains CLOSED; R2-P0 READ_ONLY; no VASP/Slurm/WSL/Japan job requested or authorized.

**Outcome:** The mathematically precise \(w_{\pm}\pmod6\) constraint is genuinely further than the earlier first-order gradient idea. A material-specific microscopic six-vortex transition is **NOT** yet established; research integrity requires the cutoff and mesh corrections above.
