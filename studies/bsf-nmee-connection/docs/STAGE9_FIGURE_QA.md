# Stage 9 figure contract and evidence

Core claim: full polar switching and local spin control are different questions; in this model the spin-response anisotropy contains information beyond the polarization sign.

Archetype: quantitative triptych, with full-path panel (a) setting the physical baseline and panels (b,c) exposing distinct controls. Python/matplotlib; editable SVG and PDF, 600-dpi PNG; 180 mm wide. This is a candidate subsection figure, not a declared PRL-ready novelty package.

Panel (a): all 17 t⊥=5 meV full-path points. Plot Pz/Pz(AB) and Cxx/Cxx(AB), not their raw units on a shared scale. Normalization compares signs and shapes only. The t⊥=20 meV path is retained as a robustness dataset in the source JSON, not hidden or fitted away.

Panel (b): all seven local-AB points for direction (0.23,0.41); plot percent dipole change and percent directional difference relative to the AB trace magnitude. Horizontal coordinate is actual displacement length with a signed direction, in pm. The second direction is an independent check retained in the JSON. No fits or invented dense samples.

Panel (c): fixed-stacking coupling values, λ=1/4. Left axis is the spin tilt from −x for E along +x and Bz positive, φ=arctan(Cyx/Cxx), since Cxx remains negative. All nine nonzero-hopping values are shown. At zero hopping the total spin vanishes and its angle is undefined; that one angular point is excluded for this physical reason. The right axis shows the electronic dipole for all ten hopping values. This dual axis compares orientation and polarization, not their amplitude or zero positions. Coupling is not relabeled as pressure.

Statistical n: not applicable; deterministic Brillouin-zone integrations. No statistical errors, p-values or biological/technical repeats. Numerical refinement is quantified separately. Source CSV includes every plotted point and its normalization; JSON retains all tensor components, gaps and parameters. Straight connecting lines are guides to the eye.

Integrity: no image processing or simulation-to-experiment substitution. Tiny floating-point residuals at symmetry zeros are plotted as computed. Manuscript field geometry is fixed Eparallel and Bz; structure is externally specified. Plot is not evidence that those fields cause sliding.

Caption: (a) Normalized electronic dipole and Cxx along AB→BA at t⊥=5 meV. (b) Dipole change and spin-response anisotropy under small shear from AB. (c) Induced spin tilt from −x and electronic dipole versus interlayer hopping at fixed λ=1/4, with E along +x and Bz positive. The spin tilts to opposite sides of −x while the dipole remains positive. Hopping amplitudes are theoretical controls. Source data and computational definitions accompany the figure.

QA: source preflight has 13 passes, zero failures and one reviewed warning (PNG preview instead of TIFF; vector PDF/SVG are the intended line-figure deliverables). The rendered panels were inspected for clipping, legend readability and correspondence to source values. The initial automated validator could not infer loop-based export suffixes or arithmetic figure dimensions; exports and dimensions are now explicit. No numerical data were changed for layout.
