# Stage 8 figure contract and QA

Claim: this mirror-paired insulating bilayer has a locally cubic electronic
layer dipole and trace of the static EB spin response, but a linear anisotropic
response. Interlayer hopping changes electronic states as well as energy weights.
This is a model result, not an experimental sensitivity or PRL novelty claim.

Backend: Python/matplotlib. Archetype: quantitative four-panel evidence figure,
183 mm wide. Source: `results/bilayer_stage8_audit.json`; every point is a computed
Hamiltonian response. Connected lines only guide the eye. No simulated stand-in
data, fitted curves, smoothing, error bars, statistical p values or cropped data.

- a: three observables vs small slide, each normalized to its own absolute value
  at u=0.1. Omit u=0 because a logarithm of zero is undefined. This compares powers,
  not physical magnitudes. All four computed nonzero points through u=0.2 shown.
- b: signed Cxx and Cyy through u=0.4, with the isotropic and anisotropic parts of
  Cxx shown separately. All use the same coefficient units; u is displacement
  amplitude along a fixed direction, not an angular variable.
- c: full, same-sublattice-only, cross-sublattice-only, and relative-sign-flipped
  hopping Hamiltonians at fixed u=0.2. Divide Cxx-Cyy by t_perp squared to test
  weak-coupling scaling. These are modified models, not unique additive channels.
- d: physical and spectrally flattened Cxx for t_perp=5 meV. Flat spectra are
  exactly +/-15 meV at all k and all slides; physical projectors and spin operators
  are retained at B=0. All five computed slides shown. Flattening is artificial.

Export: editable SVG, PDF, 600 dpi PNG; common sans-serif fonts, minimum 7 pt.
The JSON holds values and parameters; the plotting script exports a long-form
CSV of every plotted series with transformations and units given by panel axes
and the plotting script (columns are panel, series, x, y).
Convergence is in the main report, not represented as statistical uncertainty.
Final source preflight and PNG visual inspection are required before delivery.

QA completed: source preflight has zero failures. The one warning is the absence
of a TIFF submission raster: this research delivery intentionally uses PNG plus
editable SVG/PDF and is not a journal submission package. The rendered PNG was
visually inspected: labels, legends, signs and panel ordering agree with the data.
