# Minimal KM-derived valence/two-empty-sector model

An unfitted 12-state effective-model hypothesis test; not a MoSSe Wannier fit, material prediction, or new DFT calculation. All work in this checkpoint ran locally on the Mac. No cluster jobs or strain scans were launched.

## Hamiltonian and scope

Keep each layer's original four-state KM Hamiltonian, positions and physical Pauli spin operators, and add one spinful orbital at sublattice A (B is a counter-control) with onsite energy Delta. The full bilayer has 12 states, four occupied. The original KM interlayer block stays fixed at 10 meV. The added orbital's interlayer hopping envelope is copied from the original microscopic displacement-dependent hopping. Select its valence coupling with weight 1 and its empty-sector coupling with relative amplitude rho, using isolated-layer zero-field spectral projectors. Derivatives of the projectors are included in the electric velocity. The added orbital has no intralayer hybridization or independent SOC.

At rho=1 this is the straightforward local-orbital extension and it does NOT reverse the chosen bottom response in the tested main case. At rho<1 it is an orbital-selective projected effective coupling, generally longer range than the original orbital hopping. Such selectivity is a physical hypothesis, not a consequence of layer symmetry, and has NOT yet been extracted from MoSSe. Do not describe this as proof that merely adding any empty band produces reversal.

Main parameters: quarter-path stacking, m=45 meV, original Rashba SOC=20 meV and intrinsic SOC=10 meV, Delta=120 meV, rho=0.1, original interlayer scale=10 meV. Probe is the fixed 60-degree mirror-normal in-plane direction, B along z. The new orbital is higher than the original low conduction sector near the valleys, not necessarily energy-ranked above it everywhere in the BZ. No orbital energies or amplitudes were fitted to MoSSe.

## Results

Fine-grid bottom C in 1e-5 (nm V T)^-1, at new coupling g=0,40,80,120 meV: -3.146, -12.454, +0.996, +21.466. Gaps remain 13.256,15.213,17.623,20.340 meV. A full-Hamiltonian VC-only control also reverses (+31.907). The effect survives rho=0.25 (+6.977), but disappears at rho=1 (-75.444), at Delta=500 meV with the same g (-14.141), and for the B-site coupling at the same energies/strengths (-41.804). Removing Rashba makes this spin response vanish. Neighboring stackings 0.20/0.30 and m=46 meV preserve reversal; m=44 meV already has positive baseline C, so there is no sign reversal there.

The weak-g=40 versus strong-g=120 meV spectral controls keep lower C signs -,-,+,+ for weak states/weak energies, weak states/strong energies, strong states/weak energies, strong states/strong energies. This reproduces the logic of the material's state-driven reversal (whose absolute signs are opposite), without requiring the same numerical response scale. It is not an additive energy/geometry decomposition. Uniformly flattening the whole occupied/empty spectrum gives positive C at both endpoints, so no spectrum-independent sign theorem follows.

## Conditional prediction

The virtual-hybridization scale g^2/Delta suggests a leading threshold trend. Three reference roots at Delta=120,250,500 meV set one median scale, 0.05444 eV. Before evaluating Delta=180/350 meV, predictions were saved to main-thresholds-predictions.json. Predicted roots 98.99/138.03 meV compare with 95.10/145.43 meV calculated roots (errors +4.09%/-5.09%). Finer integration independently preserves opposite signs at 0.9 and 1.1 times each prediction.

This is interpolation within a fixed effective-model family, not an exact perturbative law, independent material prediction, or universal threshold. Actual denominators are Delta-e_v(k), while orbital and spin matrix elements and competing channels also matter. The three reference g_c^2/Delta values themselves vary about 0.0511–0.0673 eV. The square-root guide is a useful conditional scale, not a fitted universal phase boundary.

## Reproduce and audit

Maintained implementation: bsf_nmee/remote_empty_km.py. Main experiment: scripts/test_remote_empty_mechanism.py with --grid 41 40 64 (all controls), or --grid 81 80 120 and --select (representative refinement). Prediction experiment: scripts/predict_remote_empty_threshold.py; the prediction manifest is written before checking either held-out energy. Exploratory scripts/data are preserved under exploration; their old isolated-new-channel thresholds differ from the main model that retains the full old KM interlayer block. Do not mix them.

- controlled.json: complete principal controls at 20301 k points. The all_paths_CC control is metallic (negative indirect gap); its fixed-four-band numerical response is INVALID, even though raw C is retained for audit. Never cite that value as an enhancement.
- fine.json: representative controls at 76899 k points; g=0 reproduces the established original KM result.
- half-step.json: magnetic derivative step halved from 0.002 to 0.001 T.
- spectral-refined.json: complete weak/strong spectral controls; velocities differentiated from each control Hamiltonian and checked with two differences. Exact degeneracies cannot receive unequal target eigenvalues.
- main-thresholds*.json and heldout-fine.json: predeclared model predictions and checks.
- figure_source.csv, plot_model.py, minimal_model_mechanism.{png,pdf,svg}: discussion figure, not final submission.
- evidence.json: hashes, tests and numerical checks.

Tests cover the original KM limit, time reversal, the physical horizontal-mirror partner, C3 covariance under simultaneous stacking rotation, Hermiticity through construction, exact projector velocities, and flat-spectrum velocities. No proprietary data or credentials are included.

The next physical task is to derive the required coupling selectivity from MoSSe's valence/two-empty-sector projection, or show it fails, before adding more parameters or asserting a strain/stacking switch. The model demonstrates possibility beyond one chemical Hamiltonian; it does not establish universality across lattice classes or experimental accessibility.
