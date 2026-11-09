# MoSSe S/S: gate efficiency of net insulating NMEE

Use only the existing fully relaxed, unstrained 7269 zero-B full-KS bilayer and lower-monolayer references from mosse-strain-20261108/e000. No new DFT, strain, structure, or material.

Observable: C_net(F), both layer C(F), gate increment C(F)-C(0), and absolute slope G=dC_net/dF. Do not divide by near-zero C_net(0) to claim enhancement. Fixed positive probe Ex Bz. Model field grid: 0, ±0.025, ±0.05, ±0.1, ±0.25, ±0.5 V/nm; accept only insulating points with gap >0.05 eV.

Testable candidate, not a result: cross-layer occupied-empty mixing increases the gate-induced net response relative to matched isolated layers and a smooth electrostatic interface reference. Bringing a selected cross-layer occupied/empty pair closer may enhance its admixture and gate sensitivity, but competing channels could prevent enhancement. If the full response slope is not larger/robust, or a static reference reproduces the change, reject the enhancement claim. A different sign or a nonzero baseline alone is insufficient.

Apply a physical-space linear one-electron potential in all systems. This is F_model, an imposed potential slope, NOT a newly screened external gate voltage or a device efficiency. A bare layer offset is only a diagnostic, not an actual gate. New independent self-consistent DFT is deferred until this model gives a clear useful trend.

The matched isolated lower layer uses exact bilayer intralayer coordinates. The upper reference is its horizontal mirror; actual upper/lower relaxed intralayer heights are checked to differ by less than 0.001 A and the measured mismatch is reported. No claim of an independently calculated exact upper monolayer. Net isolated reference adds lower(F) and minus lower(-F). Static reference constructs separate lower and reflected upper interface potentials, using two spatial cutoffs. Both references include direct intralayer Stark response. This matters in an insulator too.

Full KS response space retained; selected layer valence/empty projectors are fixed at zero gate and zero probe B. Channel deletions act on H0+F V consistently, including selected gate off-diagonal blocks. They are nonadditive Hamiltonian interventions, not unique response contributions. Check 96/112 full states, three selected windows, 48/56/64 monolayer states, two electrostatic cutoffs, probe-step halving and the invariance under a rigid isolated-layer energy shift. Small total signals require absolute rather than relative checks.

All Japan-side operator exports and model calculations run in PBS compute nodes, 4 CPU/16GB/45min; no VASP invocation. Local-only input preparation and later small-data summaries/plots. Existing large wavefunctions are read in place, not copied. Output expected below 100 MB. Only known small files are transported through the frontend.
