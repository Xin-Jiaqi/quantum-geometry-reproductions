# Interlayer valence/empty mixing controls: MoSSe and KM

Research mechanism checkpoint, not a new self-consistent DFT prediction.

## Definitions

At zero field diagonalize the two isolated intralayer blocks of the existing Hamiltonian. Their occupied projectors form Pv; Pc=1-Pv. For the interlayer block T, define VV=Pv T Pv, CC=Pc T Pc and VC=Pv T Pc+Pc T Pv. Freeze these projectors when applying the Zeeman probe. These are Hamiltonian interventions; the resulting responses are NOT additive contributions. The subspaces are independent of rotations within each complete valence/empty sector, but remain conditional on the chosen material projection and layer partition.

Material: 14 occupied and 8 empty states per layer in the 44-state d+p model. Embed each intervention into the original zero-field 96/112-state KS response space, retaining physical spin and Bloch-overlap operators. The weak baseline removes only this selected d+p interlayer block; it is not a separately relaxed or self-consistent isolated bilayer. Units are 1e-9 (nm V T)^-1; b in b sigma_z is in eV. Source cache: ../zero-field-refined-20261107/full-controls-cache.npz in the remote material project. No finite-B reference is used.

KM: existing quarter-path, m=45 meV case, probe along the 60-degree mirror-normal direction. Two occupied states per layer. Modified electric velocities include the derivatives of the zero-field projectors. The eight-state model was not extended or fitted. KM units differ from the material table; raw JSON uses (nm V T)^-1.

## Findings

At window8/bands96, MoSSe bottom C is +13.498 without selected interlayer mixing, +14.305 with VV only, +10.863 with CC only, -9.093 with VC only, +12.941 with VV+CC, and -6.602 with all channels. All six window/cutoff combinations retain these signs. The VC-only gap remains about 1.399 eV, versus 1.398 eV without mixing. Conversely VV alone reduces the gap to about 0.842 eV without reversing C. A smaller gap is therefore neither sufficient nor required for this intervention-induced reversal.

The valence side can be restricted to the top two eigenstates of each isolated layer block and reversal remains. These are local layer eigenstates 13–14 (one-based, counting 22 d+p spin states), NOT fixed band indices of the complete bilayer. They need not be degenerate at general k; the executed detail script's phrase “Kramers pair” refers to the two-state sector, not a claim of pointwise degeneracy. The minimum sampled separation from deeper valence states is about 0.155 eV.

Further split the empty sector into the two lowest layer eigenstates (15–16) and the other six (17–22). Top-valence to lowest-conduction mixing alone leaves C positive (+2.538); top-valence to higher-conduction mixing alone makes C negative (-9.010). Both signs persist in all six representations. In the full background, removing either of these top-valence cross-gap channels restores positive C (+13.005 or +13.428). Thus the higher empty sector matters, but is not a unique additive negative contribution; channel cooperation is essential to the full result. The smallest sampled separation between the two conduction sectors is 0.032 eV. These subdivisions have not been checked on a denser zero-field k mesh.

KM differs microscopically: at t=20 meV, bottom C/1e-5 is -7.413 without interlayer mixing, +35.522 with VV only, +11.229 with CC only, -47.878 with VC only, and +1.077 with all channels. The original full-model endpoints -3.146 (10 meV) and +1.077 (20 meV) are reproduced. The isolated VC channel reinforces the original sign. Do not claim that the existing KM parameterization already reproduces MoSSe's cross-gap-driven reversal.

## Minimal-model implication and limits

There is no evidence here that adding a second deep valence sector is necessary. A reasonable next candidate has one spinful valence sector and two distinguishable empty sectors per layer, with symmetry-compatible, orbital-selective interlayer coupling. If all sectors are represented by spin pairs, this is a 12-state bilayer candidate, NOT a demonstrated minimal model. Alternatively higher empty states might be downfolded into effective couplings and transformed spin/electric operators; deleting them while leaving the operators unchanged is not equivalent. The channel controls do not yet establish a specific atomic-orbital minimal basis or a predictive threshold for experimental strain/gating.

The existing 6x6 zero-field cache and its finite empty-state limitations remain. 96/112-state absolute amplitudes differ; 112-state B-even residuals remain. All present interventions remain insulating and the +/-b half-step checks preserve signs. Window and band-cutoff agreement supports this causal classification, not final quantitative convergence. No strain, new stacking, or new gate DFT was run. Future symmetry/strain predictions must account for changes in both orbital mixing and energy denominators, and be verified independently.

## Files and replay

- material-channels.json: all eight VV/CC/VC subsets.
- material-channel-detail.json: VC-strength controls and top/deeper valence separation.
- material-edge-channels.json: lowest/higher conduction controls.
- km-coarse.json, km-fine.json: identical channel interventions in the current KM model.
- channel_summary.csv: compact main controls; each system's units are explicit.
- evidence.json: numerical checks, hashes and finished PBS jobs (548191, 548192, 548193 on csuper001; exits 0).
- scripts/: executed material analysis and PBS inputs, plus the local model dependency snapshot. The legacy model.py comment about removing a B-prepared seed does not describe this run's source: the supplied cache is genuinely zero-field. Transport helpers, credentials, PAW data and large caches are excluded.

Material scripts expect the established remote parent model.py and zero-field-refined-20261107 cache; place mixing_channels.py in their working directory and submit the PBS file. This archive is not a standalone distribution of the licensed DFT inputs. The maintained KM entry point is studies/bsf-nmee-connection/scripts/compare_mixing_channels.py; run with --grid 81 80 120 --out OUTPUT.json. Shared projectors and derivatives live in bsf_nmee/mixing_channels.py. Tests check within-layer basis covariance, exact Hamiltonian reconstruction and analytic derivatives against finite differences.
