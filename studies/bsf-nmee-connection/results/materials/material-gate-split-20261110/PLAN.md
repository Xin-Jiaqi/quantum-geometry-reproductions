# Screened material gate-response split, before results

2026-11-10. Reuse the existing fixed-ion zero and ±0.10 V/nm SCF states of
MoSSe S/S and InSe A; no new DFT, material, structure, strain, or model state.
Question: is the small lower-layer gate change in MoSSe intrinsically small,
or the remainder of larger transition-energy and matrix-element changes?
Contrast InSe without forcing a single cross-layer VC interpretation.

Use the SAME squared-energy-denominator convention as the preceding KM test.
The archived covariant electric dipole r_cv defines the corresponding electric
vertex V_cv=-i(E_v-E_c)r_cv. Then the existing material kernel 2Re(S_vc r_cv)/D
is identically -2Im(S_vc V_cv)/D^2. This algebra fixes the split convention;
using a single-denominator split instead would move a term between parts.
The matrix part includes ALL changes of spin/electric vertices, including those
from independent screened gate states. No frozen-operator approximation.

At each magnetic probe sign, apply the exact symmetric product identity
between the independent positive/negative gate states, then take the existing
B-odd difference. The two parts must reproduce the archived net and layer gate
slopes. Retain MoSSe96/112 and InSe128/144/160 cutoffs, B steps0.002/0.001 eV.
Track bands via existing PAW cross-field overlaps, separately within occupied
and empty spaces. Compare zero-field-anchored tracking, direct tracking and
energy ordering; monitor overlaps and near degeneracies. If the attribution
is matching-sensitive, report that instead of claiming a unique mechanism.

Local synthetic/KM identity and phase/permutation checks precede submission.
Remote data loading, eigensystems and matrix processing run only in PBS compute
nodes, one4CPU/16GB allocation, maximum30min initially. Extra output <0.1GB;
no wavefunctions copied and no other shared-account job changed. Japan is
appropriate because both existing caches and validated PAW overlaps are there;
no new resource migration. Quota/queue preflight recorded separately.

Before submission: also retain the algebraically equivalent single-denominator
partition as a convention-sensitivity diagnostic. A velocity matrix reconstructed
from Delta*r includes changes of Delta itself; its derivative cannot be called
pure electronic-state mixing. Only the squared-denominator partition is directly
compared to the KM convention, and neither is an independently observable causal
partition. A claimed large internal cancellation must be assessed against this
convention dependence as well as state matching.
