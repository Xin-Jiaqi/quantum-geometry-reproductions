# Unified Physics Interface v1 — Protocol-first design

**Scope:** a versioned interface and metadata-only registry. **Not yet a cross-paper numerical engine.**

The existing four paper modules have different physics, units, magnetic group treatments, source conventions and validation ranges. The interface must enforce these distinctions rather than expose one misleading universal QGT number.

## Architecture

~~~text
paper PDF / Supplemental Material + paper README and acceptance logs
        |
        v
source formula + units + occupation + perturbation + response protocol
        |
HamiltonianProvider: H(k), dH/dk, basis, spin/position/layer observables
        |
protocol-aware adapter [FUTURE: one tested paper at a time]
        |
ResponseResult: value + tensor convention + units + source SHA + diagnostics
        |
optional independent group-theory symmetry and real-material workflow
~~~

**Now implemented:** integration/physics_interface/protocol_registry.json, registry.py (validated model/entrypoint discovery), contracts.py (typed design for BlochProvider/ResponseRequestV1/ResponseResultV1), cli.py (read-only inspection and planning), tests, GitHub CI.

**Not implemented:** a dispatcher that numerically evaluates four paper responses. Calling contracts.evaluate intentionally raises NotImplementedError. The registry's METADATA_VALID_ONLY means file/format checks, not scientific acceptance.

## 1. Input contract

A future numerical adapter accepts a provider that returns a finite Hermitian **H(k)** in a declared orbital basis with a stable basis_id, an explicitly defined k coordinate convention and energies. Also required according to the selected response:

- dH/dk in a consistent gauge, or a validated velocity/position/Berry-connection prescription if the Bloch basis itself depends on k.
- Physical spin S_i or Pauli sigma_i **with the physical hbar factor specified**, spin-current J^a_i, magnetization, and optional layer projectors. All operators are in the **same declared basis**.
- Magnetic moment / Zeeman coupling for B-driven calculations, spinful SOC matrices, and optional layer exchange for stacking.
- Lattice dimensions, momentum metric, Brillouin-zone measure, energy unit, 2D-sheet vs 3D-bulk output normalization.

Do not assume J^a_i = {v_a,sigma_i}/4 equals {v_a,S_i}/2 until the definition S=hbar*sigma/2 and the source paper's normalization have both been checked. Same-gauge spin and position operator data cannot be recovered automatically from wannier90_hr.dat alone.

## 2. ResponseProtocol contract

Every calculation binds these **independent** parameters: chosen paper method and stable protocol ID; observable and its tensor indices; perturbation type (E, B, E+B, or second-order E); frequencies and Fourier convention; order of perturbations/field permutation; occupation and chemical potential; temperature; eta interpretation and zero-frequency/zero-eta order of limits; degeneracy policy; magnetic state and external-field/response parity; SI or model units.

The specific stable IDs and existing real CLI entry points are listed in integration/physics_interface/protocol_registry.json, and may be inspected by:

~~~bash
python -m integration.physics_interface.cli list
python -m integration.physics_interface.cli inspect xiang2025-spin-qgt
python -m integration.physics_interface.cli plan xiang2025-spin-qgt spinqgt.finite_eta
python -m integration.physics_interface.cli validate
~~~

The **plan** output is JSON describing command argument vectors, working directory and caveats. It **never executes** commands or asserts scientific PASS. Run the paper's own verifier separately in an isolated environment.

### Do not conflate these protocols

| Stable protocol | It means | It does NOT automatically mean |
| --- | --- | --- |
| spinqgt.reduced_eq13 | two-band reduced paper formula | general metallic intrinsic tau^0 |
| spinqgt.finite_eta | full adiabatic finite-eta Liouville with Drude poles | a calibrated collision model |
| spinqgt.laurent_C | scheme-specific coefficient of eta^0 | physically unique driven DC intrinsic |
| spinqgt.sm_eq2 | explicit nondegenerate four-source identity | PT-degenerate band formula |
| nmee.printed | printed model expression | full equilibrium static response |
| nmee.full_static | independent static/equilibrium branch | paper finite-T driven DC |
| igmc.full_frequency | defined clean complex AC susceptibility | general static magnetoelectric transport |
| qgmap.two_color_eb | prescribed dynamic E/B ordering | steady-state EB PHE |

When definitions differ, return an explicit BLOCKED_PROTOCOL status or separate response results, **not** an averaged or force-normalized number. Unexpected degeneracies must be handled with verified subspace projectors or rejected.

## 3. Proposed response object (example only; no physics calculated)

~~~json
{
  "schema_version": "1.0.0",
  "source_git_sha": "FULL_COMMIT_SHA_RECORDED_BY_STUDY",
  "module_id": "xiang2025-spin-qgt",
  "protocol_id": "spinqgt.finite_eta",
  "observable": {"name": "layer_staggered_spin", "basis_id": "A"},
  "input_units": {"energy": "eV", "momentum": "1/angstrom"},
  "frequency": {"value": 0.0, "unit": "rad/s"},
  "broadening": {"eta": 0.01, "type": "adiabatic_regulator"},
  "result": null,
  "output_unit": "model",
  "scientific_status": "BLOCKED_PROTOCOL",
  "evidence": [],
  "convergence": {},
  "caveats": ["EXAMPLE ONLY; NOT A CALCULATED VALUE"]
}
~~~

Future adapters should include the full physical operator matrices' names and provenance, k meshes, dimensionless/SI factor breakdown, occupied-band convention, requested susceptibility tensor index order, degeneracy/gauge tests, comparison with the paper baseline, and all code/data hashes.

## 4. Symmetry-operator handoff

Use the [Group Theory Operations Toolkit](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit) APIs, not a copied magnetic-group table. Known v0.18.0 repository SHA: 551a0daf8d4ba5a5f24b2811694bfe438616de48.

For a spatially unitary operation R, a polar spin-current flow index a, an **axial spin** index i, and two polar electric fields b,c, the *spatial* response representation has the skeleton:

\[
\chi'^{a i}_{bc}
=\det(R) R_{aa'}R_{ii'}R_{bb'}R_{cc'}\,
\chi^{a'i'}_{b'c'}.
\]

This is **not yet** the full magnetic/symmetry rule. Add actual layer permutations for layer-resolved spin observables, and separately derive antiunitary time-reversal factors, complex conjugation, frequency exchange and protocol-dependent magnetic-domain preparation. A first-order Edelstein tensor is not the same representation as a second-order spin-current response.

For sliding/stacking, T(k, delta) should be built from physically displaced orbital hopping and embeddings, respecting covariance under the material's actual atomic motif. A global interlayer phase that only changes embedding is not automatically new stacking physics. Any intercalation materials claim requires electronic charge transfer and structural/energetic input.

## 5. Acceptance criteria for a future working adapter

1. No edits to frozen papers/ computation, source figure data or reference snapshots.
2. Reproduce each selected paper's own command and CI on the actual pinned SHA.
3. Show independent solver or formula agreement on the same H(k), same operator and the **same** physical protocol.
4. Validate Hermiticity, gauge covariance, PT/other degeneracies, layer exchange, time parity, expected zeroes and exact limiting cases.
5. Perform k, eta, temperature and frequency scans with explicit stopping tolerance, dimensional/SI audit, and printed source-equation comparison.
6. Freeze study-level results and source/module/toolkit SHAs; document differences as OPEN_DISCREPANCY, not hidden corrections.

**Current milestone:** versioned discover/plan/spec/test infrastructure. **Next milestone:** read-only thin adapters and their cross-paper physical regression tests, developed one protocol at a time under studies/ or a new scoped integration directory.
