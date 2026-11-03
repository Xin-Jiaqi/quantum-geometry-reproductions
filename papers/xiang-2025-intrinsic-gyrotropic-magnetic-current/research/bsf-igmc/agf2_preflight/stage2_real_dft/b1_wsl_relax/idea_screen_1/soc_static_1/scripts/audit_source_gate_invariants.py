#!/usr/bin/env python3
"""Source-integrity preflight for PR #38 synthetic four-band IGMC gate.

Do not run VASP/PBS/Wannier. It inspects existing synthetic fixtures and
records where 'direct gap' and 'true indirect gap' differ. It is deliberately
read-only and prints JSON. It may report FAIL for the pre-fix NC4 fixture.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_multiband_igmc_gate import build_no_soc, build_v2, grid_q
from wannier_operator_contract_pinned26 import WannierOperatorBundle, UNITS

SAMPLE_Q = ([0.0, 0.0], [0.173, 0.317], [0.419, 0.227])


def commutator_max(bundle):
    maximum = 0.
    for q in SAMPLE_Q:
        d = bundle.fourier(q)
        h, sx = d["H"], d["spin"][0]
        maximum = max(maximum, float(np.linalg.norm(h @ sx - sx @ h, 2)))
    return maximum


def corrected_no_soc_fixture():
    """Fix orbital first, spin second: H=h_orb x I_spin, S=I_orb x sigma."""
    old = build_no_soc()
    h = np.zeros_like(old.h_r)
    for i in range(len(old.translations)):
        orbital_h = old.h_r[i, :2, :2]
        h[i] = np.kron(orbital_h, np.eye(2, dtype=complex))
    return WannierOperatorBundle(
        old.translations.copy(), old.lattice.copy(), h,
        old.spin_r.copy(), old.position_r.copy(), old.centers.copy(),
        dict(old.provenance), dict(UNITS),
    ).validate()


def direct_indirect_gap(bundle, n=21, b=0.0):
    direct_min = float("inf")
    valence_max = -float("inf")
    conduction_min = float("inf")
    for q in grid_q(n):
        d = bundle.fourier(q)
        eig = np.linalg.eigvalsh(d["H"] - b * d["spin"][0])
        direct_min = min(direct_min, float(eig[2] - eig[1]))
        valence_max = max(valence_max, float(eig[1]))
        conduction_min = min(conduction_min, float(eig[2]))
    return {
        "direct_min_gap_eV": direct_min,
        "valence_max_eV": valence_max,
        "conduction_min_eV": conduction_min,
        "true_indirect_gap_eV": conduction_min - valence_max,
        "sampled_k_count": n*n,
    }


def main():
    original = build_no_soc()
    fixed = corrected_no_soc_fixture()
    bad_comm = commutator_max(original)
    good_comm = commutator_max(fixed)
    positive = build_v2()
    results = {
        "scope": "synthetic code integrity; not materials",
        "source_nc4_commutator_norm_eV": bad_comm,
        "corrected_nc4_commutator_norm_eV": good_comm,
        "source_nc4_is_actually_spin_conserving": bad_comm < 1e-10,
        "corrected_nc4_spin_conserving": good_comm < 1e-10,
        "positive_fixture_gap_checks": {
            str(b): direct_indirect_gap(positive,21,b)
            for b in (0., -1e-4, +1e-4, -1e-3, +1e-3)
        },
        "source_reported_NC2_signs": {
            "response_cov": -0.025357180978182244,
            "response_bare": +0.03429374753149531,
            "note": "Committed JSON/source labels; Markdown signs reversed.",
        },
    }
    print(json.dumps(results, indent=2))
    # This program is an independent preflight that highlights the old NC4 bug;
    # its own success means the diagnostic detected the error correctly.
    assert 1.19 < bad_comm < 1.21, "Original NC4 commutator unexpectedly changed"
    assert good_comm < 1e-10, "Corrected model must conserve spin"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
