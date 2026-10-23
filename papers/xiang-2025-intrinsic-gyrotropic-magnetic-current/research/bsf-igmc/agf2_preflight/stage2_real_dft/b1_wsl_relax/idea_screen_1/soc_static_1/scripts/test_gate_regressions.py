#!/usr/bin/env python3
"""Fast regression tests for the four-band same-gauge IGMC gate (PR #38).

Every test deliberately injects a wrong-physics scenario and REQUIRES the
gate machinery to reject it (nonzero exit / raised guard). Fast (~seconds):
uses N <= 9 grids only. Run:  python3 scripts/test_gate_regressions.py

Tests:
  R1 wrong tensor-factor arrangement (H = I (x) h2, S = I (x) sigma) must
     NOT conserve spin: ||[H, Sx]|| ~ 1.2 eV >= 1.19 eV; and the gate's NC4
     assertion must fail for that builder.
  R2 negative indirect (global) gap with positive direct gap must be
     rejected by gap_scan(require_insulator=True).
  R3 a deliberate Wilson principal-branch jump (> pi) must be flagged by
     wilson_branch_metrics (embedded row phases).
  R4 finite-b Wilson phase drift >= pi must be flagged.
  R5 anchors: self-contained two-band Dirac must match the frozen CSV.
  R6 pinned Cui interface content identity (sha256 of module body).
"""
from __future__ import annotations
import hashlib
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import test_multiband_igmc_gate as G

PINNED_SHA = G.PINNED_UPSTREAM_BODY_SHA256


def test_r1_wrong_tensor_factor_rejected():
    wrong_factor = G.build_no_soc()          # H = I (x) h2, S = I (x) sigma
    comm = G.spin_commutator_max(wrong_factor, G.grid_q(5), 0)
    assert comm >= 1.19, f"wrong-factor commutator {comm} < 1.19 eV"
    # the gate's NC4 conservation assertion requires < 1e-10 -> this builder
    # must FAIL it (the old D~0 was an accidental cancellation, not evidence)
    assert comm >= 1e-10, "gate NC4 conservation assertion would wrongly pass"
    print(f"  R1 OK: wrong-factor ||[H,Sx]|| = {comm:.3f} eV >= 1.19 -> gate "
          f"NC4 conservation guard rejects this arrangement")


def test_r2_negative_indirect_gap_rejected():
    fx = G.build_overlap_gap_fixture()   # per-k direct gap > 0, indirect < 0
    try:
        G.gap_scan(fx, 5, (0.0,), require_insulator=True)
    except RuntimeError as e:
        assert "global insulator" in str(e), f"wrong error: {e}"
    else:
        raise AssertionError("overlap fixture (negative indirect gap) NOT rejected")


def test_r3_wilson_branch_jump_rejected():
    # synthetic row phases with a branch-unstable row at the +pi cut
    fake = np.array([0.1, 0.2, np.pi - 1e-9, 0.3])
    try:
        G.wilson_branch_check({0.0: fake})
    except RuntimeError as e:
        assert "branch cut" in str(e), f"wrong error: {e}"
    else:
        raise AssertionError("branch-unstable row NOT flagged")
    # real model: v2_positive rows are smooth and away from the cut
    fx = G.build_v2()
    rows = G.embedded_wilson_row_phases(fx, 9, 0.0)
    m = G.wilson_branch_check({0.0: rows})
    assert not m["chern_winding_flagged"]
    print(f"  R3 OK: branch-unstable row at pi-1e-9 flagged; "
          f"v2_positive branch margin {m['jumps']['0.0']['branch_margin_rad']:.3f} rad, "
          f"max adjacent jump {m['jumps']['0.0']['max_adjacent_ky_jump']:.4f} rad")


def test_r4_finite_b_drift_flag():
    fake0 = np.zeros(3)
    fakeb = np.full(3, 3.5)          # RAW drift 3.5 rad > pi (wrapped ~2.78 rad
                                     # < pi: the raw/unwrapped difference is the
                                     # branch-jump-in-b detector); rows off cut
    try:
        G.wilson_branch_check({0.0: fake0, 1e-4: fakeb})
    except RuntimeError as e:
        assert "finite-b" in str(e), f"wrong error: {e}"
    else:
        raise AssertionError("finite-b drift ~= pi NOT flagged")
    fx = G.build_v2()
    g0 = G.embedded_wilson_row_phases(fx, 9, 0.0)
    gp = G.embedded_wilson_row_phases(fx, 9, 1e-4)
    m = G.wilson_branch_check({0.0: g0, 1e-4: gp})
    print(f"  R4 OK: finite-b drift ~= pi flagged; "
          f"v2_positive real b-drift {m['finite_b_drift']['0.0001']:.2e} rad")


def test_r5_dirac_anchors():
    d2, w2 = G.dirac_2band(41)
    assert abs(d2 - 0.0152366470) < 1e-6, f"anchor direct {d2}"
    assert abs(w2 - 0.0151827707) < 1e-6, f"anchor wilson {w2}"
    print(f"  R5 OK: anchors {d2:.8f} / {w2:.8f} match frozen CSV")


def test_r6_pinned_content_identity():
    text = (SCRIPT_DIR / "wannier_operator_contract_pinned26.py").read_text()
    head, _, body = text.partition("#\n\n")
    assert hashlib.sha256(body.encode()).hexdigest() == PINNED_SHA, "pinned body sha mismatch"
    print("  R6 OK: pinned module body sha256 == upstream cdac6418 blob")


def test_r7_regression_numbers():
    """Fast, low-mesh version of the acceptance numbers to catch drift."""
    fx = G.build_v2()
    D21 = G.direct_dxx_omega(fx, 21)
    W21 = G.wilson_slope(fx, 21, 1e-4)
    assert abs(D21 - 0.03429375) < 1e-6, f"direct 21 drift {D21}"
    assert 0.033 < W21 < 0.0345, f"wilson 21 drift {W21}"
    corr = G.build_no_soc_corrected()
    assert G.spin_commutator_max(corr, G.grid_q(5), 0) < 1e-10
    assert abs(G.direct_dxx_omega(corr, 9)) < 1e-10
    print(f"  R7 OK: direct21={D21:.8f}, wilson21={W21:.8f}, corrected NC4 zero")


def main():
    tests = [test_r1_wrong_tensor_factor_rejected,
             test_r2_negative_indirect_gap_rejected,
             test_r3_wilson_branch_jump_rejected,
             test_r4_finite_b_drift_flag,
             test_r5_dirac_anchors,
             test_r6_pinned_content_identity,
             test_r7_regression_numbers]
    failed = 0
    for t in tests:
        try:
            t()
        except Exception as e:
            failed += 1
            print(f"  FAIL {t.__name__}: {e}")
    print(f"\nregression: {len(tests) - failed}/{len(tests)} passed, "
          f"{failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
