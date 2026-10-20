"""Independent regression tests for the axial-vector magnetic connector.

PR #8 second-round review (BLOCKER 1): real spin is an AXIAL vector; under a
spatial operation R with time-reversal tau the moment transforms as

    m' = (-1)**tau * det(R) * (R @ m).

For the AgF2 connector m_001 (R = diag(1,1,-1), det(R) = -1) and collinear
moments along [100] this means:

    unitary   m_001  (tau=0):  S_x -> -S_x
    antiunitary m'_001 (tau=1): S_x -> +S_x

An earlier draft applied a polar-vector (z-flip-only) rule, giving exactly the
opposite mapping and false-positive PASSes. These tests pin the correct rule
(a non-collinear probe with S_x,S_y,S_z != 0), re-verify the B1->B2 connector
for both types with the paper's definitions (type-I = interlayer Neel vectors
OPPOSITE, unitary connector; type-II = Neel vectors SAME, antiunitary
connector), and assert the discriminating power (swapped assignments or the
old polar rule must FAIL).

Stdlib-only (no numpy/spglib/toolkit) so it runs on the plain CI runners of
bsf-igmc-agf2-stage6.yml (python -m unittest discover -s stage6_harness/tests).
"""
from __future__ import annotations

import json
import math
import pathlib
import unittest

BASE = pathlib.Path(__file__).resolve().parents[1]          # stage6_harness
REVIEW = BASE / "stage1_review"
B1_POSCAR = REVIEW / "structures" / "constructed_B1_candidate" / "POSCAR"
B2_POSCAR = REVIEW / "structures" / "constructed_B2_candidate" / "POSCAR"

TOL = 1e-6


def wrapf(v: float) -> float:
    w = v % 1.0
    return 0.0 if (w > 1.0 - 1e-6 or w < 1e-6) else w


def read_poscar(path: pathlib.Path):
    """Return (cell_diag, numbers, frac) for the orthorhombic POSCARs."""
    lines = path.read_text(encoding="utf-8").splitlines()
    cell = []
    for row in lines[2:5]:
        cell.append([float(x) for x in row.split()])
    # assert orthorhombic diagonal cell (a,b,c; all off-diagonal ~0)
    for i in range(3):
        for j in range(3):
            if i != j:
                assert abs(cell[i][j]) < 1e-6, "expected orthorhombic cell"
    diag = [cell[0][0], cell[1][1], cell[2][2]]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    n = sum(counts)
    coords = [[float(x) for x in lines[8 + i].split()] for i in range(n)]
    kind = lines[7].strip().lower()
    if kind.startswith("car"):
        frac = [[c[0] / diag[0], c[1] / diag[1], c[2] / diag[2]] for c in coords]
    else:
        frac = coords
    numbers = []
    for s, c in zip(species, counts):
        z = 47 if s == "Ag" else 9
        numbers += [z] * c
    return diag, numbers, frac


def assign_moments(frac, numbers, t_o, type_="I"):
    """Paper-consistent moments (see verify_magnetic_symmetry.py):
    bottom Ag(0,0)=+x, Ag(1/2,1/2)=-x; top images at t_o and (1/2+t_o).
    type-I: top images FLIP sign (Neel vectors OPPOSITE);
    type-II: top images KEEP sign (Neel vectors SAME)."""
    topA = (wrapf(t_o[0]), wrapf(t_o[1]))
    topB = (wrapf(0.5 + t_o[0]), wrapf(0.5 + t_o[1]))
    sign = -1.0 if type_ == "I" else +1.0
    out = []
    for i, n in enumerate(numbers):
        m = [0.0, 0.0, 0.0]
        if n == 47:
            x, y = wrapf(frac[i][0]), wrapf(frac[i][1])
            if abs(x) < 1e-4 and abs(y) < 1e-4:
                m[0] = 0.5
            elif abs(x - 0.5) < 1e-4 and abs(y - 0.5) < 1e-4:
                m[0] = -0.5
            elif abs(x - topA[0]) < 1e-4 and abs(y - topA[1]) < 1e-4:
                m[0] = 0.5 * sign
            elif abs(x - topB[0]) < 1e-4 and abs(y - topB[1]) < 1e-4:
                m[0] = -0.5 * sign
            else:
                raise AssertionError(f"unrecognized Ag site {x},{y}")
        out.append(m)
    return out


def axial_transform(m, R, det, tau):
    """m' = (-1)**tau * det(R) * (R @ m) -- stdlib implementation."""
    rm = [sum(R[i][j] * m[j] for j in range(3)) for i in range(3)]
    f = ((-1.0) ** tau) * det
    return [f * v for v in rm]


def connector_map(f1, m1, f2, m2, tau_frac, tau):
    """Apply {m_001|tau_frac} to B1 atoms+spins (axial rule) and compare with
    B2 modulo lattice. Returns (position_failures, moment_mismatches)."""
    R = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]
    det = -1.0
    pos_bad = 0
    mom_bad = 0
    for i in range(len(f1)):
        x, y, z = f1[i]
        xp, yp, zp = wrapf(x + tau_frac[0]), wrapf(y + tau_frac[1]), wrapf(2 * 0.5 - z)
        best, bd = None, 1e9
        for j in range(len(f2)):
            d = (abs(xp - wrapf(f2[j][0])) + abs(yp - wrapf(f2[j][1]))
                 + abs(zp - wrapf(f2[j][2])))
            if d < bd:
                bd, best = d, j
        if bd > 1e-5:
            pos_bad += 1
            continue
        mm = axial_transform(m1[i], R, det, tau)
        if any(abs(mm[k] - m2[best][k]) > 1e-9 for k in range(3)):
            mom_bad += 1
    return pos_bad, mom_bad


class TestAxialSpinRule(unittest.TestCase):
    """Independent non-collinear probe of the axial-vector rule."""

    def test_axial_rule_probe_all_components(self):
        m = [0.3, 0.7, -0.5]
        R = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]
        det = -1.0
        # unitary m_001 (tau=0): m' = -R m = (-0.3, -0.7, -0.5)
        self.assertTrue(all(
            abs(a - b) < TOL for a, b in
            zip(axial_transform(m, R, det, 0), [-0.3, -0.7, -0.5])))
        # antiunitary m'_001 (tau=1): m' = +R m = (0.3, 0.7, 0.5)
        self.assertTrue(all(
            abs(a - b) < TOL for a, b in
            zip(axial_transform(m, R, det, 1), [0.3, 0.7, 0.5])))
        # the polar-vector guess (flip z only) is NOT what spins do
        self.assertTrue(axial_transform(m, R, det, 0)[0] != m[0]
                        and axial_transform(m, R, det, 0)[1] != m[1])

    def test_unitary_m001_flips_inplane_axial(self):
        # S || x : unitary m_001 -> -S_x ; antiunitary -> +S_x
        self.assertEqual(axial_transform([1.0, 0.0, 0.0],
                                         [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
                                         -1.0, 0), [-1.0, 0.0, 0.0])
        self.assertEqual(axial_transform([1.0, 0.0, 0.0],
                                         [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
                                         -1.0, 1), [1.0, 0.0, 0.0])
        # out-of-plane S || z : unitary m_001 keeps S_z (and det flips z twice)
        self.assertEqual(axial_transform([0.0, 0.0, 1.0],
                                         [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
                                         -1.0, 0), [0.0, 0.0, 1.0])


class TestB1B2Connector(unittest.TestCase):
    """Full per-atom B1->B2 magnetic connection with the paper's definitions."""

    def setUp(self):
        self.d1, self.n1, self.f1 = read_poscar(B1_POSCAR)
        self.d2, self.n2, self.f2 = read_poscar(B2_POSCAR)
        self.tau_frac = (0.5, 0.16)

    def test_typeI_unitary_connector_exact(self):
        # type-I: Neel vectors OPPOSITE; unitary m_001 (tau=0) flips S_x
        m1 = assign_moments(self.f1, self.n1, (0.5, -0.16), "I")
        m2 = assign_moments(self.f2, self.n2, (0.5, 0.16), "I")
        pos, mom = connector_map(self.f1, m1, self.f2, m2, self.tau_frac, tau=0)
        self.assertEqual((pos, mom), (0, 0))
        # sanity: type-I top layer moments are opposite to the bottom ones
        self.assertEqual(m1[0], [-0.5, 0.0, 0.0])   # bottom Ag(1/2,1/2)
        self.assertEqual(m1[1], [0.5, 0.0, 0.0])    # bottom Ag(0,0)
        self.assertEqual(m1[2], [0.5, 0.0, 0.0])    # top image of B -> flipped
        self.assertEqual(m1[3], [-0.5, 0.0, 0.0])   # top image of A -> flipped

    def test_typeII_antiunitary_connector_exact(self):
        # type-II: Neel vectors SAME; antiunitary m'_001 (tau=1) keeps S_x
        m1 = assign_moments(self.f1, self.n1, (0.5, -0.16), "II")
        m2 = assign_moments(self.f2, self.n2, (0.5, 0.16), "II")
        pos, mom = connector_map(self.f1, m1, self.f2, m2, self.tau_frac, tau=1)
        self.assertEqual((pos, mom), (0, 0))

    def test_discriminating_power_swapped_tau_fails(self):
        # the connector tau must match the configuration: swapping must FAIL
        m1 = assign_moments(self.f1, self.n1, (0.5, -0.16), "I")
        m2 = assign_moments(self.f2, self.n2, (0.5, 0.16), "I")
        pos, mom = connector_map(self.f1, m1, self.f2, m2, self.tau_frac, tau=1)
        self.assertGreater(mom, 0)
        m1 = assign_moments(self.f1, self.n1, (0.5, -0.16), "II")
        m2 = assign_moments(self.f2, self.n2, (0.5, 0.16), "II")
        pos, mom = connector_map(self.f1, m1, self.f2, m2, self.tau_frac, tau=0)
        self.assertGreater(mom, 0)

    def test_old_polar_rule_is_rejected(self):
        # the historical buggy rule (m' = m with z flipped; -m for type II)
        # must NOT pass the corrected type-I check (it produced false positives)
        m1 = assign_moments(self.f1, self.n1, (0.5, -0.16), "I")
        m2 = assign_moments(self.f2, self.n2, (0.5, 0.16), "I")
        pos_bad, mom_bad = 0, 0
        for i in range(len(self.f1)):
            x, y, z = self.f1[i]
            xp, yp, zp = (wrapf(x + 0.5), wrapf(y + 0.16), wrapf(1.0 - z))
            best, bd = None, 1e9
            for j in range(len(self.f2)):
                d = (abs(xp - wrapf(self.f2[j][0])) + abs(yp - wrapf(self.f2[j][1]))
                     + abs(zp - wrapf(self.f2[j][2])))
                if d < bd:
                    bd, best = d, j
            mm = list(m1[i])
            mm[2] *= -1.0        # old polar-vector rule: flip z only
            if any(abs(mm[k] - m2[best][k]) > 1e-9 for k in range(3)):
                mom_bad += 1
        self.assertGreater(mom_bad, 0, "old polar rule must not pass type-I")


class TestReportConsistency(unittest.TestCase):
    """The generated report must reflect the corrected axial rule."""

    def test_report_uses_axial_rule_and_paper_labels(self):
        rep = json.loads((REVIEW / "results" / "magnetic_symmetry_report.json")
                         .read_text(encoding="utf-8"))
        for t in ("I", "II"):
            c = rep[f"connector_B1_to_B2_type{t}_magnetic"]
            self.assertEqual(c["moment_mismatches"], 0)
            self.assertTrue(c["positions_match_mod_lattice"])
            self.assertIn("det(R)", c["spin_transform_rule"])
            self.assertEqual(len(c["per_atom"]), 12)
        self.assertIn("OPPOSITE", rep["connector_summary"]["type_I"])
        self.assertIn("SAME", rep["connector_summary"]["type_II"])
        self.assertIn("m'_100", rep["B1_typeI"]["internal_magnetic_point_group"])


if __name__ == "__main__":
    unittest.main()
