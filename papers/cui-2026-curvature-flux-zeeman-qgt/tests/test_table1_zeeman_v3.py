"""Round-36: independent Table-1 audit and response oracle for the v3
Zeeman-QGT classification (Issue #105).

Covers the machinery in flux_lab.zeeman_v3:
  * sector_transformation_rule  - clean +1 covariant rule for the unitary
    generators; NO clean +-1 rule for the antiunitary generators (mixing),
    independently confirming the Round-14 conditional-divergence finding.
  * table1_independent_audit    - 31-group unitary-part projection with
    machine-readable disagreement reasons (antiunitary / 2band-class / etc.).
  * response_symmetry_oracle    - microscopic supp:eq:sigma channels for the
    representative 2-band models, incl. the recorded m' gN discrepancy and
    the Cn>=3 / m'm2' class limitations.
"""
import unittest
import numpy as np
from flux_lab.zeeman_v3 import (
    sector_transformation_rule, table1_independent_audit,
    response_symmetry_oracle, MODELS, _probe_models, _c3_d, _c3_dk,
    _sectors_at_point, generator_2d, spinor_rep,
)
from flux_lab.magnetic_symmetry import load_groups


class SectorTransformationRule(unittest.TestCase):
    def test_unitary_generators_clean_plus_one(self):
        rule = sector_transformation_rule()
        for g in ('C2', 'C3', 'Mx'):
            d = rule[g]
            for sec, sgn in d['signs'].items():
                self.assertEqual(sgn, 1, g)
                self.assertLess(d['residuals'][sec], 1e-9, g)

    def test_antiunitary_generators_have_no_clean_sign_rule(self):
        # the antiunitary ops mix the sectors: neither +1 nor -1 fits
        rule = sector_transformation_rule()
        for g in ("C2'", "Mx'"):
            d = rule[g]
            self.assertEqual(d['signs'], {k: 0 for k in d['signs']}, g)
            for sec, res in d['residuals'].items():
                self.assertGreater(res, 1e-3, (g, sec))

    def test_c3_model_is_invariant(self):
        k = np.array([0.37, -0.52])
        R = generator_2d('C3')
        R3 = np.block([[R, np.zeros((2, 1))], [np.zeros((1, 2)), 1.]])
        d1 = np.asarray(MODELS['c3'][0](k[0], k[1], 0.3))
        d2 = np.asarray(MODELS['c3'][0]((R @ k)[0], (R @ k)[1], 0.3))
        self.assertLess(np.linalg.norm(d2 - R3 @ d1), 1e-12)
        # velocity covariance under the chain rule
        dk1 = np.asarray(MODELS['c3'][1](k[0], k[1], 0.3))
        dk2 = np.asarray(MODELS['c3'][1]((R @ k)[0], (R @ k)[1], 0.3))
        self.assertLess(np.max(np.abs(dk2 - R @ dk1 @ R3.T)), 1e-12)


class Table1Audit(unittest.TestCase):
    def setUp(self):
        self.tab = {r['group']: r for r in table1_independent_audit()}

    def test_all_31_groups_covered(self):
        self.assertEqual(len(self.tab), 31)

    def test_unitary_only_groups_have_no_disagreement(self):
        # groups without antiunitary ops reproduce the paper rows exactly
        for name, row in self.tab.items():
            if row['n_antiunitary'] == 0:
                self.assertEqual(row['disagreements'], [], name)

    def test_gray_and_t_families_flagged_antiunitary(self):
        # groups with antiunitary ops cannot be decided by a clean sector rule
        for name in ("11'", "m1'", "21'", "mm21'", "31'", "41'", "61'"):
            row = self.tab[name]
            self.assertTrue(row['n_antiunitary'] > 0, name)
            self.assertTrue(row['disagreements'], name)

    def test_paper_rows_match_source_text(self):
        # spot-check the printed Table-I rows (transcribed from dualQGT7.tex)
        for name in ('1', 'm', '2', "2'", 'mm2', "m'm2'", '3', '4', '6'):
            row = self.tab[name]
            for sec in ('gN', 'gA', 'omegaN', 'omegaA'):
                self.assertTrue(row['paper'][sec], (name, sec))


class ResponseOracle(unittest.TestCase):
    def setUp(self):
        self.oracle = response_symmetry_oracle(0.3, Lam=6., N=160, M=320)

    def test_blue_all_four_groups_realize_all_channels(self):
        for name in ('1', '2', 'm', "2'", 'mm2'):
            self.assertEqual(
                self.oracle[name]['realized'],
                {'omegaA': True, 'gN': True, 'omegaN': True, 'gA': True},
                name)

    def test_m_prime_gN_contradicts_paper_row(self):
        # recorded discrepancy: Table-I gN row omits m' but the microscopic
        # m'-invariant 2-band model realizes a nonzero gN channel
        o = self.oracle["m'"]
        self.assertTrue(o['realized']['gN'])
        self.assertIn('gN', o['disagreements'])
        self.assertEqual(o['disagreements']['gN'], 'realized-but-not-listed')

    def test_c3_omegaN_gA_class_zero(self):
        # Cn>=3 radial-in-plane theorem: <Omega^N_z>, <g^A_xy> identically
        # zero in the 2-band class; Omega^A and gN stay nonzero
        o = self.oracle['3']
        self.assertFalse(o['realized']['omegaN'])
        self.assertFalse(o['realized']['gA'])
        self.assertTrue(o['realized']['omegaA'])
        self.assertTrue(o['realized']['gN'])

    def test_m_m2prime_dz_zero_class(self):
        # m'm2' forces d_z = 0, hence Omega^A and gN vanish identically
        o = self.oracle["m'm2'"]
        self.assertFalse(o['realized']['omegaA'])
        self.assertFalse(o['realized']['gN'])
        self.assertTrue(o['realized']['omegaN'])
        self.assertTrue(o['realized']['gA'])


class CnRadialTheorem(unittest.TestCase):
    def test_radial_inplane_integral_vanishes(self):
        # for any C3-symmetric 2-band model the d_x+i d_y part is z f(k^2)
        # (covariance under the C3 rotation forces this), which makes the
        # epsilon^2-weighted <Omega^N_z> and <g^A_xy> integrals vanish.
        from flux_lab.zeeman_v3 import _spinor_fields, polar_grid
        kx, ky, dA = polar_grid(8., 240, 480)
        w = dA / (2 * np.pi) ** 2
        d, E, dh, r, s = _spinor_fields(kx, ky, MODELS['c3'], 0.3)
        x = np.einsum('a...,b...->ab...', r, s)
        OmN = -np.imag(x[0, 1] - x[1, 0])
        gA = (x[0, 1].imag + x[1, 0].imag) / 2
        eps = 2 * E
        self.assertLess(abs(float(np.sum(OmN / eps ** 2 * w))), 1e-7)
        self.assertLess(abs(float(np.sum(gA / eps ** 2 * w))), 1e-7)


if __name__ == '__main__':
    unittest.main()
