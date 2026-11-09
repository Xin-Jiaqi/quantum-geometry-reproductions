"""Independent geometry, symmetry, and source-equation acceptance tests."""
import unittest
import numpy as np
from altermagnet_bpve import (Params, band_observables, contour_injection,
                            analytic_xxx_eq17, grid_injection,
                            analytic_xyy_printed_eq19)


class EzawaTests(unittest.TestCase):
    def test_metric_is_transition_matrix_element(self):
        sx = np.array([[0, 1], [1, 0]], complex)
        sy = np.array([[0, -1j], [1j, 0]], complex)
        sz = np.diag([1, -1]).astype(complex)
        p = Params(altermagnet=.1, gap=.1)
        for x, y in [(.19, .31), (-.27, .11), (.55, -.20)]:
            o = band_observables(x, y, p)
            hx, hy, hz = o['spin_field']
            H = hx*sx + hy*sy + hz*sz
            _, U = np.linalg.eigh(H)
            dHx = (p.rashba+2*p.altermagnet*x)*sy
            dHy = -p.rashba*sx-2*p.altermagnet*y*sy
            for a, D1 in enumerate((dHx, dHy)):
                for b, D2 in enumerate((dHx, dHy)):
                    # r_{-+}^a r_{+-}^b = v_{-+}^a v_{+-}^b / DeltaE^2
                    g = np.real((U[:, 0].conj() @ D1 @ U[:, 1])
                                *(U[:, 1].conj() @ D2 @ U[:, 0])) / o['gap']**2
                    self.assertAlmostEqual(float(o['metric'][a, b]), float(g), places=12)

    def test_eq17_weak_j_without_fitted_factor(self):
        for E in (.25, .4, .8, 1.2):
            p = Params(altermagnet=.005, gap=.1)
            exact = contour_injection(E, p, ntheta=720)
            approx = analytic_xxx_eq17(E, p)
            self.assertLess(abs(exact-approx)/abs(approx), 2e-5)

    def test_rashba_only_zero_and_neel_reversal(self):
        for E in (.25, .4, .8):
            self.assertLess(abs(contour_injection(E, Params(altermagnet=0.0))), 1e-12)
            plus = contour_injection(E, Params(altermagnet=+.10))
            minus = contour_injection(E, Params(altermagnet=-.10))
            self.assertLess(abs(plus+minus), 1e-11)

    def test_inversion_zero_and_threshold(self):
        self.assertEqual(contour_injection(.18, Params()), 0.0)
        with self.assertRaises(ValueError):
            contour_injection(.4, Params(gap=0.0))
        with self.assertRaises(ValueError):
            contour_injection(3.9, Params(), kmax=.7)

    def test_rashba_zero_inversion_forbids_dipole_bpve(self):
        p = Params(rashba=0.0, altermagnet=.1, gap=.1)
        self.assertLess(abs(grid_injection(.4, p, nk=301, gamma=.02)), 1e-13)

    def test_mu_pauli_blocking(self):
        p = Params(mu=-20)
        self.assertAlmostEqual(contour_injection(.4, p), 0, places=14)

    def test_convergence_and_independent_grid(self):
        p = Params(altermagnet=.03)
        E = .4
        sharp_1 = contour_injection(E, p, ntheta=180)
        sharp_2 = contour_injection(E, p, ntheta=720)
        self.assertLess(abs(sharp_1-sharp_2), 2e-8)
        grid = grid_injection(E, p, nk=651, kmax=1.2, gamma=.01)
        self.assertLess(abs(grid-sharp_2)/abs(sharp_2), .05)

    def test_xyy_source_discrepancy_recorded_not_fitted(self):
        p = Params(altermagnet=.005)
        E = .4
        numeric = contour_injection(E, p, optical='yy')
        printed = analytic_xyy_printed_eq19(E, p)
        expected_from_metric = -p.altermagnet/(8*p.rashba)*(
            1-8*(p.gap/E)**2+16*(p.gap/E)**4)
        self.assertLess(abs(numeric-expected_from_metric), 2e-8)
        self.assertGreater(abs(numeric-printed)/abs(printed), .30)


if __name__ == '__main__':
    unittest.main()
