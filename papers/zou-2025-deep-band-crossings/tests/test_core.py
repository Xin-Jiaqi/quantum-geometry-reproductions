"""Independent checks of SM (S6), (S10), (S12), (S19) and Table S1."""
import unittest
import numpy as np
from deep_dbc import (EPSILONS, SOURCE_E0, SOURCE_VRATIO, matrices,
                      nodal_parameters, band_data, covariant_derivative,
                      hermitian_connection, singular_connection, shg_pieces,
                      benchmark)


class DeepBandCrossingTests(unittest.TestCase):
    def test_source_hamiltonian_mirror_and_spinless_tr(self):
        x, y, ep = .32, 1.1, 1.03
        mx = np.diag([1., -1., 1.])
        h, v, w = matrices(x, y, ep)
        np.testing.assert_allclose(h, h.conj().T, atol=1e-14)
        np.testing.assert_allclose(matrices(-x, y, ep)[0], mx @ h @ mx, atol=1e-14)
        np.testing.assert_allclose(matrices(-x, -y, ep)[0], h.conj(), atol=1e-14)
        # Independent centered finite difference of analytical dH and d²H.
        delta = 1e-5
        for axis in (0, 1):
            dx, dy = (delta, 0) if axis == 0 else (0, delta)
            hp = matrices(x + dx, y + dy, ep)[0]
            hm = matrices(x - dx, y - dy, ep)[0]
            np.testing.assert_allclose((hp-hm)/(2*delta), v[axis], atol=2e-10)
            np.testing.assert_allclose((hp-2*h+hm)/delta**2,
                                       w[(0, 2)[axis]], atol=1e-5)

    def test_table_s1_all_six_parameters(self):
        for ep, e0ref, vref in zip(EPSILONS, SOURCE_E0, SOURCE_VRATIO):
            row = nodal_parameters(ep)
            self.assertLess(abs(row['epsilon0'] - e0ref), 0.001)
            self.assertLess(abs(row['vv_over_vc'] - vref), 0.001)
            energies = np.linalg.eigvalsh(matrices(0, row['ky'], ep)[0])
            self.assertLess(energies[2] - energies[1], 1e-12)

    def test_covariant_derivative_vs_parallel_transport_finite_difference(self):
        # Independently calculate r_nm(k+dk) with state phases aligned to k.
        x, y, ep, n, m = .29, 1.29, .9, 1, 0
        d = band_data(x, y, ep)
        H, _, _ = matrices(x, y, ep)
        _, u0 = np.linalg.eigh(H)
        dk = 2e-5
        for deriv, optical in ((0, 0), (0, 1), (1, 0), (1, 1)):
            vals = []
            for sign in (+1, -1):
                k = (x + sign*dk, y) if deriv == 0 else (x, y + sign*dk)
                h, v, _ = matrices(*k, ep)
                e, u = np.linalg.eigh(h)
                phases = np.angle(np.sum(u0.conj()*u, axis=0))
                u = u * np.exp(-1j*phases)[None, :]
                vnm = u[:, n].conj() @ v[optical] @ u[:, m]
                vals.append(-1j * vnm/(e[n]-e[m]))
            finite_difference = (vals[0]-vals[1])/(2*dk)
            exact = covariant_derivative(d, n, m, deriv, optical)
            self.assertLess(abs(finite_difference-exact), 2e-7)

    def test_s12_leading_part_vs_full_s10(self):
        ep = 1.03
        y0 = nodal_parameters(ep)['ky']
        previous = None
        for q in (1e-3, 1e-4, 1e-5):
            d = band_data(q, y0+.37*q, ep)
            exact = hermitian_connection(d, 0, 1, 0, 0, 1)
            singular = singular_connection(d.berry, 0, 1, 0, 0, 1)
            error = abs(exact-singular)/abs(singular)
            self.assertLess(error, 4*q)
            if previous is not None: self.assertLess(error, previous)
            previous = error

    def test_shg_cancel_double_and_single_frequency_residues(self):
        for ep in (.75, 1.03, 1.5):
            e0 = nodal_parameters(ep)
            omega = .7 * e0['epsilon0']
            residuals = []
            inters = []
            for q in (1e-3, 1e-4, 1e-5):
                d = band_data(q, e0['ky'] + .37*q, ep)
                p = shg_pieces(d.energies, d.berry, omega, (0,0,1))
                inters.append(abs(p.interband))
                residuals.append(abs(p.residual))
                self.assertLess(abs(p.double_residual), 8e-5)
                self.assertLess(abs(p.single_residual), 3e-5)
            # Both individual parts grow like 1/q, while the total remains finite.
            self.assertGreater(inters[-1]/inters[0], 90)
            self.assertLess(max(residuals)/min(residuals), 1.1)
            self.assertLess(residuals[-1]/inters[-1], 7e-5)

    def test_interband_and_mixed_scaling_separately(self):
        d = [benchmark(1.03, q) for q in (1e-3, 1e-4, 1e-5)]
        for name in ('inter_abs','mixed_abs'):
            self.assertAlmostEqual(d[0][name]*1e-3, d[1][name]*1e-4, delta=2e-7)
            self.assertAlmostEqual(d[1][name]*1e-4, d[2][name]*1e-5, delta=2e-8)
        self.assertLess(d[-1]['residual_over_inter'], d[0]['residual_over_inter']/90)

    def test_optical_index_exchange_band_swap_and_u1_gauge(self):
        ep = .9; n = nodal_parameters(ep); q=1e-4
        d = band_data(q, n['ky'] +.37*q, ep)
        omega = .7*n['epsilon0']
        p = shg_pieces(d.energies,d.berry,omega,(0,0,1))
        pg = shg_pieces(d.energies,d.berry,omega,(0,1,0))
        self.assertLess(abs(p.interband-pg.interband),1e-12)
        self.assertLess(abs(p.mixed_singular-pg.mixed_singular),1e-12)
        phases=np.array([.22,1.29,-.66]); factors=np.exp(-1j*phases[:,None]+1j*phases[None,:])
        r=d.berry*factors[...,None]
        pp=shg_pieces(d.energies,r,omega,(0,0,1))
        for field in p.__dataclass_fields__:
            self.assertLess(abs(getattr(p,field)-getattr(pp,field)),1e-12)
        ix=[0,2,1]
        swapped=shg_pieces(d.energies[ix],d.berry[ix,:,:][:,ix,:],omega,(0,0,1))
        for field in p.__dataclass_fields__:
            self.assertLess(abs(getattr(p,field)-getattr(swapped,field)),1e-10)

    def test_independent_nonzero_tensor_directions(self):
        n=nodal_parameters(1.03);q=2e-5
        d=band_data(q,n['ky']+.37*q,1.03)
        for axes in ((1,0,0),(0,1,1)):
            p=shg_pieces(d.energies,d.berry,.7*n['epsilon0'],axes)
            self.assertGreater(abs(p.interband),.03)
            self.assertLess(abs(p.residual)/abs(p.interband),1e-4)


if __name__ == '__main__': unittest.main()

