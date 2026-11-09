"""Physical consistency of the Hamiltonian controls, independent of a fitted sign."""
import unittest
from types import SimpleNamespace

import numpy as np

from scripts.compare_reversal_mechanism import (
    AuditBilayer, AuditSpec, ab_ba_shift, adjoint, km_response, transplant,
)
from scripts.run_fixed_stack_layer_cxx import response


class TestSpectralControl(unittest.TestCase):
    def setUp(self):
        self.k = np.array([[.13, .29], [8.1, -3.2], [-5.8, 12.2]])
        self.axis = np.array([.5, np.sqrt(3)/2])
        self.models = [AuditBilayer(AuditSpec(shift=ab_ba_shift(.25),
                       staggered_eV=.045, tunneling_eV=t)) for t in (.010, .020)]

    def test_velocity_against_independent_finite_difference(self):
        for donor, target_model in ((0, 1), (1, 0)):
            a, b = self.models[donor], self.models[target_model]
            e, u = np.linalg.eigh(a.matrices(self.k))
            eq, uq = np.linalg.eigh(b.matrices(self.k))
            v = sum(self.axis[j]*a.matrices(self.k, axis=j) for j in range(2))
            vq = sum(self.axis[j]*b.matrices(self.k, axis=j) for j in range(2))
            slopes = np.diagonal(adjoint(uq)@vq@uq, axis1=1, axis2=2).real
            _, analytic = transplant(e, u, eq, v, slopes)
            hs = []
            for sign in (1, -1):
                k = self.k + sign*1e-5*self.axis
                ep, up = np.linalg.eigh(a.matrices(k))
                hs.append(transplant(ep, up, np.linalg.eigvalsh(b.matrices(k))))
            np.testing.assert_allclose(analytic, (hs[0]-hs[1])/2e-5, atol=1e-8)

    def test_degenerate_rotation_and_rejection(self):
        e = np.array([[-3, -3, -1, -1, 1, 1, 3, 3]], dtype=float)
        u = np.eye(8, dtype=complex)[None]
        rotation = u.copy()
        block = np.array([[1, 1j], [1j, 1]])/np.sqrt(2)
        for start in range(0, 8, 2):
            rotation[0, start:start+2, start:start+2] = block
        np.testing.assert_allclose(transplant(e, u, e*2),
                                   transplant(e, rotation, e*2), atol=1e-14)
        split = e.copy()
        split[0, 0] += .1
        with self.assertRaisesRegex(ValueError, 'degenerate donor'):
            transplant(e, u, split)

    def test_response_matches_existing_kernel_and_basis_rotation(self):
        model = self.models[0]
        unit = AuditBilayer(AuditSpec(shift=ab_ba_shift(.25), staggered_eV=.045,
                                    tunneling_eV=1.))
        spins = [sum(self.axis[j]*model.spins[name+'_'+axis]
                 for j, axis in enumerate(('x', 'y')))
                 for name in ('bottom', 'top', 'total')]
        sz = model.spins['total_z']
        weights = np.array([.2, .3, .5])
        h = model.matrices(self.k)
        v = sum(self.axis[j]*model.matrices(self.k, axis=j) for j in range(2))
        step = .03125
        actual = km_response(h, v, spins, sz, weights, step, model.mono.g_factor)[0]
        spin_map = dict(zip(('bottom_x', 'top_x', 'total_x'), spins), total_z=sz)
        v_unit = sum(self.axis[j]*unit.matrices(self.k, axis=j) for j in range(2))
        expected = response(SimpleNamespace(mono=model.mono, spins=spin_map),
                            unit.matrices(self.k), v_unit, weights, 10., 0., step)
        np.testing.assert_allclose(actual, [expected[n] for n in ('bottom', 'top', 'total')], atol=1e-12)
        rng = np.random.default_rng(17)
        q, _ = np.linalg.qr(rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8)))
        rotated = km_response(adjoint(q)@h@q, adjoint(q)@v@q,
                  [adjoint(q)@s@q for s in spins], adjoint(q)@sz@q,
                  weights, step, model.mono.g_factor)[0]
        np.testing.assert_allclose(actual, rotated, atol=1e-12)
        scaled = km_response(2*h, 2*v, spins, sz, weights, 2*step, model.mono.g_factor)[0]
        np.testing.assert_allclose(scaled, actual/4, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
