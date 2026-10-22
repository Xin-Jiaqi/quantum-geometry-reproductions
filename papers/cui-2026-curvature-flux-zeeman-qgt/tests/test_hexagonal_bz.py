"""Round17 triangular reciprocal lattice and first-hexagon BZ quadrature."""
import unittest
import numpy as np
from flux_lab.magnetic_symmetry import matrices
from flux_lab.group_microscopic import build_group,project_frequency_tensor
from flux_lab.hexagonal_bz import (
    BASIS,fractional_action,hexagonal_compatible,
    hexagonal_symmetrized_model,hexagonal_bz_quadrature,
    hexagonal_local_residual,hexagonal_bz_response)

class HexagonalBrillouinZone(unittest.TestCase):
    def test_c6_and_mirror_are_integer_lattice_automorphisms(self):
        for name in ("6","6mm","6m'm'","3m'","m'm'2"):
            self.assertTrue(hexagonal_compatible(name),name)
        self.assertFalse(hexagonal_compatible("4mm"))

    def test_fractional_rotation_closes_modulo_reciprocal_periods(self):
        g=next(g for g in build_group("6") if not g.anti and np.linalg.norm(g.R-np.eye(2))>1)
        A=fractional_action(g)
        np.testing.assert_allclose(np.linalg.matrix_power(A,6),np.eye(2),atol=0)

    def test_hexagonal_quadrature_area_and_normalized_weights(self):
        for order in (2,3,5):
            k,w,area=hexagonal_bz_quadrature(order)
            self.assertEqual(k.shape,(12*order*order,2))
            self.assertAlmostEqual(float(w.sum()),1.,places=13)
            self.assertAlmostEqual(area,2*np.sqrt(3)*np.pi**2,places=12)
            self.assertGreater(float(np.min(w)),0.)

    def test_quadrature_is_exactly_d6_symmetric_as_a_weighted_set(self):
        k,w,_=hexagonal_bz_quadrature(4)
        for R in (matrices()["C6"],np.diag([-1.,1.])):
            kp=k@R.T
            for j,p in enumerate(kp):
                distances=np.linalg.norm(k-p,axis=1);i=int(np.argmin(distances))
                self.assertLess(distances[i],2e-13)
                self.assertAlmostEqual(float(w[j]),float(w[i]),places=13)

    def test_hamiltonian_and_cartesian_velocity_are_reciprocal_periodic(self):
        q=np.array([.37,-.81]);k=BASIS@q
        H,V,S=hexagonal_symmetrized_model(k,"6mm")
        for shift in (np.array([2*np.pi,0]),np.array([0,2*np.pi]),np.array([2*np.pi,-2*np.pi])):
            Hp,Vp,_=hexagonal_symmetrized_model(k+BASIS@shift,"6mm")
            np.testing.assert_allclose(Hp,H,atol=2e-12)
            np.testing.assert_allclose(Vp,V,atol=2e-12)

    def test_local_hexagonal_kubo_covariance_for_rotations_mirrors_and_antiunitaries(self):
        for name in ("6","6mm","6m'm'","3m'"):
            with self.subTest(group=name):
                r=hexagonal_local_residual(name)
                self.assertLess(max(r[x] for x in ("H","velocity","chi","alpha")),2e-10,r)

    def test_first_bz_points_are_inside_all_six_ws_halfplanes(self):
        k,w,_=hexagonal_bz_quadrature(5)
        reciprocal=2*np.pi*BASIS
        nearest=[reciprocal[:,0],reciprocal[:,1],reciprocal[:,1]-reciprocal[:,0]]
        for G in nearest:
            self.assertLess(float(np.max(k@G)),np.dot(G,G)/2+1e-12)
            self.assertGreater(float(np.min(k@G)),-np.dot(G,G)/2-1e-12)

    def test_bz_integral_respects_even_odd_magnetic_projectors(self):
        for name in ("6","6m'm'"):
            with self.subTest(group=name):
                r=hexagonal_bz_response(name,order=3)
                np.testing.assert_allclose(r["even"],project_frequency_tensor(name,r["even"],"even"),atol=3e-10,rtol=3e-10)
                np.testing.assert_allclose(r["odd"],project_frequency_tensor(name,r["odd"],"odd"),atol=3e-10,rtol=3e-10)
                self.assertAlmostEqual(float(r["weights"].sum()),1.,places=12)

    def test_square_lattice_group_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"triangular reciprocal"):
            hexagonal_bz_response("4mm",order=3)

    def test_invalid_quadrature_order_is_rejected(self):
        with self.assertRaises(ValueError):hexagonal_bz_quadrature(1)

if __name__=="__main__":unittest.main()
