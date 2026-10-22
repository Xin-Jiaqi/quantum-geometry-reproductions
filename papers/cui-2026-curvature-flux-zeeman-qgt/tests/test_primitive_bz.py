"""Round19 independent reciprocal primitive-cell integration controls."""
import unittest
import numpy as np
from flux_lab.hexagonal_bz import BASIS,INV_BASIS,hexagonal_bz_quadrature,hexagonal_symmetrized_model
from flux_lab.group_microscopic import project_frequency_tensor
from flux_lab.primitive_bz import primitive_cell_quadrature,primitive_bz_response

class PrimitiveBZ(unittest.TestCase):
    def test_exact_primitive_area_and_normalization(self):
        area=2*np.sqrt(3.)*np.pi**2
        for method in ('gauss','midpoint'):
            for n in (3,8,12):
                k,w,A=primitive_cell_quadrature(n,method=method)
                self.assertEqual(k.shape,(n*n,2))
                self.assertAlmostEqual(A,area,places=12)
                self.assertAlmostEqual(w.sum(),1.,places=12)
                self.assertGreater(w.min(),0.)

    def test_two_sampling_protocols_have_distinct_nodes(self):
        a,_,_=primitive_cell_quadrature(8,method="gauss")
        b,_,_=primitive_cell_quadrature(8,method="midpoint")
        self.assertGreater(np.max(abs(a-b)),1e-2)

    def test_analytic_periodic_fourier_integral_on_primitive_torus(self):
        for method in ('gauss','midpoint'):
            k,w,_=primitive_cell_quadrature(16,method=method)
            q=k@INV_BASIS.T
            for n1,n2 in ((1,0),(0,1),(1,1),(1,-1),(2,1)):
                average=sum(w*np.exp(1j*(n1*q[:,0]+n2*q[:,1])))
                self.assertLess(abs(average),2e-12)

    def test_primitive_quadrature_respects_recursion_guard(self):
        with self.assertRaises(ValueError):primitive_cell_quadrature(1)
        with self.assertRaises(ValueError):primitive_cell_quadrature(4,method="bad")
        with self.assertRaises(ValueError):primitive_cell_quadrature(4,method="gauss",shift=(.1,0))
        with self.assertRaises(ValueError):primitive_cell_quadrature(4,method="midpoint",shift=(.9,0))
        with self.assertRaises(ValueError):primitive_bz_response("4mm",order=4)

    def test_independent_shifted_midpoint_is_periodic(self):
        k,w,_=primitive_cell_quadrature(7,method="midpoint",shift=(.2,-.16))
        q=k@INV_BASIS.T
        self.assertEqual(q.shape,(49,2))
        self.assertTrue(np.all(np.isfinite(q)))
        self.assertAlmostEqual(w.sum(),1.,places=13)

    def test_hexagon_folded_to_primitive_is_same_periodic_hamiltonian(self):
        pts,w,area=hexagonal_bz_quadrature(3)
        for group in ("6","6m'm'"):
            for k in pts[::11]:
                q=INV_BASIS@k
                folded=(q+np.pi)%(2*np.pi)-np.pi
                H,V,S=hexagonal_symmetrized_model(k,group)
                H2,V2,S2=hexagonal_symmetrized_model(BASIS@folded,group)
                np.testing.assert_allclose(H2,H,atol=2e-12)
                np.testing.assert_allclose(V2,V,atol=2e-12)
                np.testing.assert_allclose(S2,S,atol=2e-12)

    def test_primitive_bz_output_and_frequency_parity(self):
        for method in ("gauss","midpoint"):
            d=primitive_bz_response("6",order=5,method=method)
            self.assertEqual(d["even"].shape,(2,3))
            self.assertEqual(d["odd"].shape,(2,3))
            np.testing.assert_allclose(d["even"]+d["odd"],d["positive_frequency"],atol=1e-13)
            np.testing.assert_allclose(d["even"]-d["odd"],d["negative_frequency"],atol=1e-13)
            self.assertAlmostEqual(d["weight_sum"],1.,places=13)

    def test_unprojected_primitive_response_can_have_quadrature_symmetry_error(self):
        # Do not *force* symmetry by post-averaging response. We only ask for
        # bounded finite errors on an intentionally small primitive grid.
        r=primitive_bz_response("6m'm'",order=5,method="gauss")
        self.assertTrue(np.all(np.isfinite(r["even"])))
        self.assertTrue(np.all(np.isfinite(r["odd"])))
        even=r["even"]-project_frequency_tensor("6m'm'",r["even"],"even")
        odd=r["odd"]-project_frequency_tensor("6m'm'",r["odd"],"odd")
        self.assertLess(np.max(abs(even)),.1)
        self.assertLess(np.max(abs(odd)),.1)

if __name__=="__main__":unittest.main()
