"""Regression of the migrated v30 four-band spin-MNHE algebra; no DC-transport claim."""
import unittest
import numpy as np
from qgmap.spin_multiband_v5 import four_band_audit,model_four_band,density_matrix_audit

class FourBandSpin(unittest.TestCase):
    def test_frozen_normalized_result(self):
        d=four_band_audit(N=23)
        self.assertLess(abs(d["total_rho"]-0.0007898128713846947),3e-13)
        self.assertLess(abs(d["direct_SM29"]-0.0007898128713851758),3e-13)
        self.assertLess(abs(d["closing_residual"]),1e-12)
    def test_phase_invariant(self):
        a=four_band_audit(N=19)
        b=four_band_audit(N=19,phase_seed=502)
        self.assertLess(abs(a["total_rho"]-b["total_rho"]),1e-11)
        self.assertLess(abs(b["closing_residual"]),1e-11)
    def test_four_terms_and_hermiticity(self):
        H,vx,vy,S=model_four_band(N=11)
        self.assertEqual(H.shape,(11,11,4,4))
        np.testing.assert_allclose(H,H.conj().swapaxes(-1,-2),atol=1e-13)
        self.assertEqual(len(density_matrix_audit(H,vx,vy,S)["j_rho_terms"]),4)
    def test_hall_antisymmetry(self):
        a=four_band_audit(N=19,a=0,b=1)
        b=four_band_audit(N=19,a=1,b=0)
        self.assertLess(abs(a["total_rho"]+b["total_rho"]),1e-10)
    def test_fails_closed(self):
        H,vx,vy,S=model_four_band(N=7)
        with self.assertRaises(ValueError):
            density_matrix_audit(H,vx,vy,S,T=0)
        with self.assertRaises(ValueError):
            density_matrix_audit(H,vx,vy,S,a=0,b=0)
        with self.assertRaises(ValueError):
            density_matrix_audit(np.zeros_like(H),vx,vy,S)
if __name__=="__main__": unittest.main()
