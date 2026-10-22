"""Eq.(27) real-spin vorticity: independent bands, k gradients, Hessian."""
import unittest
import numpy as np
from flux_lab.core import PAULI
from flux_lab.spin_current import spin_orbital_benchmark, spin_current_derivative
from flux_lab.spin_vorticity import (spin_vorticity_spectral, vorticity_with_hessian,
                                      vorticity_finite_difference)

class SpinVorticity(unittest.TestCase):
    def test_periodic_curl_spectral_vs_independent_finite_difference(self):
        for x,y in [(.31,-.22),(-.48,.27),(.63,.11)]:
            h,v,K,s=spin_orbital_benchmark(x,y)
            w=spin_vorticity_spectral(h,v,s)['bar_omegaA_z']
            np.testing.assert_allclose(w,vorticity_finite_difference(x,y),atol=3e-9)
    def test_linear_same_identity(self):
        for x,y in [(.31,-.22),(.16,.49)]:
            h,v,K,s=spin_orbital_benchmark(x,y,periodic=False)
            w=spin_vorticity_spectral(h,v,s)['bar_omegaA_z']
            np.testing.assert_allclose(w,vorticity_finite_difference(x,y,periodic=False),atol=3e-9)
    def test_full_derivative_curl_is_spectral_vorticity(self):
        for x,y in [(.31,-.22),(-.45,.13)]:
            h,v,K,s=spin_orbital_benchmark(x,y)
            d=vorticity_with_hessian(h,v,K,s)
            np.testing.assert_allclose(d['curl_half_from_full_derivative'],d['bar_omegaA_z'],atol=3e-13)
    def test_explicit_hessian_curl_cancels_but_symmetric_correction_does_not(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        d=vorticity_with_hessian(h,v,K,s)
        self.assertLess(np.max(abs(d['vertex_curl_half'])),1e-14)
        self.assertGreater(np.max(abs(d['symmetric_hessian_correction'])),.1)
    def test_four_sector_reconstruction(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        d=spin_vorticity_spectral(h,v,s)
        np.testing.assert_allclose(d['intrinsic'],d['intrinsic_from_sectors'],atol=1e-13)
        np.testing.assert_allclose(d['reduced_RTA'],d['reduced_RTA_from_sectors'],atol=1e-13)
    def test_real_spin_rotational_basis_gauge(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        rng=np.random.default_rng(43)
        z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        U,_=np.linalg.qr(z)
        tr=lambda M:U.conj().T@M@U
        a=spin_vorticity_spectral(h,v,s)
        b=spin_vorticity_spectral(tr(h),np.stack([tr(w) for w in v]),tr(s))
        for key in ('bar_omegaA','bar_gN','intrinsic','reduced_RTA'):
            np.testing.assert_allclose(a[key],b[key],atol=3e-13)
    def test_energy_shift_does_not_change_vorticity(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        a=spin_vorticity_spectral(h,v,s)
        b=spin_vorticity_spectral(h+3.8*np.eye(4),v,s)
        np.testing.assert_allclose(a['bar_omegaA_z'],b['bar_omegaA_z'],atol=1e-13)
    def test_hessian_required_only_for_full_derivative(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        spin_vorticity_spectral(h,v,s)
        with self.assertRaises(ValueError):vorticity_with_hessian(h,v,None,s)
    def test_round14_frozen_vorticity_and_pointwise_errors(self):
        import json
        from pathlib import Path
        root=Path(__file__).resolve().parents[1]
        ref=json.loads((root/'results/frozen/round14_vorticity_mpg_reference.json').read_text())
        row=ref['local_spin_vorticity']['+0.31_-0.22']
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        a=vorticity_with_hessian(h,v,K,s)
        np.testing.assert_allclose(a['bar_omegaA_z'],row['band_omega_bar_A_z'],atol=2e-13)
        self.assertLess(np.max(abs(a['bar_omegaA_z']-vorticity_finite_difference(.31,-.22))),1e-10)
        self.assertEqual(ref['magnetic_group_count'],31)
        self.assertEqual(ref['projected_candidate_spaces'],124)
        self.assertEqual(ref['source_table_mismatch_count'],6)

    def test_periodic_band_vorticity_integrates_to_zero_on_torus(self):
        # Smooth, nondegenerate, periodic spin velocity: Stokes gives
        # an unoccupied full-BZ band integral of curl equal to zero.
        for nmesh in (12,24):
            total=np.zeros(4)
            for kx in -np.pi+(np.arange(nmesh)+.5)*2*np.pi/nmesh:
                for ky in -np.pi+(np.arange(nmesh)+.5)*2*np.pi/nmesh:
                    h,v,_,s=spin_orbital_benchmark(kx,ky)
                    total+=spin_vorticity_spectral(h,v,s)['bar_omegaA_z']
            self.assertLess(np.max(abs(total/nmesh**2)),1e-11)
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        self.assertGreater(np.max(abs(spin_vorticity_spectral(h,v,s)['bar_omegaA_z'])),1e-3)

    def test_zero_spin_or_synthetic_two_band_not_conflate_spin(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        with self.assertRaises(ValueError):spin_vorticity_spectral(h,v,np.eye(2))
    def test_analytic_derivative_vs_hessian_periodic(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        d=spin_current_derivative(h,v,K,s)
        q=spin_vorticity_spectral(h,v,s)
        np.testing.assert_allclose(q['reduced_RTA'],d['geometry_only'],atol=3e-13)

if __name__=='__main__':unittest.main()
