"""Round6 eigenvector-free identities, nonresonant complex response and TR."""
import unittest
import numpy as np
from igmc_lab.d_vector_independent import (dvector_dirac,pauli_trace_geometry,
                                           midpoint_dirac,aggregate_response)
from igmc_lab.core import (dirac_hamiltonian,pair_kernels,dirac_response,
                           energy_from_frequency)

OMEGA=energy_from_frequency(1e12,protocol='cyclic_hz')


class DVectorIndependentTests(unittest.TestCase):
    def test_eigenvector_free_pointwise_all_components_both_domains(self):
        rng=np.random.default_rng(325)
        x=rng.uniform(-3.0,3.0,(9,14));y=rng.uniform(-3.0,3.0,(9,14))
        for sign in (1,-1):
            H,dh=dirac_hamiltonian(x,y,mass_sign=sign)
            d,dd=dvector_dirac(x,y,mass_sign=sign)
            for alpha in ('x','y'):
                for spin in ('x','y','z'):
                    with self.subTest(sign=sign,alpha=alpha,spin=spin):
                        old=pair_kernels(H,dh,alpha,spin)
                        new=pauli_trace_geometry(d,dd,alpha,spin)
                        self.assertLess(np.max(abs(old['Z01']-new['Z01'])),2e-13)
                        self.assertLess(np.max(abs(2*old['Q01']/old['gap']-new['D_per_omega01'])),2e-13)

    def test_complex_kubo_full_frequency_both_domains(self):
        for sign in (1,-1):
            for spin in ('x','y','z'):
                with self.subTest(sign=sign,spin=spin):
                    x=midpoint_dirac(51,0.,OMEGA,alpha='x',spin=spin,mass_sign=sign)
                    y=dirac_response(51,0.,omega_eV=OMEGA,alpha='x',beta=spin,mass_sign=sign)
                    for key in ('sigma_C','sigma_D_per_omega','Kubo_real','Kubo_imag'):
                        self.assertLess(abs(x[key]-y[key]),2e-13)

    def test_tr_pairs_conjugate_full_response(self):
        for spin in ('x','y','z'):
            for axis in ('x','y'):
                for w in (OMEGA,3*OMEGA):
                    x=midpoint_dirac(41,0.,w,alpha=axis,spin=spin,mass_sign=1)
                    y=midpoint_dirac(41,0.,w,alpha=axis,spin=spin,mass_sign=-1)
                    with self.subTest(spin=spin,axis=axis,w=w):
                        self.assertLess(abs(x['Kubo_real']-y['Kubo_real']),1e-12)
                        self.assertLess(abs(x['Kubo_imag']+y['Kubo_imag']),1e-12)

    def test_insulator_dperomega_gauss_reference(self):
        # This constant is a regression target from independent Gauss-Legendre
        # quadrature, not fitted to Fig.2(d). Uniform midpoint at N101 matches.
        value=midpoint_dirac(101,0.,OMEGA,alpha='x',spin='x')
        self.assertAlmostEqual(value['sigma_D_per_omega'],0.015234078015630925,places=7)

    def test_midpoint_displacement_odd_under_mass_reversal(self):
        a=midpoint_dirac(101,0.,OMEGA,alpha='x',spin='x',mass_sign=1)
        b=midpoint_dirac(101,0.,OMEGA,alpha='x',spin='x',mass_sign=-1)
        self.assertAlmostEqual(a['sigma_D'],-b['sigma_D'],places=13)

    def test_dvector_bad_gap_and_indices_fail_closed(self):
        with self.assertRaises(ValueError):pauli_trace_geometry(np.zeros((2,3)),{'x':np.ones((2,3))})
        d=np.ones((2,3)); grad={'x':np.ones((2,3))}
        with self.assertRaises(ValueError):pauli_trace_geometry(d,grad,alpha='z')
        with self.assertRaises(ValueError):pauli_trace_geometry(d,grad,beta='invalid')
        with self.assertRaises(ValueError):dvector_dirac(0.,0.,a=0.)
        with self.assertRaises(ValueError):aggregate_response({'gap':np.array([1.]),'Z01':np.array([1.]),'D_per_omega01':np.array([1.])},0.,1.1)

    def test_finite_frequency_quadrature_not_manually_doubled(self):
        # Omega->0 consistency and exact 1/(1-(omega/gap)^2) correction.
        x=midpoint_dirac(91,0.,OMEGA,alpha='x',spin='x')
        self.assertGreater(x['Kubo_imag'],0.)
        self.assertGreater(x['Kubo_imag']/x['sigma_D'],1.)
        self.assertLess(x['Kubo_imag']/x['sigma_D'],1.001)

if __name__=='__main__':unittest.main()
