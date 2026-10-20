"""Independently check the origin and normalization of Fig 2(d)'s displacement term."""
import unittest
import numpy as np
from igmc_lab.independent_checks import finite_field_wilson_derivative,orbital_liouville_kernel
from igmc_lab.core import dirac_response,energy_from_frequency

class TestIndependentDisplacement(unittest.TestCase):
    def test_wilson_loop_finite_field_derivative_matches_main_eq12(self):
        # Wilson loop depends solely on finite-field eigenspinor overlaps.
        # Link error decreases as the kx mesh is refined.
        expected=dirac_response(151,mu=0.,alpha='x',beta='x')['sigma_D_per_omega']
        x=finite_field_wilson_derivative(151,b_eV=1e-4)
        self.assertAlmostEqual(x/expected,1.,delta=3e-4)
        coarse=finite_field_wilson_derivative(51,b_eV=1e-4)
        self.assertLess(abs(x-expected),abs(coarse-expected))

    def test_wilson_derivative_field_amplitude_stability(self):
        x=finite_field_wilson_derivative(101,b_eV=1e-4)
        y=finite_field_wilson_derivative(101,b_eV=1e-3)
        self.assertLess(abs(x-y),2e-8)

    def test_direct_orbital_liouville_complex_kubo(self):
        w=energy_from_frequency(1e12,protocol='cyclic_hz')
        benchmark=dirac_response(51,0,omega_eV=w,alpha='x',beta='x')
        K=orbital_liouville_kernel(51,omega_eV=w,alpha='x',spin='x')
        self.assertAlmostEqual(K.real,benchmark['Kubo_real'],delta=2e-11)
        self.assertAlmostEqual(K.imag,benchmark['Kubo_imag'],delta=2e-11)
        self.assertGreater(K.imag,1e-5)

    def test_independent_liouville_frequency_even_odd(self):
        w=energy_from_frequency(1e12,protocol='cyclic_hz')
        plus=orbital_liouville_kernel(31,omega_eV=w)
        minus=orbital_liouville_kernel(31,omega_eV=-w)
        self.assertAlmostEqual(plus.real,minus.real,delta=1e-12)
        self.assertAlmostEqual(plus.imag,-minus.imag,delta=1e-12)

    def test_sm_eq34_local_intermediate_algebra_double_count(self):
        # SM p.6 after Eq (34) explicitly says Q=2 Re(z), but SM's
        # intermediate cos coefficient already has a full +-omega sum.
        # For real z and finite omega the result is exactly twofold.
        omega=.004135667694;gap=.6;real_z=.37
        from_first_line=real_z*(1/(omega-gap)-1/(omega+gap))
        middle_as_printed=(2*real_z)*(1/(omega-gap)-1/(omega+gap))
        self.assertAlmostEqual(middle_as_printed/from_first_line,2.,places=12)
        # The correct main Eq.(5) Q=Re(z) agrees with the first line.
        main_q=real_z
        self.assertAlmostEqual(main_q*(1/(omega-gap)-1/(omega+gap)),from_first_line)

if __name__ == '__main__': unittest.main()
