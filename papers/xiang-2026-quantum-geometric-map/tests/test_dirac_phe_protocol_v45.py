"""v45 tests target the ORIGINAL symmetric EB PHE, not QWZ antisymmetric Hall."""
import unittest
import numpy as np
from qgmap.dirac_phe_protocol_v45 import (
    hamiltonian,shift_identity,zeeman_curvature,printed_eq11_geometry,
    angular_eq8_geometry,equilibrium_drude_xy,drude_delta_derivative_terms,audit)


class OriginalDiracPHEV45(unittest.TestCase):
    def test_two_helicities_pure_momentum_translation(self):
        for h in ('main','supplementary'):
            for x,y,d in ((.31,-.7,.12),(.41,.3,-.21),(-.1,.05,0.)):
                self.assertLess(shift_identity(x,y,.125,d,h),1e-13)

    def test_main_zeeman_curvature_eigenvectors_match_printed_source(self):
        for x,y in ((.3,.4),(.11,-.28),(-.27,.41)):
            r=np.hypot(x,y)
            for band,sgn in ((1,1),(0,-1)):
                zxx,zyx=zeeman_curvature(x,y,band,'main')
                self.assertAlmostEqual(zxx,sgn*x*y/r**3,delta=1e-12)
                self.assertAlmostEqual(zyx,-sgn*x*x/r**3,delta=1e-12)

    def test_eq8_angular_quadrature_vs_printed_eq11_sign_branch(self):
        for beta in (.08,.125,.3,.6):
            for sign in (-1,1):
                eq8,err=angular_eq8_geometry(beta,sign)
                eq11=sign*printed_eq11_geometry(beta)
                self.assertLess(abs(eq8+eq11),2e-12)
                self.assertGreater(abs(eq8),1e-5)
                self.assertLess(err,1e-10)

    def test_chemical_potential_odd_parity_and_beta_zero_limit(self):
        for beta in (.025,.125,.45):
            plus,_=angular_eq8_geometry(beta,1)
            minus,_=angular_eq8_geometry(beta,-1)
            self.assertAlmostEqual(plus,-minus,delta=1e-12)
        self.assertLess(abs(printed_eq11_geometry(.000001)),3e-8)

    def test_full_plane_finite_B_symmetry(self):
        for mu in (.3,-.3):
            for helicity in ('main','supplementary'):
                vals=[equilibrium_drude_xy(mu=mu,delta=d,helicity=helicity,grid=301)
                      for d in (-.1,0.,.1)]
                self.assertLess(max(abs(z) for z in vals),1e-7)
                self.assertLess(abs(vals[0]-vals[2]),2e-7)

    def test_cutoff_artifact_decreases_when_extended(self):
        short=equilibrium_drude_xy(mu=.3,thermal=.08,delta=.24,cutoff=.45,grid=301)
        long=equilibrium_drude_xy(mu=.3,thermal=.08,delta=.24,cutoff=2.5,grid=301)
        self.assertGreater(abs(short),1e-4)
        self.assertLess(abs(long),abs(short)*1e-4)

    def test_drift_velocity_and_occupation_contributions_cancel(self):
        z=drude_delta_derivative_terms(grid=600)
        self.assertGreater(z['velocity_change'],.004)
        self.assertLess(z['equilibrium_occupation_change'],-.004)
        self.assertLess(abs(z['full_derivative']),1e-6)

    def test_finite_difference_and_derivative_consistent(self):
        z=audit(grid=301)
        self.assertEqual(z['scientific_status'],'PARTIAL')
        self.assertLess(abs(z['finite_B_central_dsigma_d_delta_per_tau']),1e-7)
        self.assertLess(abs(z['eq8_plus_eq11_sign_residual']),1e-11)
        self.assertIn('protocol',z)
        self.assertIn('scientific_open',z)

    def test_magnetic_and_current_conventions_explicit(self):
        h=hamiltonian(.21,.24,delta=.05)
        z=hamiltonian(.21,.29,delta=0.)
        np.testing.assert_allclose(h,z,atol=1e-13)
        self.assertNotAlmostEqual(printed_eq11_geometry(.125),0.,delta=1e-5)

    def test_reject_invalid_or_undefined_cases(self):
        with self.assertRaises(ValueError):hamiltonian(0,0,helicity='other')
        with self.assertRaises(ValueError):zeeman_curvature(0,0)
        with self.assertRaises(ValueError):angular_eq8_geometry(.125,0)
        with self.assertRaises(ValueError):printed_eq11_geometry(1.)
        with self.assertRaises(ValueError):drude_delta_derivative_terms(grid=301)
        with self.assertRaises(ValueError):equilibrium_drude_xy(thermal=0)

if __name__=='__main__':unittest.main()
