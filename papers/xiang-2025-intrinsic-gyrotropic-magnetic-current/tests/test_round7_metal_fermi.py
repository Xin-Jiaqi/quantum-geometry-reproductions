"""Metal contour regression; not a claim of matching original published pixels."""
import unittest
import numpy as np
from igmc_lab.fermi_contour import (contour_roots_x,y_topology_splits,integrate_contour,
                                    integrate_conduction_fermi_line)
from igmc_lab.d_vector_independent import midpoint_dirac

class MetalFermiTests(unittest.TestCase):
    def test_two_root_nonmonotone_fermi_contour(self):
        roots=contour_roots_x(.5,1.)
        self.assertEqual(len(roots),2)
        for x in roots:
            v=.825;b=.25;M=.1;y=.5
            E=np.sqrt(v*v*(np.sin(x)**2+np.sin(y)**2)+(b+2*M*(2-np.cos(x)-np.cos(y)))**2)
            self.assertLess(abs(E-1.),1e-13)

    def test_topology_segments_partition(self):
        for mu in (.8,1.,1.05,1.2):
            edges=y_topology_splits(mu)
            self.assertEqual(edges[0],0.);self.assertEqual(edges[-1],float(np.pi))
            self.assertEqual(edges,sorted(edges))
            self.assertTrue(all(0<=x<=np.pi for x in edges))

    def test_independent_fermi_line_eq13_at_mu1(self):
        sea=integrate_contour(1.,outer_order=128,inner_order=64)
        fs=integrate_conduction_fermi_line(1.,order=256)
        self.assertLess(abs(sea['C_xy']-fs['C_xy_eq13']),2e-9)

    def test_chemical_potential_particle_hole_parity(self):
        for mu in (.8,1.,1.2):
            a=integrate_contour(mu,outer_order=64,inner_order=64)
            b=integrate_contour(-mu,outer_order=64,inner_order=64)
            self.assertAlmostEqual(a['C_xy'],b['C_xy'],places=14)
            self.assertAlmostEqual(a['D_peromega_xx'],b['D_peromega_xx'],places=14)

    def test_metal_finite_step_convergence(self):
        for mu in (.95,1.,1.05,1.1):
            a=integrate_contour(mu,outer_order=96,inner_order=64)
            b=integrate_contour(mu,outer_order=160,inner_order=64)
            self.assertLess(abs(a['C_xy']-b['C_xy']),3e-8)
            self.assertLess(abs(a['D_peromega_xx']-b['D_peromega_xx']),1e-8)

    def test_smooth_intervals_inner_gauss_convergence(self):
        a=integrate_contour(1.,outer_order=96,inner_order=32)
        b=integrate_contour(1.,outer_order=96,inner_order=64)
        self.assertLess(abs(a['C_xy']-b['C_xy']),1e-10)
        self.assertLess(abs(a['D_peromega_xx']-b['D_peromega_xx']),1e-10)

    def test_gapped_conduction_zero_displacement_nonzero(self):
        r=integrate_contour(0.,outer_order=128,inner_order=64)
        self.assertLess(abs(r['C_xy']),1e-11)
        self.assertGreater(r['D_peromega_xx'],.015)

    def test_nonresonant_frequency_and_full_Kubo(self):
        w=.0041356676943898555
        r=integrate_contour(1.,outer_order=96,inner_order=64,photon_eV=w)
        self.assertLess(abs(r['K_real_xy']/r['C_xy']-1),1e-4)
        self.assertLess(abs(r['K_imag_xx']/(r['D_peromega_xx']*w)-1),1e-4)

    def test_domain_odd_displacement_even_conduction(self):
        pos=integrate_contour(1.05,outer_order=64,mass_sign=1)
        neg=integrate_contour(1.05,outer_order=64,mass_sign=-1)
        self.assertAlmostEqual(pos['C_xy'],neg['C_xy'],places=12)
        self.assertAlmostEqual(pos['D_peromega_xx'],-neg['D_peromega_xx'],places=12)

    def test_uniform_grid_metal_not_used_as_convergence_gold(self):
        ref=integrate_contour(1.,outer_order=160)['D_peromega_xx']
        naive=midpoint_dirac(201,1.,spin='x')['sigma_D_per_omega']
        self.assertGreater(abs(ref-naive),1e-7)

if __name__=='__main__':unittest.main()
