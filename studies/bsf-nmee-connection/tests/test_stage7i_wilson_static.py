"""Regression: Stage7i Wilson occupied-bundle curvature vs Kubo static response."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingBilayer,SlidingSpec
from bsf_nmee.mixed_wilson_geometry import (
    _field_loop_phase, wilson_static_kernel,
    compare_valley_wilson_kubo)


def spec_u(u,t):
    return replace(SlidingSpec(),shift=(.23*u,.41*u),tunneling_eV=t)


class WilsonMixedCurvatureTests(unittest.TestCase):
    def test_local_plaquette_steps_converge_quadratically(self):
        model=SlidingBilayer(spec_u(.4,.005))
        k=np.array([[1/3,2/3]])@model.model().reciprocal+np.array([[.05,.02]])
        from bsf_nmee.valley_kp import static_kernel
        H=model.matrices(k)
        V=np.stack([model.matrices(k,axis=j) for j in (0,1)])
        exact=static_kernel(H,V,[model.spins['total_x'],model.spins['total_y']],
                            model.spins['total_z'],model.spec.mu_eV,
                            g_factor=model.mono.g_factor)[0,0,0]
        errs=[]
        for dk,dl in ((.02,.0005),(.01,.00025),(.005,.000125)):
            ker,gap=wilson_static_kernel(model,k,dk=dk,dl=dl)
            errs.append(abs(ker[0]-exact))
            self.assertGreater(gap,.003)
        self.assertGreater(errs[0]/errs[1],3.5)
        self.assertGreater(errs[1]/errs[2],3.5)
        self.assertLess(errs[2]/abs(exact),.001)

    def test_adaptive_two_valley_integral_agreement(self):
        for u,t,sgn in ((.1,.005,-1),(.4,.005,1),(.4,.02,1)):
            r=compare_valley_wilson_kubo(spec_u(u,t),
                      radius=.5,radial_n=16,angular_n=40,
                      dk=.0025,dl=.0000625)
            self.assertLess(r['relative_difference'],.0006)
            self.assertGreater(sgn*r['wilson'],0)
            self.assertGreater(sgn*r['kubo'],0)
            self.assertGreater(r['minimum_sampled_mu_gap_eV'],.002)

    def test_mirror_sliding_partner_is_opposite(self):
        a=compare_valley_wilson_kubo(spec_u(.4,.005),
                radius=.5,radial_n=12,angular_n=32)
        b=compare_valley_wilson_kubo(spec_u(-.4,.005),
                radius=.5,radial_n=12,angular_n=32)
        self.assertLess(abs(a['wilson']+b['wilson']),1e-9)

    def test_y_velocity_agrees(self):
        r=compare_valley_wilson_kubo(spec_u(.4,.005),spin_axis='y',
                k_axis=1,radius=.5,radial_n=12,angular_n=32)
        self.assertLess(r['relative_difference'],.002)


if __name__=='__main__':
    unittest.main()
