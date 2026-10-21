"""v43: noncrossing matched SCBA Hall I(a) only; NO full Hall/PHE acceptance."""
import unittest
import numpy as np
from qgmap.streda_clean_sea_v38 import qwz_model,I
from qgmap.scba_streda_v39 import SCBAConfig,solve_self_energy
from qgmap.current_vertex_v40 import pauli_vertex
from qgmap.exact_hall_v41 import torus
from qgmap.exact_streda_v42 import direct_surface_trace
from qgmap.scba_surface_v43 import scba_surface,exact_ensemble_surface,comparison

class MatchedSurfaceV43(unittest.TestCase):
    def test_frozen_v40_antisymmetric_RA_values(self):
        r=scba_surface(mesh=21, W=.12,eta=.12,mu=1.1)
        self.assertAlmostEqual(r['bare_surface_Ia'],-1.2788621328465,delta=3e-8)
        self.assertAlmostEqual(r['one_ladder_surface_Ia'],-1.4310682340406,delta=3e-8)
        self.assertLess(r['ra_density_resolvent_ward_residual'],1e-8)

    def test_one_sided_left_right_BSE_same_observable(self):
        for W in (0.,.03,.12):
            r=scba_surface(mesh=21,W=W)
            self.assertLess(abs(r['right_left_ladder_difference']),1e-9)
            self.assertLess(r['max_current_bse_residual'],1e-9)
            self.assertLess(r['ra_density_resolvent_ward_residual'],1e-8)

    def test_neumann_iteration_independent_ladder_solver(self):
        W=.02;mesh=17
        H,vx,vy=qwz_model(mesh,-1.,offset=(0.,0.))
        cfg=SCBAConfig(mesh=mesh,disorder_variance=W,eta=.12)
        gr=solve_self_energy(H,1.1,cfg)['G'];ga=gr.conj().swapaxes(-1,-2)
        solved=pauli_vertex(gr,ga,vy,W)
        C=np.zeros((2,2),complex)
        for _ in range(250):
            C=W*(gr@(vy+C)@ga).mean((0,1))
        self.assertLess(np.max(abs(C-solved['constant'])),1e-10)

    def test_clean_twist_mesh_reference(self):
        v=comparison(mesh=35,L=4,W=0.,mu=1.1,eta=.12,
                     seeds=(5,8,12),ntheta=7)
        self.assertEqual(v['exact_finite_ensemble']['number_independent_disorder_seeds'],1)
        self.assertEqual(v['exact_finite_ensemble']['seed_standard_error'],0.)
        self.assertAlmostEqual(v['scba']['one_ladder_surface_Ia'],v['scba']['bare_surface_Ia'],delta=1e-12)
        self.assertLess(v['absolute_error_single_ladder'],.004)

    def test_exact_disorder_surface_trace_same_v42_convention(self):
        W=.12;mu=1.1;eta=.12;theta=(0.,np.pi)
        H,x,y=torus(L=3,mass=-1.,W=W,seed=1947,theta=theta)
        direct=direct_surface_trace(H,x,y,mu,eta)['surface_Ia_hall_odd']
        values=[]
        for i in range(2):
            for j in range(2):
                h,a,b=torus(3,-1.,W,1947,(i*np.pi,j*np.pi))
                values.append(direct_surface_trace(h,a,b,mu,eta)['surface_Ia_hall_odd'])
        avg=exact_ensemble_surface(L=3,W=W,mu=mu,eta=eta,seeds=(1947,),ntheta=2)
        self.assertAlmostEqual(avg['exact_surface_Ia_mean'],np.mean(values),delta=1e-12)
        self.assertTrue(np.isfinite(direct))

    def test_weak_disorder_ladder_comparison(self):
        out=comparison(mesh=35,L=4,W=.03,mu=1.1,eta=.12,
                       seeds=tuple(range(1947,1971)),ntheta=7)
        self.assertTrue(out['ladder_improves_this_sample'])
        self.assertLess(out['absolute_error_single_ladder'],.09)
        self.assertGreater(out['exact_finite_ensemble']['seed_standard_error'],0.)

    def test_moderate_disorder_incomplete_agreement_is_recorded(self):
        out=comparison(mesh=35,L=4,W=.12,mu=1.1,eta=.12,
                       seeds=tuple(range(1947,1971)),ntheta=7)
        self.assertTrue(out['ladder_improves_this_sample'])
        self.assertGreater(out['absolute_error_single_ladder'],.03)
        self.assertEqual(out['scientific_status'],'PARTIAL')
        self.assertTrue(any('Streda II' in x for x in out['scba']['missing']))

    def test_scba_grid_offset_small(self):
        x=scba_surface(mesh=31,W=.12,offset=(0.,0.))
        y=scba_surface(mesh=31,W=.12,offset=(.17,.31))
        self.assertLess(abs(x['one_ladder_surface_Ia']-y['one_ladder_surface_Ia']),.01)

    def test_one_ladder_only_and_no_false_total(self):
        r=scba_surface(mesh=21,W=.12)
        self.assertLess(r['one_ladder_surface_Ia'],r['bare_surface_Ia'])
        self.assertNotIn('full_hall_conductivity',r)
        self.assertEqual(r['scientific_status'],'PARTIAL')

    def test_sample_error_uses_seeds_not_twists(self):
        a=exact_ensemble_surface(L=3,W=.12,seeds=(1947,1948,1949,1950),ntheta=3)
        self.assertEqual(a['number_independent_disorder_seeds'],4)
        self.assertEqual(a['twist_grid_per_seed'],9)
        self.assertAlmostEqual(a['seed_standard_error'],
            np.std(a['per_seed_twist_averages'],ddof=1)/2,delta=1e-12)

    def test_fail_closed(self):
        for args in ({'W':-.2},{'eta':0.},{'mesh':6},{'W':float('nan')}):
            with self.assertRaises(ValueError):scba_surface(**args)
        for args in ({'W':-.2},{'eta':0.},{'ntheta':1},{'seeds':(1,1)},{'seeds':()}):
            with self.assertRaises(ValueError):exact_ensemble_surface(**args)

if __name__=='__main__':unittest.main()
