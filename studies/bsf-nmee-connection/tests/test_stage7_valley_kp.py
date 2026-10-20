"""Stage-7 scientific regression: 8-band local valley Taylor models."""
import json
from pathlib import Path
import unittest
import numpy as np
from dataclasses import replace
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.valley_kp import ValleyKP,compare_patch,static_kernel
from bsf_nmee.valley_quadrature import valley_quadrature
from bsf_nmee.symmetry_scaling import static_inplane_matrix
from scripts.run_stage3_sign_reversal import spec_at

ROOT=Path(__file__).resolve().parents[1]

class Stage7Physics(unittest.TestCase):
    def setUp(self):
        self.obj=SlidingBilayer(replace(spec_at(.1),tunneling_eV=.005))
        self.center=np.array([1/3,2/3])@self.obj.model().reciprocal
        self.kp=ValleyKP.from_microscopic(self.obj,self.center)

    def test_quadratic_H_is_more_accurate_than_linear_near_valley(self):
        k=self.center+np.array([.1,.06]);exact=self.obj.matrices(k[None,:])[0]
        e1=np.max(abs(exact-self.kp.matrices(k,order=1)[0]))
        e2=np.max(abs(exact-self.kp.matrices(k,order=2)[0]))
        self.assertLess(e2,e1*.02)
        self.assertLess(e2,1e-6)
        self.assertLess(np.max(abs(self.kp.h2-self.kp.h2.swapaxes(0,1))),1e-12)

    def test_kp_velocity_is_exact_derivative_of_kp_H(self):
        point=self.center+np.array([.04,-.07]);eps=1e-5
        for axis in range(2):
            delta=np.eye(2)[axis]*eps
            fd=(self.kp.matrices(point+delta)[0]-self.kp.matrices(point-delta)[0])/(2*eps)
            self.assertLess(np.max(abs(fd-self.kp.matrices(point,axis=axis)[0])),1e-10)

    def test_hessian_robust_to_finite_difference_step(self):
        kp2=ValleyKP.from_microscopic(self.obj,self.center,fd=5e-4)
        self.assertLess(np.max(abs(self.kp.h2-kp2.h2)),1e-8)

    def test_exact_kernel_matches_previous_implementation_on_same_disk(self):
        g=self.obj.model().reciprocal;centers=np.array([[1/3,2/3],[2/3,1/3]])@g
        k,w,meta=valley_quadrature(g,exterior_n=7,radius=.2,radial_n=10,angular_n=24)
        n=meta['n_valley_total']//2;k=k[-2*n:];w=w[-2*n:]
        older=static_inplane_matrix(self.obj.spec,k,w)['gamma']
        newer=compare_patch(self.obj.spec,k,w,centers,n,order=2)
        self.assertLess(np.max(abs(older-newer['full'])),5e-13)
        self.assertLess(abs(newer['kp'][0,0]-newer['full'][0,0]),1e-9)

    def test_no_fictitious_gap_closing_from_kp(self):
        data=json.loads((ROOT/'results/sliding_stage7_capture.json').read_text())
        self.assertTrue(all(min(row['min_mu_distance_eV'])>0 for row in data))

    def test_valley_capture_and_second_order_accuracy_frozen(self):
        data=json.loads((ROOT/'results/sliding_stage7_capture.json').read_text())
        selected=[x for x in data if x['radius_inv_nm']==.75]
        self.assertEqual(len(selected),4)
        self.assertTrue(all(abs(v['disk_fraction_signed_xx']-1)<.002 for v in selected))
        self.assertTrue(all(v['kp2_vs_disk_relative_error_xx']<1e-4 for v in selected))

    def test_valley_model_captures_interference_sign(self):
        data=json.loads((ROOT/'results/sliding_stage7_channels.json').read_text())
        self.assertEqual(len(data),3)
        for row in data:
            self.assertLess(row['BZ_mixed_coeff'],0)
            self.assertLess(row['disk_kp2_mixed_coeff'],0)
            self.assertGreater(row['BZ_even_coeff'],0)
            self.assertGreater(row['disk_kp2_even_coeff'],0)
            self.assertLess(abs(row['BZ_mixed_coeff']-row['disk_kp2_mixed_coeff']),.00012)

    def test_root_and_trigonal_scaling_frozen(self):
        data=json.loads((ROOT/'results/sliding_stage7_root_and_trace.json').read_text())
        self.assertLess(abs(data['estimated_zero_kp2']-.1974),.001)
        d={r['u']:r for r in data['trace_scan']}
        self.assertAlmostEqual(d[.1]['kp2_disk_trace']/d[.05]['kp2_disk_trace'],8,delta=.3)
        self.assertAlmostEqual(d[.1]['kp2_disk_nematic']/d[.05]['kp2_disk_nematic'],2,delta=.1)

    def test_model_parameter_robustness_is_finite_and_gapped(self):
        dat=json.loads((ROOT/'results/sliding_stage7_robust.json').read_text())['values']
        self.assertEqual(len(dat),7)
        for scenario in dat.values():
            first,last=scenario['data']
            self.assertLess(first['gamma_xz_x'],0)
            self.assertGreater(last['gamma_xz_x'],0)
            self.assertGreater(min(first['min_mu_distance_eV'],last['min_mu_distance_eV']),.002)

if __name__=='__main__':unittest.main()
