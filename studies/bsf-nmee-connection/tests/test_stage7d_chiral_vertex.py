"""Stage 7d physical gates: chiral valley vertex and spin-rotation sign."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingBilayer, SlidingSpec
from bsf_nmee.valley_vertex_symmetry import (
    chiral_valley_summary, orbital_hopping_vertices, rashba_spin_rotation_check)
from bsf_nmee.channel_spectral_mechanism import disk_channel_sectors


class Stage7dPhysicsTests(unittest.TestCase):
    def test_valley_chiral_hopping_and_velocity_identity(self):
        spec=replace(SlidingSpec(),shift=(0.,0.),tunneling_eV=.005)
        vals=chiral_valley_summary(spec)
        self.assertEqual(len(vals),2)
        for i,r in enumerate(vals):
            self.assertLess(r['hopping_AB_meV'],1e-8)
            self.assertGreater(r['velocity_AB_x_meVnm'],.3)
            self.assertLess(r['slide_velocity_identity_residual'],1e-11)
            self.assertLess(r['microscopic_velocity_residual'],1e-11)
            self.assertAlmostEqual(r['vy_over_vx_AB'][1],1 if i==0 else -1,places=8)
            self.assertAlmostEqual(r['vy_over_vx_BA'][1],-1 if i==0 else 1,places=8)

    def test_slide_derivative_against_displaced_structures(self):
        spec=replace(SlidingSpec(),shift=(0.,0.),tunneling_eV=.005)
        obj=SlidingBilayer(spec)
        k=chiral_valley_summary(spec)[0]['center']
        _,_,J=orbital_hopping_vertices(obj,k)
        delta=1.e-5
        for axis in (0,1):
            d=np.eye(2)[axis]*delta
            shift_plus=d@np.linalg.inv(obj.a)
            shift_minus=-shift_plus
            plus=SlidingBilayer(replace(spec,shift=tuple(shift_plus))).tunneling(
                np.asarray(k)[None,:])[0][0,2]
            minus=SlidingBilayer(replace(spec,shift=tuple(shift_minus))).tunneling(
                np.asarray(k)[None,:])[0][0,2]
            self.assertLess(abs((plus-minus)/(2*delta)-J[axis,0,1]),1e-8)

    def test_rashba_spin_rotation(self):
        for u in (.1,.4):
            spec=replace(SlidingSpec(),shift=(u*.23,u*.41),
                         tunneling_eV=.005)
            g=SlidingBilayer(spec).model().reciprocal
            k=np.array([1/3,2/3])@g+np.array([.13,.07])
            r=rashba_spin_rotation_check(spec,k)
            self.assertLess(max(r.values()),1e-11)

    def test_tunable_channel_ratio_sign_bracket(self):
        spec=replace(SlidingSpec(),shift=(.4*.23,.4*.41),
                     tunneling_eV=.0025)
        values=[]
        for ratio in (.50,.57):
            r=disk_channel_sectors(spec,cross_scale=ratio,radius=.75,
                                  radial_n=14,angular_n=48)
            self.assertGreater(r['sampled_min_mu_distance_eV'],.002)
            values.append(r['total'])
        self.assertLess(values[0],0)
        self.assertGreater(values[1],0)


if __name__ == '__main__':
    unittest.main()
