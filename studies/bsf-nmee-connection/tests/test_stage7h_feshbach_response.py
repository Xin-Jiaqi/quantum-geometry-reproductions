"""Stage 7h: exact response-preserving Feshbach identity physics tests."""
from dataclasses import replace
import unittest
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingBilayer, SlidingSpec
from bsf_nmee.feshbach_response import (
    reference_basis,feshbach,inserted_source,feshbach_static_at_k)


class Stage7HExactResponseTests(unittest.TestCase):
    def test_complex_resolvent_and_all_source_contacts(self):
        for t in (.005,.020):
            spec=replace(SlidingSpec(),shift=(.4*.23,.4*.41),
                         tunneling_eV=t)
            obj=SlidingBilayer(spec)
            g=obj.model().reciprocal
            k=np.array([1/3,2/3])@g+np.array([.05,.02])
            frame=reference_basis(obj,k)
            H=obj.matrices(k[None,:])[0]
            vx=obj.matrices(k[None,:],axis=0)[0]
            for z in (.01+.02j,-.1+.04j):
                d=feshbach(z,H,frame)
                full=np.linalg.inv(z*np.eye(8)-H)
                self.assertLess(np.max(abs(full-d['G'])),5e-12)
                for O in (obj.spins['total_x'],obj.spins['total_z'],vx):
                    sectors,eff=inserted_source(d,O)
                    self.assertEqual(eff.shape,(4,4))
                    self.assertLess(np.max(abs(sum(sectors.values())
                                                -full@O@full)),3e-9)

    def test_contour_recovers_eigen_kubo_at_K(self):
        for u,t in ((.1,.005),(.4,.005),(.4,.02)):
            spec=replace(SlidingSpec(),shift=(u*.23,u*.41),
                         tunneling_eV=t)
            g=SlidingBilayer(spec).model().reciprocal
            for center in ([1/3,2/3],[2/3,1/3]):
                k=np.asarray(center)@g
                res=feshbach_static_at_k(spec,k,contour_n=48)
                self.assertLess(abs(res['gamma_eigen']-
                                    res['gamma_feshbach_exact']),1e-8)
                self.assertLess(max(v['projector_error'] for v in
                                    res['plus_minus']),1e-10)

    def test_missing_contacts_can_reverse_local_sign(self):
        spec=replace(SlidingSpec(),shift=(.4*.23,.4*.41),
                     tunneling_eV=.020)
        k=np.array([1/3,2/3])@SlidingBilayer(spec).model().reciprocal
        res=feshbach_static_at_k(spec,k,contour_n=48)
        self.assertLess(res['gamma_feshbach_exact'],0)
        self.assertGreater(res['gamma_dressed_only_diagnostic'],0)
        self.assertLess(res['gamma_missing_contacts'],0)
        self.assertGreater(abs(res['gamma_missing_contacts']),
                           abs(res['gamma_feshbach_exact']))

    def test_contour_converges_and_field_step_sign_is_stable(self):
        spec=replace(SlidingSpec(),shift=(.4*.23,.41*.4),
                     tunneling_eV=.020)
        k=np.array([1/3,2/3])@SlidingBilayer(spec).model().reciprocal
        coarse=feshbach_static_at_k(spec,k,contour_n=24)
        fine=feshbach_static_at_k(spec,k,contour_n=64)
        self.assertLess(abs(coarse['gamma_feshbach_exact']-
                            fine['gamma_feshbach_exact']),1e-8)
        for step in (.125,.5):
            r=feshbach_static_at_k(spec,k,step_T=step,contour_n=48)
            self.assertLess(r['gamma_feshbach_exact'],0)
            self.assertGreater(r['gamma_dressed_only_diagnostic'],0)


if __name__=='__main__':
    unittest.main()
