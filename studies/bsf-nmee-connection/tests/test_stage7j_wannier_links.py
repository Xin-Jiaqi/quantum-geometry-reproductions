"""Stage7j: gauge stress tests for future spinor Wannier operator pipelines."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingSpec,SlidingBilayer
from bsf_nmee.wannier_gauge_links import (
    embedding_unitary,physical_link,gauge_stress_case)


def spec_for(u,t):
    return replace(SlidingSpec(),shift=(.23*u,.41*u),tunneling_eV=t)


class Stage7JGaugeLinksTests(unittest.TestCase):
    def test_pure_reembedding_preserves_static_kernel(self):
        for u,t in ((.1,.005),(.4,.005),(.4,.02)):
            spec=spec_for(u,t)
            k=np.array([1/3,2/3])@SlidingBilayer(spec).model().reciprocal
            k=k+np.array([.05,.02])
            row=gauge_stress_case(spec,k)
            self.assertLess(row['physical_gauge_error'],1e-6)
            self.assertGreater(abs(row['naive_artifact']),.01)
            self.assertGreater(row['minimum_mu_gap_eV'],.003)
            for v in row['derivative_errors'].values():
                self.assertLess(v['covariant_error'],2e-7)
                self.assertGreater(v['naive_error'],.01)

    def test_naive_parameter_gauge_can_reverse_local_sign(self):
        spec=spec_for(.1,.005)
        k=np.array([1/3,2/3])@SlidingBilayer(spec).model().reciprocal
        row=gauge_stress_case(spec,k+np.array([.05,.02]))
        self.assertLess(row['original'],0)
        self.assertLess(row['corrected'],0)
        self.assertGreater(row['naive'],0)

    def test_link_is_invariant_under_independent_U4_frame_rotation(self):
        rng=np.random.default_rng(20261020)
        Ua=embedding_unitary(np.array([.13,.07]),+.0001)
        Ub=embedding_unitary(np.array([.19,.04]),-.0001)
        Y=[]
        for i in range(2):
            Q,_=np.linalg.qr(rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)))
            Y.append(Q[:,:4])
        Qs=[]
        for i in range(2):
            Q,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
            Qs.append(Q)
        original=physical_link(Y[0],Y[1],Ua,Ub)
        transformed=physical_link(Y[0]@Qs[0],Y[1]@Qs[1],Ua,Ub)
        np.testing.assert_allclose(
            transformed,Qs[0].conj().T@original@Qs[1],atol=2e-14)

    def test_spins_and_source_derivatives_require_same_frame(self):
        spec=spec_for(.4,.005)
        model=SlidingBilayer(spec)
        k=np.array([1/3,2/3])@model.model().reciprocal+np.array([.05,.02])
        row=gauge_stress_case(spec,k)
        self.assertGreater(row['derivative_errors']['lambda']['naive_error'],1.)
        self.assertLess(row['derivative_errors']['lambda']['covariant_error'],1e-7)


if __name__=='__main__':
    unittest.main()
