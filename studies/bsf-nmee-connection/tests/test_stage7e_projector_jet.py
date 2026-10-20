"""Physical gates for exact weak-t projector response coefficients."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingSpec, SlidingBilayer
from bsf_nmee.projector_jet_response import (
    analytic_valley_ABC, projector_coefficient)


class Stage7EProjectorPhysicsTests(unittest.TestCase):
    def test_symmetry_forbidden_at_reference_slide(self):
        spec=replace(SlidingSpec(),shift=(0.,0.))
        response=analytic_valley_ABC(spec,nr=10,nt=24)
        for key in ("A_same2","B_cross2","C_mixed"):
            self.assertLess(abs(response["coefficients"][key]["total"]),3e-7)

    def test_velocity_pathways_have_opposite_physical_roles(self):
        u=.4
        spec=replace(SlidingSpec(),shift=(.23*u,.41*u))
        r=analytic_valley_ABC(spec,nr=14,nt=36)["coefficients"]
        A,B,C=(r["A_same2"],r["B_cross2"],r["C_mixed"])
        self.assertGreater(A["total"],0)
        self.assertGreater(B["total"],0)
        self.assertLess(C["total"],0)
        self.assertGreater(B["v0"]/B["total"],.99)
        self.assertGreater(C["vy"]/C["total"],.65)
        self.assertLess(C["v0"],0)
        self.assertGreater(C["vx"],0)
        self.assertAlmostEqual(
            A["total"]+B["total"]+C["total"],
            sum(x["total"] for x in (A,B,C)),places=12)

    def test_global_unitary_covariance(self):
        # No k-dependent reembedding: a constant unitary is an exact basis change.
        spec=replace(SlidingSpec(),shift=(.23*.4,.41*.4))
        no_h=SlidingBilayer(replace(spec,tunneling_eV=0.))
        same=SlidingBilayer(replace(spec,tunneling_eV=1.,cross_sublattice=False))
        full=SlidingBilayer(replace(spec,tunneling_eV=1.,cross_sublattice=True))
        reciprocal=no_h.model().reciprocal
        k=np.array([[1/3,2/3]])@reciprocal+np.array([[.05,.02]])
        h0=no_h.matrices(k)
        hs=same.matrices(k)
        hf=full.matrices(k)
        vx=no_h.matrices(k,axis=0)
        wx=same.matrices(k,axis=0)-vx
        wy=full.matrices(k,axis=0)-same.matrices(k,axis=0)
        Sx=no_h.spins["total_x"];Sz=no_h.spins["total_z"]
        original=projector_coefficient(h0,hs-h0,hf-hs,vx,wx,wy,Sx,Sz,spec.mu_eV,
                                       g_factor=no_h.mono.g_factor)
        rng=np.random.default_rng(128)
        z=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8))
        Q,_=np.linalg.qr(z)
        transform=lambda m:Q.conj().T@m@Q
        rotated=projector_coefficient(*[transform(m) for m in
            (h0,hs-h0,hf-hs,vx,wx,wy,Sx,Sz)],
            mu=spec.mu_eV,g_factor=no_h.mono.g_factor)
        np.testing.assert_allclose(original,rotated,rtol=1e-8,atol=1e-8)

    def test_analytic_coefficient_matches_weak_hopping(self):
        # Exact direct response at finite hopping compared with analytic
        # projector derivative. No finite-difference sources in analytic code.
        from bsf_nmee.channel_spectral_mechanism import disk_channel_sectors
        u=.4
        spec=replace(SlidingSpec(),shift=(.23*u,.41*u),
                     tunneling_eV=.000625)
        radius=.65;nr=14;nt=36
        coefficient=analytic_valley_ABC(spec,radius=radius,nr=nr,nt=nt)
        c=coefficient["coefficients"]["C_mixed"]["total"]
        plus=disk_channel_sectors(spec,radius=radius,radial_n=nr,angular_n=nt,
                                 same_scale=1.,cross_scale=1.)["total"]
        minus=disk_channel_sectors(spec,radius=radius,radial_n=nr,angular_n=nt,
                                  same_scale=1.,cross_scale=-1.)["total"]
        finite=(plus-minus)/(2*spec.tunneling_eV**2)
        self.assertLess(abs(finite-c)/abs(c),.005)


if __name__=="__main__":
    unittest.main()
