"""Round28 independent charge-conserving toy relaxation and gauge tests."""
import unittest
import numpy as np
from flux_lab.finite_q_gauge_ward import soc_torus,site_phase,pure_gauge_response
from flux_lab.conserving_collision import (
    local_block,collision_projector,pure_gauge_drive,zeeman_vertex,collision_response)

Z=.42+.13j
RATE=.12
def setup():
    t=soc_torus()
    L=site_phase(t,qx=1,qy=1,offset=.21)
    K,Vscalar=pure_gauge_drive(t,L,Z)
    vals=np.diagonal(L).real[::2]
    Vz=zeeman_vertex(t,vals,spin_axis="y",kappa=.7)
    return t,L,K,Vscalar,Vz

class LocalConservingCollision(unittest.TestCase):
    def test_site_collision_projection_conserves_charge_exactly(self):
        t,*_=setup()
        rng=np.random.default_rng(28)
        X=rng.normal(size=t.H.shape)+1j*rng.normal(size=t.H.shape)
        Q=collision_projector(X,t.projectors)
        np.testing.assert_allclose(np.einsum('aij,ji->a',t.projectors,Q),0,atol=2e-14)
        np.testing.assert_allclose(collision_projector(Q,t.projectors),Q,atol=2e-14)
        np.testing.assert_allclose(local_block(Q,t.projectors),0,atol=2e-14)

    def test_collision_projector_local_gauge_covariant(self):
        t,L,*_=setup()
        rng=np.random.default_rng(29)
        X=rng.normal(size=t.H.shape)+1j*rng.normal(size=t.H.shape)
        U=np.diag(np.exp(1j*np.diagonal(L)))
        rot=lambda M:U@M@U.conj().T
        np.testing.assert_allclose(collision_projector(rot(X),t.projectors),
                                   rot(collision_projector(X,t.projectors)),atol=5e-14)

    def test_collisionless_limit_exactly_recovers_round27_pure_gauge(self):
        t,L,K,Vs,_=setup()
        old=pure_gauge_response(t,L,z=Z)
        got=collision_response(t,K+Vs,G=L,gamma=0,z=Z)
        np.testing.assert_allclose(got["dr"],old["dr_gauge"],atol=4e-13)
        np.testing.assert_allclose(got["current"],old["current_gauge"],atol=4e-13)

    def test_nonzero_gamma_pure_gauge_remains_exact(self):
        t,L,K,Vs,_=setup()
        for gamma in (.03,.12,.25):
            out=collision_response(t,K+Vs,G=L,gamma=gamma,z=Z)
            expected=1j*(L@out["rho"]-out["rho"]@L)
            np.testing.assert_allclose(out["dr"],expected,atol=3e-12)
            np.testing.assert_allclose(out["density"],0,atol=3e-12)
            np.testing.assert_allclose(out["current"],0,atol=3e-12)
            np.testing.assert_allclose(out["collision"],0,atol=3e-12)
            np.testing.assert_allclose(out["ward"],0,atol=3e-12)

    def test_vector_and_scalar_separate_nonzero_sum_pure_gauge(self):
        t,L,K,Vs,_=setup()
        vec=collision_response(t,K,G=L,z=Z)
        scalar=collision_response(t,Vs,z=Z)
        full=collision_response(t,K+Vs,G=L,z=Z)
        self.assertGreater(np.max(abs(vec["density"])),.04)
        self.assertGreater(np.max(abs(scalar["density"])),.04)
        np.testing.assert_allclose(vec["dr"]+scalar["dr"],full["dr"],atol=2e-12)
        np.testing.assert_allclose(vec["current"]+scalar["current"],full["current"],atol=2e-12)

    def test_true_longitudinal_scalar_response_nonzero_and_conserved(self):
        t,L,_,_,_=setup()
        out=collision_response(t,.17*L,z=Z)
        self.assertGreater(np.max(abs(out["density"])),.004)
        self.assertGreater(np.max(abs(out["current"])),.005)
        np.testing.assert_allclose(out["ward"],0,atol=3e-12)
        np.testing.assert_allclose(out["collision_site"],0,atol=5e-14)

    def test_nontrivial_collision_rate_changes_scalar_response(self):
        t,L,*_=setup()
        baseline=collision_response(t,.17*L,gamma=0,z=Z)
        dressed=collision_response(t,.17*L,gamma=.23,z=Z)
        self.assertGreater(np.max(abs(baseline["density"]-dressed["density"])),1e-4)
        np.testing.assert_allclose(dressed["ward"],0,atol=4e-12)

    def test_physical_spin_zeeman_drive_nonzero_and_continuous(self):
        t,_,_,_,Vz=setup()
        out=collision_response(t,Vz,z=Z)
        self.assertGreater(np.max(abs(out["dr"])),.005)
        np.testing.assert_allclose(out["ward"],0,atol=3e-12)
        np.testing.assert_allclose(out["collision_site"],0,atol=3e-14)
        self.assertGreater(np.max(abs(t.H@Vz-Vz@t.H)),.01)

    def test_zeeman_current_invariant_under_added_em_pure_gauge(self):
        t,L,K,Vs,Vz=setup()
        plain=collision_response(t,Vz,z=Z)
        transformed=collision_response(t,Vz+K+Vs,G=L,z=Z)
        expected=1j*(L@plain["rho"]-plain["rho"]@L)
        np.testing.assert_allclose(transformed["dr"]-plain["dr"],expected,atol=2e-12)
        np.testing.assert_allclose(transformed["density"],plain["density"],atol=2e-12)
        np.testing.assert_allclose(transformed["current"],plain["current"],atol=2e-12)

    def test_naive_global_rta_fails_sitewise_charge_continuity(self):
        t,L,_,_,_=setup()
        bad=collision_response(t,.17*L,z=Z,conserve_local_charge=False,gauge_reference=False)
        good=collision_response(t,.17*L,z=Z)
        self.assertGreater(np.max(abs(bad["ward"])),1e-3)
        np.testing.assert_allclose(good["ward"],0,atol=2e-12)

    def test_ungauged_reference_fails_pure_gauge_while_conserving_charge(self):
        t,L,K,Vs,_=setup()
        bad=collision_response(t,K+Vs,G=L,z=Z,gauge_reference=False)
        self.assertGreater(np.max(abs(bad["dr"]-1j*(L@bad["rho"]-bad["rho"]@L))),1e-3)
        self.assertGreater(np.max(abs(bad["current"])),1e-3)
        np.testing.assert_allclose(bad["collision_site"],0,atol=3e-14)
        np.testing.assert_allclose(bad["ward"],0,atol=3e-12)

    def test_independent_dense_liouville_agrees_with_spectral_iteration(self):
        t,L,K,Vs,Vz=setup()
        for V,G in ((.17*L,None),(Vz+K+Vs,L)):
            a=collision_response(t,V,G=G,z=Z,method="spectral")
            b=collision_response(t,V,G=G,z=Z,method="liouville")
            np.testing.assert_allclose(a["dr"],b["dr"],atol=4e-12)
            np.testing.assert_allclose(a["current"],b["current"],atol=4e-12)

    def test_bad_parameters_and_pseudospin_not_accepted(self):
        t,L,_,_,_=setup()
        for bad_gamma in (-.01,float("nan")):
            with self.assertRaises(ValueError):collision_response(t,L,gamma=bad_gamma)
        with self.assertRaisesRegex(ValueError,"retarded"):collision_response(t,L,z=.4)
        with self.assertRaisesRegex(ValueError,"Peierls"):collision_response(t,L,G=np.ones(t.H.shape))
        with self.assertRaisesRegex(ValueError,"method"):collision_response(t,L,method="wrong")
        with self.assertRaisesRegex(ValueError,"physical spin axis"):
            zeeman_vertex(t,np.ones(9),spin_axis='orbital')
        with self.assertRaises(ValueError):zeeman_vertex(t,np.ones(8))

if __name__=="__main__":unittest.main()
