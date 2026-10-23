"""Round25 independent static Ward k-dependent spin/Bloch-gauge validation."""
import unittest
import numpy as np
from flux_lab.covariant_spin_ward import (
    periodic_spin_embedding,covariant_static_ward,covariant_velocities,
    fermi_frechet,finite_zeeman_current_derivative,finite_k_spin_density_gradient)
from flux_lab.degenerate_static_ward import degenerate_static_ward
from flux_lab.hexagonal_bz import BASIS

K=np.array([.36,-.29])
def run(k=K,variation=.16,rotate=True):
    H,V,S,dS,A=periodic_spin_embedding(k,effective_spin_variation=variation,rotate=rotate)
    return covariant_static_ward(H,V,S,dS,A), (H,V,S,dS,A)

class CovariantSpinWard(unittest.TestCase):
    def test_fixed_spin_reference_recovers_exact_kramers_static_ward(self):
        d,(H,V,S,dS,A)=run(np.array([0.,0.]),variation=0.,rotate=False)
        old=degenerate_static_ward(H,V,S,mu=.15,temperature=.22,kappa=.7)
        self.assertEqual(old["cluster_dimensions"],[2,2])
        np.testing.assert_allclose(d["full_static"],old["full_static"],atol=3e-13)
        np.testing.assert_allclose(d["spin_vertex_contact"],0.,atol=1e-13)

    def test_physical_full_response_is_gauge_invariant(self):
        for k in (K,np.array([0.,0.]),np.array([.71,.23])):
            fixed,_=run(k,variation=.16,rotate=False)
            rotated,_=run(k,variation=.16,rotate=True)
            for key in ("full_static","occupation","spin_vertex_contact","spin_density"):
                np.testing.assert_allclose(fixed[key],rotated[key],atol=3e-12,rtol=3e-12)

    def test_independent_finite_B_current_is_complete_with_spin_contact(self):
        for k in (K,np.array([0.,0.])):
            d,(H,V,S,dS,A)=run(k,variation=.16,rotate=True)
            finite=finite_zeeman_current_derivative(H,V,S,dS,A,
                                 mu=.15,temperature=.22,kappa=.7,step=1e-5)
            np.testing.assert_allclose(d["full_static"],finite,atol=5e-10)

    def test_independent_momentum_gradient_matches_complete_static_response(self):
        for k in (K,np.array([0.,0.])):
            d,_=run(k,variation=.16,rotate=True)
            finite=finite_k_spin_density_gradient(k,effective_spin_variation=.16,
                           rotate=True,mu=.15,temperature=.22,kappa=.7,step=1e-5)
            np.testing.assert_allclose(d["full_static"],finite,atol=6e-10)

    def test_nonzero_spin_vertex_contact_is_physically_required(self):
        d,(H,V,S,dS,A)=run(K,variation=.16,rotate=True)
        self.assertGreater(np.max(abs(d["spin_vertex_contact"])),1e-4)
        finite=finite_k_spin_density_gradient(K,effective_spin_variation=.16)
        np.testing.assert_allclose(d["full_static"],finite,atol=6e-10)
        self.assertGreater(np.max(abs(d["occupation"]-finite)),1e-4)

    def test_connection_changes_velocity_vertex_without_changing_observable(self):
        d,(H,V,S,dS,A)=run(K,variation=0.,rotate=True)
        self.assertGreater(np.max(abs(d["covariant_velocity"]-V)),1e-3)
        unrot,_=run(K,variation=0.,rotate=False)
        np.testing.assert_allclose(d["full_static"],unrot["full_static"],atol=2e-12)
        np.testing.assert_allclose(d["covariant_spin_gradient"],0,atol=1e-13)

    def test_frechet_kernel_is_unitary_covariant_at_kramers_point(self):
        H,V,S,dS,A=periodic_spin_embedding([0,0],rotate=False, effective_spin_variation=0.)
        rng=np.random.default_rng(25)
        Z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        U,_=np.linalg.qr(Z)
        r=lambda X:U.conj().T@X@U
        direct=fermi_frechet(H,S[0],mu=.15,temperature=.22)
        after=fermi_frechet(r(H),r(S[0]),mu=.15,temperature=.22)
        np.testing.assert_allclose(after,r(direct),atol=1e-12)

    def test_reciprocal_periodic_gauge_spin_and_ward(self):
        k=K
        ref,_=run(k,variation=.16,rotate=True)
        for dn in ([1,0],[0,1],[1,-1]):
            shifted=k+2*np.pi*BASIS@np.asarray(dn,float)
            got,_=run(shifted,variation=.16,rotate=True)
            np.testing.assert_allclose(got["full_static"],ref["full_static"],atol=2e-12)
            np.testing.assert_allclose(got["spin_density"],ref["spin_density"],atol=2e-12)

    def test_temperature_variation_does_not_change_covariant_algebra(self):
        H,V,S,dS,A=periodic_spin_embedding(K)
        for T in (.08,.17,.4):
            d=covariant_static_ward(H,V,S,dS,A,temperature=T)
            finite=finite_zeeman_current_derivative(H,V,S,dS,A,temperature=T)
            np.testing.assert_allclose(d["full_static"],finite,atol=6e-10)

    def test_bad_spin_derivative_and_connection_fail_closed(self):
        H,V,S,dS,A=periodic_spin_embedding(K)
        with self.assertRaisesRegex(ValueError,"dS="):
            covariant_static_ward(H,V,S,np.zeros((2,4,4)),A)
        with self.assertRaisesRegex(ValueError,"anti-Hermitian"):
            covariant_static_ward(H,V,S,dS,np.ones((2,4,4)))
        with self.assertRaisesRegex(ValueError,"T>0"):
            covariant_static_ward(H,V,S,dS,A,temperature=0)
        with self.assertRaises(ValueError):
            periodic_spin_embedding([0.,0.,0.])

if __name__=="__main__":unittest.main()
