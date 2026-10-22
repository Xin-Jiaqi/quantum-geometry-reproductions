"""Round23 independently check static equilibrium Ward and missing Fermi term."""
import unittest
import numpy as np
from flux_lab.core import PAULI
from flux_lab.causal_response import interband_kubo
from flux_lab.hexagonal_bz import hexagonal_symmetrized_model
from flux_lab.equilibrium_ward import (
    static_zeeman_current_ward,equilibrium_current,
    zeeman_finite_difference,spin_k_gradient_fd,static_grid_audit)

I=np.eye(2,dtype=complex)
K=np.array([.47,-.31])

def smooth_two_band(k):
    x,y=k
    d=np.array([.8+.15*np.cos(x),.55*np.sin(y),.7*np.sin(x)+.1*np.cos(y)])
    dv=np.array([[-.15*np.sin(x),0.,.7*np.cos(x)],
                 [0.,.55*np.cos(y),-.1*np.sin(y)]])
    c=.2*np.sin(x)+.08*np.cos(y)
    dc=np.array([.2*np.cos(x),-.08*np.sin(y)])
    H=c*I+np.einsum("i,ijk->jk",d,PAULI)
    V=np.stack([dc[a]*I+np.einsum("i,ijk->jk",dv[a],PAULI) for a in range(2)])
    return H,V,PAULI.copy()

def gauge_transform(data,U):
    H,V,S=data
    t=lambda x:U@x@U.conj().T
    return t(H),np.stack([t(v) for v in V]),np.stack([t(s) for s in S])

class StaticEquilibriumWard(unittest.TestCase):
    def test_interband_static_matches_existing_directed_kubo(self):
        for model in (smooth_two_band,lambda k:hexagonal_symmetrized_model(k,"1")):
            H,V,S=model(K)
            r=static_zeeman_current_ward(H,V,S,mu=.42,temperature=.28)
            ref=interband_kubo(H,V,S,z=0.,mu=.42,temperature=.28)[0]
            np.testing.assert_allclose(r["interband"],ref,atol=2e-12,rtol=1e-11)

    def test_full_static_matches_independent_finite_B_current(self):
        for model in (smooth_two_band,lambda k:hexagonal_symmetrized_model(k,"1")):
            H,V,S=model(K)
            r=static_zeeman_current_ward(H,V,S,mu=.6,temperature=.3)
            fd=zeeman_finite_difference(H,V,S,mu=.6,temperature=.3,field_step=1e-5)
            np.testing.assert_allclose(r["full_static"],fd,atol=3e-9,rtol=1e-7)

    def test_static_derivative_equals_negative_spin_density_gradient(self):
        for model in (smooth_two_band,lambda k:hexagonal_symmetrized_model(k,"1")):
            H,V,S=model(K)
            r=static_zeeman_current_ward(H,V,S,mu=.6,temperature=.3)
            fd=spin_k_gradient_fd(model,K,mu=.6,temperature=.3,k_step=1e-5)
            np.testing.assert_allclose(r["full_static"],fd,atol=5e-9,rtol=3e-7)

    def test_fermi_occupation_derivative_is_nonzero_and_necessary(self):
        H,V,S=smooth_two_band(K)
        r=static_zeeman_current_ward(H,V,S,mu=.42,temperature=.3)
        self.assertGreater(np.max(abs(r["occupation_derivative"])),1e-4)
        fd=zeeman_finite_difference(H,V,S,mu=.42,temperature=.3)
        self.assertLess(np.max(abs(r["full_static"]-fd)),1e-8)
        self.assertGreater(np.max(abs(r["interband"]-fd)),1e-4)

    def test_low_temperature_occupation_derivative_suppressed_in_gap(self):
        H=1.1*PAULI[2]+.21*PAULI[0]
        V=np.stack([PAULI[0],PAULI[1]])
        S=PAULI.copy()
        r=static_zeeman_current_ward(H,V,S,mu=0.,temperature=.015)
        self.assertLess(np.max(abs(r["occupation_derivative"])),1e-20)

    def test_all_static_terms_invariant_under_constant_orbital_basis_change(self):
        H,V,S=smooth_two_band(K)
        rng=np.random.default_rng(72)
        X=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        U,_=np.linalg.qr(X)
        a=static_zeeman_current_ward(H,V,S,mu=.42,temperature=.3)
        b=static_zeeman_current_ward(*gauge_transform((H,V,S),U),
                                     mu=.42,temperature=.3)
        for key in ("interband","occupation_derivative","full_static"):
            np.testing.assert_allclose(a[key],b[key],atol=2e-12)

    def test_finite_B_current_is_real_in_equilibrium(self):
        H,V,S=smooth_two_band(K)
        for b in (.01,-.01):
            j=equilibrium_current(H-b*S[2],V,mu=.42,temperature=.3)
            self.assertLess(np.max(abs(j.imag)),2e-14)

    def test_periodic_grid_telescope_vanishes_exactly(self):
        a=static_grid_audit("1",order=4,mu=2.3,temperature=.3)
        self.assertEqual(a["nodes"],16)
        self.assertLess(np.max(abs(a["telescoping_derivative"])),1e-12)
        np.testing.assert_allclose(a["interband"]+a["occupation_derivative"],
                                   a["full_static"],atol=1e-12)

    def test_finite_grid_spectral_current_need_not_be_zero(self):
        # No false claim: telescoping FD is exactly zero at every mesh
        # while sampled spectral derivative generally is not.
        a=static_grid_audit("1",order=4,mu=2.3,temperature=.3)
        self.assertTrue(np.all(np.isfinite(a["full_static"])))
        self.assertIn("NOT bounded",a["status"])

    def test_nonconstant_spin_or_bad_inputs_do_not_silently_pass(self):
        H,V,S=smooth_two_band(K)
        with self.assertRaises(ValueError):
            static_zeeman_current_ward(H,V,S,temperature=0.)
        with self.assertRaises(ValueError):
            static_zeeman_current_ward(np.zeros((2,2)),V,S)
        with self.assertRaises(ValueError):
            static_zeeman_current_ward(H,V,S[:2])
        with self.assertRaises(ValueError):
            zeeman_finite_difference(H,V,S,field_step=0.)
        with self.assertRaises(ValueError):
            static_grid_audit("1",order=4,method="gauss")

if __name__=="__main__":unittest.main()
