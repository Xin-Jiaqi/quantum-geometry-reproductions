"""Round24 independent non-Abelian finite-T static Zeeman Ward controls."""
import unittest
import numpy as np
from flux_lab.core import PAULI
from flux_lab.equilibrium_ward import static_zeeman_current_ward,spin_density
from flux_lab.hexagonal_bz import BASIS,hexagonal_symmetrized_model
from flux_lab.degenerate_static_ward import (
    fermi_divided_difference,degenerate_static_ward,projector_static_ward,
    kramers_periodic_model,finite_zeeman_static,finite_k_spin_gradient)

K=np.array([.31,-.22])
PARAMS=dict(mu=.15,temperature=.22,kappa=.7)

def random_unitary(n,seed):
    rng=np.random.default_rng(seed)
    Z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    return np.linalg.qr(Z)[0]

def rotate(H,V,S,U):
    cv=lambda A:U@A@U.conj().T
    return cv(H),np.stack([cv(x) for x in V]),np.stack([cv(x) for x in S])

class DegenerateStaticWard(unittest.TestCase):
    def setUp(self):
        self.H,self.V,self.S=kramers_periodic_model(np.zeros(2))

    def test_exact_kramers_doublet_spectrum(self):
        E=np.linalg.eigvalsh(self.H)
        self.assertLess(abs(E[0]-E[1]),1e-13)
        self.assertLess(abs(E[2]-E[3]),1e-13)
        D=degenerate_static_ward(self.H,self.V,self.S,**PARAMS)
        self.assertEqual(D["cluster_dimensions"],[2,2])
        self.assertGreater(np.max(abs(D["intracluster"])),1e-3)

    def test_time_reversal_acts_on_full_H_velocity_and_spin(self):
        UT=np.kron(np.eye(2),1j*PAULI[1])
        T=lambda A:UT@A.conj()@UT.conj().T
        np.testing.assert_allclose(T(self.H),self.H,atol=3e-14)
        for v in self.V:np.testing.assert_allclose(T(v),-v,atol=3e-14)
        for s in self.S:np.testing.assert_allclose(T(s),-s,atol=3e-14)

    def test_divided_difference_at_exact_repeated_energy_matches_fermi_prime(self):
        E=np.array([-1.25,-1.25,.43,.43])
        q=fermi_divided_difference(E,mu=.15,temperature=.22)
        f=1/(1+np.exp((E-.15)/.22))
        prime=-f*(1-f)/.22
        self.assertAlmostEqual(q[0,1],prime[0],delta=2e-14)
        self.assertAlmostEqual(q[2,3],prime[2],delta=2e-14)
        self.assertTrue(np.all(np.isfinite(q)))
        np.testing.assert_allclose(q,q.T,atol=0)

    def test_near_degenerate_hyperbolic_formula_is_continuous(self):
        E=np.array([-1.2,-1.2+1e-10,.7,.7+1e-10])
        a=fermi_divided_difference(E,mu=.15,temperature=.22)
        b=fermi_divided_difference(np.array([-1.2,-1.2,.7,.7]),mu=.15,temperature=.22)
        np.testing.assert_allclose(a,b,atol=2e-10,rtol=0)

    def test_full_static_matches_field_and_spin_momentum_finite_differences_at_TRIM(self):
        r=degenerate_static_ward(self.H,self.V,self.S,**PARAMS)
        a=finite_zeeman_static(self.H,self.V,self.S,**PARAMS,step=1e-5)
        b=finite_k_spin_gradient([0.,0.],**PARAMS,step=1e-5)
        np.testing.assert_allclose(r["full_static"],a,atol=3e-9)
        np.testing.assert_allclose(r["full_static"],b,atol=3e-9)

    def test_projector_contraction_matches_band_pair_at_exact_degeneracy(self):
        q=degenerate_static_ward(self.H,self.V,self.S,**PARAMS)
        p=projector_static_ward(self.H,self.V,self.S,**PARAMS)
        for x in ("full_static","intracluster","intercluster"):
            np.testing.assert_allclose(p[x],q[x],atol=3e-12)

    def test_entire_result_invariant_under_arbitrary_orbital_spin_unitary(self):
        u=random_unitary(4,202)
        r=degenerate_static_ward(self.H,self.V,self.S,**PARAMS)
        x=degenerate_static_ward(*rotate(self.H,self.V,self.S,u),**PARAMS)
        for key in ("full_static","intracluster","intercluster"):
            np.testing.assert_allclose(r[key],x[key],atol=3e-12)

    def test_within_doublet_rotation_changes_naive_diagonal_occupation_only(self):
        e,U=np.linalg.eigh(self.H)
        r=random_unitary(2,30);s=random_unitary(2,31)
        W=np.zeros((4,4),complex);W[:2,:2]=r;W[2:,2:]=s
        U1=U@W
        # Explicit manual diagonal-only expectation, no arbitrarily assigned
        # orbital axes allowed as a physical degenerate-sector response.
        def naive(X):
            v=np.stack([X.conj().T@v@X for v in self.V])
            sp=np.stack([X.conj().T@sp@X for sp in self.S])
            fp=-1/(4*PARAMS["temperature"]*np.cosh((e-PARAMS["mu"])/(2*PARAMS["temperature"]))**2)
            return -PARAMS["kappa"]*np.einsum("an,bn,n->ab",np.diagonal(v,axis1=-2,axis2=-1),
                    np.diagonal(sp,axis1=-2,axis2=-1),fp)
        self.assertGreater(np.max(abs(naive(U)-naive(U1))),1e-4)
        np.testing.assert_allclose(degenerate_static_ward(self.H,self.V,self.S,**PARAMS)["full_static"],
                                   degenerate_static_ward(self.H,self.V,self.S,**PARAMS)["intracluster"]+
                                   degenerate_static_ward(self.H,self.V,self.S,**PARAMS)["intercluster"],atol=1e-12)

    def test_old_nondegenerate_formula_agrees_away_from_TRIM(self):
        H,V,S=hexagonal_symmetrized_model([.47,-.31],"6")
        old=static_zeeman_current_ward(H,V,S,mu=2.3,temperature=.3)
        new=degenerate_static_ward(H,V,S,mu=2.3,temperature=.3)
        np.testing.assert_allclose(new["full_static"],old["full_static"],atol=2e-12)
        np.testing.assert_allclose(new["intracluster"],old["occupation_derivative"],atol=2e-12)

    def test_old_nondegenerate_formula_rejects_TRIM_but_new_succeeds(self):
        with self.assertRaisesRegex(ValueError,"nondegenerate"):
            static_zeeman_current_ward(self.H,self.V,self.S,**PARAMS)
        d=degenerate_static_ward(self.H,self.V,self.S,**PARAMS)
        self.assertTrue(np.all(np.isfinite(d["full_static"])))

    def test_analytic_cartesian_velocities_match_independent_finite_k_derivative(self):
        k=K;H,V,S=kramers_periodic_model(k)
        for a in range(2):
            dk=np.zeros(2);dk[a]=2e-5
            H1=kramers_periodic_model(k+dk)[0];H2=kramers_periodic_model(k-dk)[0]
            np.testing.assert_allclose(V[a],(H1-H2)/(4e-5),atol=2e-10)

    def test_periodic_reciprocal_lattice_and_time_reversal_at_generic_k(self):
        k=K
        H,V,S=kramers_periodic_model(k)
        for direction in (np.array([2*np.pi,0]),np.array([0,2*np.pi])):
            x=kramers_periodic_model(k+BASIS@direction)
            np.testing.assert_allclose(x[0],H,atol=2e-12)
            np.testing.assert_allclose(x[1],V,atol=2e-12)
        UT=np.kron(np.eye(2),1j*PAULI[1])
        H2,V2,S2=kramers_periodic_model(-k)
        np.testing.assert_allclose(H2,UT@H.conj()@UT.conj().T,atol=2e-13)

    def test_energy_shift_with_matching_mu_preserves_static_response(self):
        d=degenerate_static_ward(self.H,self.V,self.S,**PARAMS)
        shifted=degenerate_static_ward(self.H+3*np.eye(4),self.V,self.S,
                            mu=PARAMS["mu"]+3,temperature=PARAMS["temperature"],kappa=.7)
        np.testing.assert_allclose(d["full_static"],shifted["full_static"],atol=2e-13)

    def test_input_validation_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"strictly positive"):
            degenerate_static_ward(self.H,self.V,self.S,temperature=0)
        with self.assertRaises(ValueError):
            degenerate_static_ward(self.H,self.V,self.S[:2])
        with self.assertRaises(ValueError):
            degenerate_static_ward(self.H+1j*np.eye(4),self.V,self.S)
        with self.assertRaises(ValueError):
            degenerate_static_ward(self.H,self.V,self.S,degeneracy_tolerance=-1)

    def test_round24_frozen_audit_static_kramers_and_generic_points(self):
        import json
        from pathlib import Path
        root=Path(__file__).resolve().parents[1]
        ref=json.loads((root/"results/frozen/round24_degenerate_static_ward_audit.json").read_text())
        self.assertEqual(ref["schema"],"cui-round24-degenerate-static-ward-v1")
        for case in ref["results"]:
            k=np.array(case["k"])
            H,V,S=kramers_periodic_model(k)
            actual=degenerate_static_ward(H,V,S,**PARAMS)
            frozen=np.array(case["full"]["real"])+1j*np.array(case["full"]["imag"])
            np.testing.assert_allclose(actual["full_static"],frozen,atol=2e-12,rtol=0)
            self.assertEqual(actual["cluster_dimensions"],case["cluster_dimensions"])

if __name__=="__main__":unittest.main()
