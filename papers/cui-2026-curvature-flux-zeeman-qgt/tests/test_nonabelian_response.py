"""Round17 non-Abelian cluster response and Kramers-gauge tests."""
import unittest
import numpy as np
from flux_lab.core import PAULI
from flux_lab.nonabelian_response import (
    cluster_structure,cluster_zeeman_qgt,cluster_interband_kubo,
    cluster_interband_liouville)

I=np.eye(2,dtype=complex)
TAU=np.stack([np.kron(t,I) for t in PAULI])
SPIN=np.stack([np.kron(I,s) for s in PAULI])

def kramers_model(split=1.2):
    # TR-even spin-orbit term tau_y sigma_x preserves two Kramers doublets.
    H=split*TAU[2]+.31*np.kron(PAULI[1],PAULI[0])
    V=np.stack([
        .55*np.kron(PAULI[0],PAULI[0])+.21*np.kron(PAULI[1],I),
        .47*np.kron(PAULI[0],PAULI[1])+.13*np.kron(PAULI[0],PAULI[2])])
    return H,V,SPIN

def random_unitary(n,seed):
    rng=np.random.default_rng(seed)
    a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    q,r=np.linalg.qr(a)
    return q

def cluster_rotation(H,seed,split=2):
    e,U=np.linalg.eigh(H);groups=[]
    labels=np.zeros(len(e),dtype=int)
    for i in range(1,len(e)):
        labels[i]=labels[i-1]+int(e[i]-e[i-1]>1e-7)
    Wb=np.zeros((len(e),len(e)),complex)
    rng=seed
    for g in range(labels[-1]+1):
        idx=np.flatnonzero(labels==g)
        Wb[np.ix_(idx,idx)]=random_unitary(len(idx),rng)
        rng+=1
    return U@Wb@U.conj().T

def rotate_data(data,U):
    H,V,S=data
    return (U@H@U.conj().T,
            np.stack([U@x@U.conj().T for x in V]),
            np.stack([U@x@U.conj().T for x in S]))

class NonAbelianClusterResponse(unittest.TestCase):
    def setUp(self):
        self.data=kramers_model()
        self.kw=dict(mu=0.,temperature=.17,z=.071+.009j,kappa=1.)

    def test_kramers_doublets_detected_as_projector_clusters(self):
        d=cluster_structure(self.data[0])
        self.assertEqual(d["dimensions"],[2,2])
        self.assertEqual(len(d["projectors"]),2)
        np.testing.assert_allclose(sum(d["projectors"]),np.eye(4),atol=1e-13)
        for P in d["projectors"]:
            np.testing.assert_allclose(P@P,P,atol=1e-13)
            self.assertAlmostEqual(float(np.trace(P).real),2.,places=12)

    def test_spectral_cluster_kubo_matches_independent_liouville(self):
        for z in (.071+.009j,-.043+.013j,.08):
            a=cluster_interband_kubo(*self.data,**{**self.kw,"z":z})
            b=cluster_interband_liouville(*self.data,**{**self.kw,"z":z})
            for x,y in zip(a,b):np.testing.assert_allclose(x,y,atol=2e-12,rtol=2e-12)

    def test_full_basis_unitary_covariance(self):
        reference=cluster_interband_kubo(*self.data,**self.kw)
        changed=rotate_data(self.data,random_unitary(4,21))
        got=cluster_interband_kubo(*changed,**self.kw)
        for x,y in zip(reference,got):np.testing.assert_allclose(x,y,atol=2e-12,rtol=2e-12)

    def test_arbitrary_rotations_inside_each_degenerate_doublet(self):
        W=cluster_rotation(self.data[0],31)
        np.testing.assert_allclose(W@self.data[0]@W.conj().T,self.data[0],atol=2e-13)
        reference=cluster_interband_kubo(*self.data,**self.kw)
        got=cluster_interband_kubo(*rotate_data(self.data,W),**self.kw)
        for x,y in zip(reference,got):np.testing.assert_allclose(x,y,atol=2e-12,rtol=2e-12)

    def test_nonabelian_qgt_trace_and_eigenvalues_are_gauge_invariant(self):
        ref=cluster_zeeman_qgt(*self.data,cluster_index=0)
        changed=cluster_zeeman_qgt(*rotate_data(self.data,random_unitary(4,51)),cluster_index=0)
        np.testing.assert_allclose(ref["trace"],changed["trace"],atol=2e-12,rtol=2e-12)
        for a in range(2):
            for b in range(3):
                x=np.sort_complex(np.linalg.eigvals(ref["qgt"][a,b]))
                y=np.sort_complex(np.linalg.eigvals(changed["qgt"][a,b]))
                np.testing.assert_allclose(x,y,atol=2e-12,rtol=2e-12)
        self.assertGreater(np.max(abs(ref["trace"])),1e-5)

    def test_qgt_projector_rotates_covariantly(self):
        U=random_unitary(4,71)
        ref=cluster_zeeman_qgt(*self.data,cluster_index=0)
        got=cluster_zeeman_qgt(*rotate_data(self.data,U),cluster_index=0)
        np.testing.assert_allclose(got["projector"],U@ref["projector"]@U.conj().T,atol=2e-12)

    def test_exact_time_reversal_at_kramers_point(self):
        H,V,S=self.data
        U_T=np.kron(I,1j*PAULI[1])
        transform=lambda A:U_T@A.conj()@U_T.conj().T
        np.testing.assert_allclose(transform(H),H,atol=1e-13)
        for v in V:np.testing.assert_allclose(transform(v),-v,atol=1e-13)
        for s in S:np.testing.assert_allclose(transform(s),-s,atol=1e-13)
        pos=cluster_interband_kubo(H,V,S,z=self.kw["z"],mu=0.,temperature=.17)
        neg=cluster_interband_kubo(H,V,S,z=-self.kw["z"],mu=0.,temperature=.17)
        np.testing.assert_allclose(pos[0],neg[0],atol=2e-12,rtol=2e-12)
        np.testing.assert_allclose(pos[1],-neg[1],atol=2e-12,rtol=2e-12)

    def test_near_degenerate_pair_is_flattened_consistently(self):
        H,V,S=self.data
        Hn=H+1e-11*np.diag([-1.,1.,-1.,1.])
        ref=cluster_interband_kubo(Hn,V,S,**self.kw,gap_tolerance=1e-8)
        W=cluster_rotation(Hn,81)
        got=cluster_interband_kubo(*rotate_data((Hn,V,S),W),**self.kw,gap_tolerance=1e-8)
        for x,y in zip(ref,got):np.testing.assert_allclose(x,y,atol=2e-12,rtol=2e-12)

    def test_invalid_cluster_inputs_fail_closed(self):
        with self.assertRaises(ValueError):cluster_zeeman_qgt(*self.data,cluster_index=4)
        with self.assertRaises(ValueError):cluster_interband_kubo(*self.data,gap_tolerance=-1)
        with self.assertRaises(ValueError):cluster_interband_kubo(self.data[0],self.data[1],self.data[2][:2])

if __name__=="__main__":unittest.main()
