"""Independent spin/orbital and projector controls for the public method subset."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from flux_lab.multiband import zeeman_qgt,covariant_h_derivative
from flux_lab.core import PAULI,model,qgt_eigen

def make4():
    orbital=np.array([[1.,0.],[0.,-1.]])
    h=np.kron(np.diag([-.7,1.1]),np.eye(2))+np.kron(np.eye(2),.24*PAULI[0]+.47*PAULI[2])
    dh=np.stack([np.kron(orbital,.15*PAULI[0]),np.kron(np.eye(2),.33*PAULI[1]),np.kron(orbital,.2*PAULI[2])])
    spin=np.stack([np.kron(np.eye(2),x) for x in PAULI])
    return h,dh,spin

class MultiBand(unittest.TestCase):
    def test_two_band_agreement_with_round1(self):
        d,dk=model(.31,-.47,.26)
        h=np.einsum('a,aij->ij',d,PAULI)
        v=np.einsum('xa,aij->xij',dk,PAULI)
        a=zeeman_qgt(h,v,PAULI)
        b=qgt_eigen(d,dk,-1)
        np.testing.assert_allclose(a['tensor'][0],b['T'],atol=1e-13)
        for key in ['gN','gA','omegaN','omegaA']:
            np.testing.assert_allclose(a[key][0],b[key],atol=1e-13)
    def test_projector_independently_matches_spectral(self):
        h,v,s=make4()
        a=zeeman_qgt(h,v,s)
        b=zeeman_qgt(h,v,s,method='projector')
        np.testing.assert_allclose(a['pairs'],b['pairs'],atol=2e-13)
    def test_same_basis_rotation(self):
        h,v,s=make4()
        U=np.diag(np.exp(1j*np.array([.1,.3,-.17,.23])))
        transform=lambda A: U.conj().T@A@U
        a=zeeman_qgt(h,v,s)
        b=zeeman_qgt(transform(h),np.stack([transform(o) for o in v]),np.stack([transform(o) for o in s]))
        np.testing.assert_allclose(a['tensor'],b['tensor'],atol=1e-13)
    def test_spectrum_and_band_pair_sum(self):
        h,v,s=make4(); z=zeeman_qgt(h,v,s)
        np.testing.assert_allclose(z['tensor'],z['pairs'].sum(axis=1))
        self.assertTrue(np.allclose(np.diagonal(z['pairs'],axis1=0,axis2=1),0))
    def test_degeneracy_fails_closed(self):
        h=np.diag([0.,0.,1.]); v=np.zeros((3,3,3));s=v.copy()
        with self.assertRaisesRegex(ValueError,'degenerate'):
            zeeman_qgt(h,v,s)
    def test_covariant_derivative_under_k_dependent_gauge(self):
        h=np.array([[.2,.5],[.5,-.1]],complex); dh=np.zeros((3,2,2),complex)
        theta=.6;k=.34
        W=lambda x: np.diag([np.exp(1j*theta*x),1.])
        H=lambda x: W(x).conj().T@h@W(x)
        diff=(H(k+1e-6)-H(k-1e-6))/(2e-6)
        A=np.zeros((3,2,2),complex);A[0]=np.diag([-theta,0.])
        d=np.zeros_like(A);d[0]=diff
        cov=covariant_h_derivative(H(k),d,A)
        np.testing.assert_allclose(cov,np.zeros_like(cov),atol=1e-9)
    def test_missing_physical_spin_forbidden(self):
        h,v,s=make4()
        with self.assertRaises(ValueError):zeeman_qgt(h,v,s[:2])
    def test_wrong_operator_gauge_detected_in_control(self):
        h,v,s=make4(); U=np.diag(np.exp(1j*np.array([0.,.5,1.1,-.2])))
        trans=lambda A:U.conj().T@A@U
        a=zeeman_qgt(h,v,s)
        wrong=zeeman_qgt(trans(h),np.stack([trans(o) for o in v]),s)
        self.assertGreater(np.max(np.abs(a['tensor']-wrong['tensor'])),1e-5)

if __name__=='__main__': unittest.main()
