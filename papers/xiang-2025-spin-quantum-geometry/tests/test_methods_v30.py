"""V30 multiband/PT: numerical identity, non-Abelian and frozen V16 checks."""
import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from methods.finite_eta import hermitian,equilibrium,resolvent,fd4,second_order
from methods.pt_doublet import (PT,STAGGERED,UNIFORM,I4,hamiltonian,
    projectors,projector_resolvent,transform_PT)

class MethodsV30Tests(unittest.TestCase):
    def test_multiband_independent_backends(self):
        for n in (2,3,4,5,6):
            rng=np.random.default_rng(20261020+n)
            A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));H=(A+A.conj().T)/2
            X=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
            for omega in (0.,.17):
                a=resolvent(H,X,.1,2,omega,'band')
                b=resolvent(H,X,.1,2,omega,'orbital')
                self.assertLess(np.max(abs(a-b)),3e-12)
    def test_liouville_equation(self):
        rng=np.random.default_rng(27)
        A=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));H=(A+A.conj().T)/2
        X=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        R=resolvent(H,X,.07,2,.12)
        lhs=(.14-.12j)*R+1j*(H@R-R@H)
        self.assertLess(np.max(abs(lhs-X)),1e-12)
    def test_equilibrium_covariance(self):
        H=np.diag([-1.,-.3,.5,2.]).astype(complex)
        rho=equilibrium(H,mu=.1,temperature=.23)
        Q,_=np.linalg.qr(np.random.default_rng(19).normal(size=(4,4))+1j*np.random.default_rng(20).normal(size=(4,4)))
        self.assertLess(np.max(abs(equilibrium(Q@H@Q.conj().T,.1,.23)-Q@rho@Q.conj().T)),2e-14)
    def test_PT_symmetry_and_rank2_projectors(self):
        self.assertLess(np.max(abs(PT@PT.conj()+I4)),1e-15)
        for k in ((0.,0.),(.38,.71),(-1.1,.9)):
            H,V=hamiltonian(*k)
            self.assertLess(np.max(abs(transform_PT(H)-H)),1e-14)
            self.assertLess(np.max(abs(transform_PT(V)-V)),1e-14)
            q=projectors(*k)
            for E,P in q:
                self.assertAlmostEqual(float(np.trace(P).real),2.,places=12)
                self.assertLess(np.max(abs(P@P-P)),1e-13)
            self.assertLess(np.max(abs(q[0][1]@q[1][1])),1e-12)
        self.assertLess(np.max(abs(transform_PT(STAGGERED)-STAGGERED)),1e-14)
        self.assertLess(np.max(abs(transform_PT(UNIFORM)+UNIFORM)),1e-14)
    def test_PT_projector_resolvent_vs_orbital(self):
        rng=np.random.default_rng(321)
        X=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        for k in ((.38,.71),(-1.1,.9)):
            H,_=hamiltonian(*k)
            for w in (0.,.17):
                a=projector_resolvent(*k,X,.08,1,w)
                b=resolvent(H,X,.08,1,w,'orbital')
                self.assertLess(np.max(abs(a-b)),1e-12)
    def test_U2_intra_pair_gauge(self):
        k=(.38,.71);H,_=hamiltonian(*k)
        rng=np.random.default_rng(16)
        X=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        e,U=np.linalg.eigh(H)
        W=np.zeros((4,4),complex)
        for start in (0,2):
            a=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
            W[start:start+2,start:start+2]=np.linalg.qr(a)[0]
        V=U@W
        b=V.conj().T@X@V
        result=V@(b/(.1+1j*(e[:,None]-e[None,:])))@V.conj().T
        oracle=projector_resolvent(*k,X,.1)
        self.assertLess(np.max(abs(result-oracle)),5e-13)
    def test_six_band_three_doublets(self):
        H4,_=hamiltonian(.38,.71)
        H=np.zeros((6,6),complex);H[:4,:4]=H4
        H[4,4]=3.15;H[5,5]=3.15
        rng=np.random.default_rng(666)
        X=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
        a=resolvent(H,X,.09,1,.13,'band')
        b=resolvent(H,X,.09,1,.13,'orbital')
        self.assertLess(np.max(abs(a-b)),2e-12)
    def test_PT_finite_eta_offline_V16_reference(self):
        d=second_order(lambda x,y:hamiltonian(x,y)[0],.38,.71,STAGGERED,
            mu=1.1,temperature=.17,eta=.1,step=.0004)
        self.assertAlmostEqual(d['value_real'],2.0232340477418247,delta=4e-8)
        self.assertLess(abs(d['value_imag']),1e-8)
        self.assertIn('NOT intrinsic',d['scientific_scope'])
    def test_PT_second_order_backends(self):
        f=lambda x,y:hamiltonian(x,y)[0]
        a=second_order(f,.38,.71,STAGGERED,mu=1.1,temperature=.17,eta=.1,step=.0004,backend='band')
        b=second_order(f,.38,.71,STAGGERED,mu=1.1,temperature=.17,eta=.1,step=.0004,backend='orbital')
        self.assertLess(abs(a['value_real']-b['value_real']),2e-7)
    def test_PT_omega_parity(self):
        f=lambda x,y:hamiltonian(x,y)[0]
        a=second_order(f,.38,.71,STAGGERED,mu=1.1,temperature=.17,eta=.1,omega=.02)
        b=second_order(f,.38,.71,STAGGERED,mu=1.1,temperature=.17,eta=.1,omega=-.02)
        self.assertLess(abs(a['value_real']-b['value_real']),1e-10)
    def test_regulator_is_not_a_tau0_response(self):
        f=lambda x,y:hamiltonian(x,y)[0]
        a=second_order(f,.38,.71,STAGGERED,mu=1.1,temperature=.17,eta=.1)
        b=second_order(f,.38,.71,STAGGERED,mu=1.1,temperature=.17,eta=.2)
        self.assertGreater(abs(a['value_real']-b['value_real']),1e-3)
        self.assertEqual(a['protocol'],'adiabatic_finite_eta_full_density_matrix_DC')
    def test_reject_bad_inputs(self):
        with self.assertRaises(ValueError):hermitian(np.ones((2,3)))
        with self.assertRaises(ValueError):equilibrium(np.eye(2),temperature=0.)
        with self.assertRaises(ValueError):resolvent(np.eye(2),np.eye(2),eta=0)
        with self.assertRaises(ValueError):resolvent(np.eye(2),np.eye(2),backend='unknown')
        with self.assertRaises(ValueError):second_order(lambda x,y:np.eye(2),0,0,np.eye(2),axis=3)

if __name__=='__main__':unittest.main()
