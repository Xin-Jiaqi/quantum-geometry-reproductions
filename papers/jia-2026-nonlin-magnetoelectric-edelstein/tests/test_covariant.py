"""Synthetic algebraic gauge tests; not DFT or full finite-T DC validation."""
import unittest
import numpy as np
from nmee_lab import (mesh, WannierTB2D, to_clean_model, evaluate_clean, CleanSpec,
                      CovariantBloch2D,evaluate_insulating_static,
                      covariant_from_wannier, covariant_derivative)
from nmee_lab.clean_core.core import MU_B
from nmee_lab.clean_core.models import SX,SY,SZ,I2


def toy():
    r=np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]])
    h0=np.array([[-.96,.08+.03j,.005,.01], [.08-.03j,-.45,.03j,.004],
                 [.005,-.03j,.40,.08],[.01,.004,.08,1.17]],complex)
    x=np.zeros((4,4),complex);x[0,2]=.027+.018j;x[1,3]=.016;x[2,0]=.006j
    y=np.zeros((4,4),complex);y[0,1]=.04;y[3,1]=.031j
    return WannierTB2D(r,np.array([h0,x,x.conj().T,y,y.conj().T]),
                np.array([[.4,0],[0,.38]]),'synthetic-gapped')


def physical_model(tb, nonzero_connection=True):
    n=tb.nwann
    s={'x':np.kron(I2,SX),'y':np.kron(I2,SY),'z':np.kron(I2,SZ)}
    ptop=np.diag([1,1,0,0])
    s['top_z']=.5*(ptop@s['z']+s['z']@ptop)
    rng=np.random.default_rng(42)
    a={}
    for b in (0,1):
        mat=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        a[b]=.012*(mat+mat.conj().T) if nonzero_connection else np.zeros((n,n),complex)
    def op(k,name):return np.broadcast_to(s[name],(len(k),n,n)).copy()
    def zero(k,b):return np.zeros((len(k),n,n),complex)
    m=CovariantBloch2D(hamiltonian=tb.matrices,derivative_h=lambda k,b:tb.matrices(k,b),
        connection=lambda k,b:np.broadcast_to(a[b],(len(k),n,n)).copy(),
        operators={x:(lambda k,x=x:op(k,x)) for x in s},
        derivative_operators={x:zero for x in s},reciprocal_nm=tb.reciprocal_nm)
    return m,s,a


def twist(m,X):
    """Exact periodic k-dependent orbital gauge U=e^{ik_x X}."""
    def W(k):
        phases=np.exp(1j*k[:,0,None]*np.diag(X)[None,:])
        return np.einsum('ki,ij->kij',phases,np.eye(len(X)))
    def tr(v,k):
        u=W(k);return u.conj().transpose(0,2,1)@v@u
    def comm(a,b):return a@b-b@a
    def h(k):return tr(m.hamiltonian(k),k)
    def dh(k,b):
        value=tr(m.derivative_h(k,b),k)
        return value+(1j*comm(h(k),X) if b==0 else 0)
    def a(k,b):
        c=tr(m.connection(k,b),k)
        return c-X if b==0 else c
    ops={x:(lambda k,x=x:tr(m.operators[x](k),k)) for x in m.operators}
    dop={x:(lambda k,b,x=x:(tr(m.derivative_operators[x](k,b),k)+
                     (1j*comm(ops[x](k),X) if b==0 else 0))) for x in m.operators}
    return CovariantBloch2D(h,dh,a,ops,dop,m.reciprocal_nm,m.g_factor,'periodic-gauge-twisted')


class TestCovariant(unittest.TestCase):
    def setUp(self):
        self.tb=toy();self.k=mesh(self.tb.reciprocal_nm,9,offset=(.31,.27))

    def test_zero_connection_matches_legacy_full_static(self):
        m,s,_=physical_model(self.tb,False)
        cov=evaluate_insulating_static(m,self.k,observable='z',magnetic='x',beta=0,n_occupied=2,B_step_T=.01)
        legacy=to_clean_model(self.tb,spin_onsite={key:s[key] for key in ('x','y','z')},acknowledge_zero_connection=True)
        old=evaluate_clean(legacy,self.k,CleanSpec(alpha='z',gamma='x',beta=0,mu_eV=0.,temperature_K=0))
        self.assertLess(abs(cov['response']-old['result']),3e-8)

    def test_finite_field_step_converges(self):
        m,_,_=physical_model(self.tb)
        vals=[evaluate_insulating_static(m,self.k,observable='top_z',magnetic='x',beta=1,n_occupied=2,B_step_T=step)['response'] for step in (.1,.05,.025)]
        self.assertLess(abs(vals[1]-vals[2]),abs(vals[0]-vals[1])*.27+1e-12)

    def test_covariant_velocity_and_operator_under_periodic_gauge(self):
        m,_,_=physical_model(self.tb)
        X=np.diag([0.,self.tb.lattice_nm[0,0],0.,self.tb.lattice_nm[0,0]])
        mt=twist(m,X)
        k=self.k
        W=np.einsum('ki,ij->kij',np.exp(1j*k[:,0,None]*np.diag(X)[None,:]),np.eye(4))
        ud=W.conj().transpose(0,2,1)
        for b in (0,1):
            v=covariant_derivative(m.hamiltonian(k),m.derivative_h(k,b),m.connection(k,b))
            vt=covariant_derivative(mt.hamiltonian(k),mt.derivative_h(k,b),mt.connection(k,b))
            self.assertLess(np.max(abs(vt-ud@v@W)),5e-14)
            ds=covariant_derivative(m.operators['x'](k),m.derivative_operators['x'](k,b),m.connection(k,b))
            dst=covariant_derivative(mt.operators['x'](k),mt.derivative_operators['x'](k,b),mt.connection(k,b))
            self.assertLess(np.max(abs(dst-ud@ds@W)),5e-14)

    def test_response_invariant_under_k_dependent_gauge(self):
        m,_,_=physical_model(self.tb)
        mt=twist(m,np.diag([0.,.4,0.,.4]))
        for comp,mag,beta in (('z','x',0),('top_z','x',1),('x','z',0)):
            a=evaluate_insulating_static(m,self.k,observable=comp,magnetic=mag,beta=beta,n_occupied=2)['response']
            b=evaluate_insulating_static(mt,self.k,observable=comp,magnetic=mag,beta=beta,n_occupied=2)['response']
            self.assertLess(abs(a-b),5e-10)

    def test_real_connection_changes_response(self):
        a,_,_=physical_model(self.tb,True)
        b,_,_=physical_model(self.tb,False)
        x=evaluate_insulating_static(a,self.k,observable='z',magnetic='x',beta=0,n_occupied=2)['response']
        y=evaluate_insulating_static(b,self.k,observable='z',magnetic='x',beta=0,n_occupied=2)['response']
        self.assertGreater(abs(x-y),1e-8)

    def test_layer_sum(self):
        m,s,a=physical_model(self.tb)
        sb=np.diag([0,0,1,1]);ob=.5*(sb@s['z']+s['z']@sb)
        m.operators['bottom_z']=lambda k:np.broadcast_to(ob,(len(k),4,4)).copy()
        m.derivative_operators['bottom_z']=lambda k,b:np.zeros((len(k),4,4),complex)
        out={s:evaluate_insulating_static(m,self.k,observable=s,magnetic='x',beta=0,n_occupied=2)['response'] for s in ('z','top_z','bottom_z')}
        self.assertLess(abs(out['top_z']+out['bottom_z']-out['z']),3e-10)

    def test_insulator_gap_gate(self):
        m,_,_=physical_model(self.tb)
        with self.assertRaises(ValueError):evaluate_insulating_static(m,self.k,observable='z',magnetic='x',beta=0,n_occupied=2,min_direct_gap_eV=100.)
        with self.assertRaises(ValueError):evaluate_insulating_static(m,self.k,observable='z',magnetic='x',beta=0,n_occupied=0)

    def test_no_implicit_position_connection(self):
        m,s,_=physical_model(self.tb)
        with self.assertRaises(ValueError):covariant_from_wannier(self.tb,position_connections={0:m.connection},
                            spin_operators=m.operators,spin_derivatives=m.derivative_operators)
        with self.assertRaises(ValueError):covariant_from_wannier(self.tb,position_connections={0:m.connection,1:m.connection},
                            spin_operators=m.operators,spin_derivatives={})
        wm=covariant_from_wannier(self.tb,position_connections={0:m.connection,1:m.connection},
                                 spin_operators=m.operators,spin_derivatives=m.derivative_operators)
        self.assertLess(abs(evaluate_insulating_static(wm,self.k,observable='z',magnetic='x',beta=0,n_occupied=2)['response']-
                            evaluate_insulating_static(m,self.k,observable='z',magnetic='x',beta=0,n_occupied=2)['response']),1e-12)

if __name__=='__main__':unittest.main()
