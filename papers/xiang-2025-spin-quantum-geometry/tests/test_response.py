import unittest
from dataclasses import replace
import numpy as np
from methods.response import SpinModel,response_at_k,spin_observable,integrate_response
from methods.sm_eq2_v32 import Model,eq2_four_sources
from methods.finite_eta import second_order

SX=np.array([[0.,1.],[1.,0.]]);SY=np.array([[0,-1j],[1j,0]]);SZ=np.diag([1.,-1.])

def model():
    H=lambda x,y:.2*np.cos(x)*np.eye(2)+.7*np.sin(x)*SX+.5*np.sin(y)*SY+(.9+.1*np.cos(y))*SZ
    return SpinModel(H,lambda x,y,a:np.zeros((2,2)),lambda x,y,c:dict(x=SX,y=SY,z=SZ)[c],2*np.pi*np.eye(2))

class ResponseTests(unittest.TestCase):
    def test_literal_sm2_four_sources_multiband(self):
        for n in (2,4,6):
            src=Model(n);d=src.data(.23)
            def H(x,y):
                t=src.data(x);return (t['U']*t['en'])@t['U'].conj().T
            m=SpinModel(H,lambda x,y,a:np.zeros((n,n)),lambda x,y,c:np.eye(n),2*np.pi*np.eye(2))
            actual=response_at_k(m,[.23,0],energy1=.02j,energy2=.02j,axis2=0,mu=.25,thermal_energy=.45,step=.004)['density']
            expected=d['U']@eq2_four_sources(src,.23,.02).sum(axis=0)@d['U'].conj().T
            np.testing.assert_allclose(actual,expected,rtol=1e-6,atol=2e-7)

    def test_legacy_quadratic_same_protocol(self):
        m=model();k=[.4,.7];O=SZ/2
        old=second_order(m.hamiltonian,*k,O,mu=.1,temperature=.2,eta=.13,omega=.08)
        new=response_at_k(m,k,energy1=.08+.13j,energy2=-.08+.13j,axis2=0,mu=.1,thermal_energy=.2)
        np.testing.assert_allclose(new['density'],old['rho2'],rtol=1e-8,atol=1e-8)

    def test_two_frequency_backend_and_permutation(self):
        m=model();args=dict(energy1=.12+.1j,energy2=.07+.15j,axis1=0,axis2=1,current_axis=0,component='y')
        a=response_at_k(m,[.3,.5],**args)
        b=response_at_k(m,[.3,.5],**args,backend='orbital')
        c=response_at_k(m,[.3,.5],**dict(args,energy1=args['energy2'],energy2=args['energy1'],axis1=1,axis2=0))
        np.testing.assert_allclose(a['response'],b['response'],rtol=2e-7,atol=1e-9)
        np.testing.assert_allclose(a['density'],c['density'],atol=1e-12)

    def test_linear_drude_sign_and_quadratic_curved_band(self):
        # Scalar band: exact SM(1) rho1=f'(epsilon)*v/(eta-i omega).
        m=SpinModel(lambda x,y:np.array([[x*x]]),lambda x,y,a:np.zeros((1,1)),lambda x,y,c:np.ones((1,1)),np.eye(2))
        x=.3;T=.2;mu=.1;z=.1+.2j;f=1/(1+np.exp((x*x-mu)/T))
        first=-f*(1-f)/T*2*x/(-1j*z)
        actual=response_at_k(m,[x,0],energy1=z,mu=mu,thermal_energy=T)
        self.assertAlmostEqual(abs(actual['response']-first/2),0,places=10)
        self.assertAlmostEqual(spin_observable(m,[x,0],current_axis=0)[0,0].real,x,places=10)
        d2=f*(1-f)*(1-2*f)/T**2*(2*x)**2-f*(1-f)/T*2
        r=response_at_k(m,[x,0],energy1=z,energy2=z,axis2=0,mu=mu,thermal_energy=T)
        np.testing.assert_allclose(r['response'],d2/((-1j*z)*(-2j*z))/2,rtol=2e-7)

    def test_k_dependent_gauge_and_linear_conjugacy(self):
        m=model();X=np.diag([.4,0.])
        def W(x):return np.diag(np.exp(1j*x*np.diag(X)))
        tr=lambda a,x:W(x).conj().T@a@W(x)
        rotated=replace(m,hamiltonian=lambda x,y:tr(m.hamiltonian(x,y),x),
            connection=lambda x,y,a:-X if a==0 else np.zeros((2,2)),spin=lambda x,y,c:tr(m.spin(x,y,c),x))
        for extra in ({},{'energy2':.07+.15j,'axis2':1}):
            args=dict(energy1=.12+.1j,current_axis=0,component='y',**extra)
            a=response_at_k(m,[.3,.5],**args)['response'];b=response_at_k(rotated,[.3,.5],**args)['response']
            np.testing.assert_allclose(a,b,rtol=1e-6,atol=1e-8)
        a=response_at_k(m,[.3,.5],energy1=.12+.1j)['response']
        b=response_at_k(m,[.3,.5],energy1=-.12+.1j)['response']
        np.testing.assert_allclose(a.conjugate(),b,atol=1e-12)

    def test_reject_invalid_inputs(self):
        with self.assertRaises(ValueError):response_at_k(model(),[0,0],energy1=.1)
        with self.assertRaises(ValueError):response_at_k(model(),[0,0],energy1=.1j,thermal_energy=0)
        with self.assertRaises(ValueError):integrate_response(model(),[],energy1=.1j)

if __name__=='__main__':unittest.main()
