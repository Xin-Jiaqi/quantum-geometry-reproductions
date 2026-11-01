"""Source v3 literal equations and physically invariant counterexamples."""
import unittest
import numpy as np
from flux_lab.core import PAULI, qgt_eigen, qgt_analytic
from flux_lab.causal_response import (source_four_sector_geometric, source_kme_expansion,
    interband_kubo, interband_liouville, static_subtracted_gyrotropic)
from flux_lab.antiunitary_response import examples, operator_residual, kubo_symmetry_residual
from test_causal_response import fixture


def appendix112(H,V,S,w):
    # Literal source p12 Eq112, before static subtraction; no response helper.
    e,U=np.linalg.eigh(H); out=np.zeros((2,3),complex)
    f=1/(1+np.exp(np.clip((e-2.3)/.025,-700,700)))
    for n in range(len(e)):
        for m in range(len(e)):
            if n==m:continue
            eps=e[n]-e[m]
            r=np.array([1j*np.vdot(U[:,n],v@U[:,m])/(-eps) for v in V])
            s=np.array([np.vdot(U[:,m],q@U[:,n]) for q in S])
            T=np.zeros((3,3),complex); T[:2]=np.outer(r,s)
            gn=(T+T.T).real/2; ga=(T+T.T).imag/2
            bxA=-(T-T.T).real; bxN=-(T-T.T).imag
            out-=f[n]*(1j*w/eps*(2*gn-bxA)+(1+w*w/eps**2)*(bxN-2*ga))[:2]
    return out

class SourceCore(unittest.TestCase):
    def test_main16_appendix112_same_after_static_subtraction(self):
        H,V,S=fixture()
        for w in (.017,-.017):
            a=source_four_sector_geometric(H,V,S,mu=2.3,temperature=.025,omega_source=w)
            np.testing.assert_allclose(a['sigma_source_omega'],appendix112(H,V,S,w)-appendix112(H,V,S,0),atol=2e-15)

    def test_kme_appendix126_matches_independent_dynamics_main19_linear_does_not(self):
        H,V,S=fixture(); errors=[]; wrong=[]
        opts=dict(mu=2.3,temperature=.025,kappa=.5)
        static=interband_kubo(H,V,S,z=0,**opts)[1]
        for w in (.04,.02,.01):
            a=source_kme_expansion(H,V,S,mu=2.3,temperature=.025,omega_source=w)
            exact=interband_liouville(H,V,S,z=-w,**opts)[1]-static
            np.testing.assert_allclose(a['static'],static,atol=1e-13)
            errors.append(float(np.max(abs(exact-a['appendix126_dynamic']))))
            wrong.append(float(np.max(abs(exact-a['main19_printed_dynamic']))))
        self.assertGreater(errors[0]/errors[1],7.8)
        self.assertGreater(errors[1]/errors[2],7.8)
        self.assertGreater(wrong[-1]/errors[-1],1000)
        self.assertLess(abs(wrong[1]/wrong[2]-2),.01)

    def test_invalid_mirror_prime_counterexample_rejected(self):
        spec=examples()['mxT']; k=np.array([.37,-.52])
        def H(k):
            x,y=k; return np.einsum('a,aij->ij',[x,x*y,.3+x*x+y*y],PAULI)
        mismatch=np.max(abs(H(spec.R@k)-spec.U@H(k).conj()@spec.U.conj().T))
        self.assertGreater(mismatch,.7)
        self.assertLess(max(operator_residual('mxT',k).values()),1e-13)
        q=kubo_symmetry_residual('mxT',k)
        self.assertLess(max(q['chi'],q['alpha']),1e-12)

    def test_tangential_C3_has_normal_curvature(self):
        theta=2*np.pi/3; Q=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
        U=np.diag(np.exp(-.5j*theta*np.array([1.,-1.])))
        dk=np.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,0.]])
        for k in (np.array([.37,-.52]),np.array([.7,.2])):
            d=lambda q:np.array([-q[1],q[0],.3])
            H=lambda q:np.einsum('a,aij->ij',d(q),PAULI)
            np.testing.assert_allclose(H(Q@k),U@H(k)@U.conj().T,atol=2e-15)
            a=qgt_eigen(d(k),dk)['omegaN'][2]
            b=qgt_analytic(d(k),dk)['omegaN'][2]
            np.testing.assert_allclose(a,b,atol=1e-14)
            self.assertGreater(abs(a),.5)

if __name__=='__main__':unittest.main()
