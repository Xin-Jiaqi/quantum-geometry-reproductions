import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from optical_design.closure import cusp_profile

class PathTests(unittest.TestCase):
    def test_control_chain_rule_against_direct_perturbation(self):
        from optical_design.path import path_stationarity
        width=.08; xi=-.1; A=2.; B=.04; K=.01; a=.06; b=.003
        slopes=np.array([.2,.006,-.0004,-.03,.002])
        z=a/np.sqrt(B*width);nu=b/(B**.25*width**.75);h=1e-4
        f=lambda z,n:cusp_profile(xi,z,n)[0]
        phi=f(z,nu);fz=(f(z+h,nu)-f(z-h,nu))/(2*h);fn=(f(z,nu+h)-f(z,nu-h))/(2*h)
        out=path_stationarity(phi,fz,fn,width=width,A=A,B=B,kernel=K,a=a,b=b,
              A_prime=slopes[0],B_prime=slopes[1],kernel_prime=slopes[2],a_prime=slopes[3],b_prime=slopes[4])
        def direct(t):
            aa,bb,kk,curv,tilt=np.array([A,B,K,a,b])+t*slopes
            return kk*width**(-.25)/((2*np.pi)**2*np.sqrt(aa)*bb**.25)*f(curv/np.sqrt(bb*width),tilt/(bb**.25*width**.75))
        numeric=(direct(h)-direct(-h))/(2*h)
        self.assertLess(abs(numeric-out['control_derivative']),2e-10)
    def test_isolated_even_optimum_is_stationary(self):
        from optical_design.path import path_stationarity
        r=path_stationarity(3.,0.,0.,width=.001,A=2.,B=.04,kernel=.01,a=.01,b=0.,a_prime=-.3)
        self.assertEqual(r['control_derivative'],0.)
    def test_path_refuses_zero_kernel(self):
        from optical_design.path import path_stationarity
        with self.assertRaises(ValueError):path_stationarity(1,0,0,width=.1,A=1,B=1,kernel=0,a=0,b=0,a_prime=1)

if __name__=='__main__':unittest.main()
