"""Independent algebra, symmetry, quadrature and scientific-scope checks."""
import unittest
from dataclasses import replace
from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from optical_design.models import *
from optical_design.lineshapes import optimum, profile_x,phi,phi_quad

class ModelChecks(unittest.TestCase):
    def test_hermitian_and_periodicity_all_paths(self):
        for p in FAMILIES.values():
            H,V,W=matrices(.27,2.18,p)
            for x in (H,V,W):np.testing.assert_allclose(x,np.swapaxes(x.conj(),-2,-1),atol=1e-14)
            np.testing.assert_allclose(H,h_only(.27+2*np.pi,2.18-2*np.pi,p),atol=4e-15)
    def test_spinless_time_reversal(self):
        for p in FAMILIES.values():
            np.testing.assert_allclose(h_only(-.37,-1.15,p),h_only(.37,1.15,p).conj(),atol=1e-14)
    def test_mass_reversal_is_mirror_related(self):
        M=np.diag([1,-1,1])
        for p in FAMILIES.values():
            np.testing.assert_allclose(h_only(-.31,1.7,replace(p,mass=-p.mass)),M@h_only(.31,1.7,p)@M,atol=1e-14)
    def test_analytic_first_derivatives(self):
        h=1e-5
        for p in FAMILIES.values():
            H,V,W=matrices(.31,2.3,p)
            for a in range(2):
                d=np.eye(2)[a]*h;k=np.array([.31,2.3])
                num=(h_only(*(k+d),p)-h_only(*(k-d),p))/(2*h)
                np.testing.assert_allclose(V[a],num,atol=1e-10,rtol=1e-9)
    def test_analytic_second_derivatives(self):
        h=1e-5;k=np.array([.31,2.3])
        for p in FAMILIES.values():
            H,V,W=matrices(*k,p)
            for a in range(2):
                d=np.eye(2)[a]*h
                num=(matrices(*(k+d),p)[1]-matrices(*(k-d),p)[1])/(2*h)
                np.testing.assert_allclose(W[a],num,atol=1e-10,rtol=1e-9)
    def test_eigen_derivatives(self):
        p=FAMILIES['mixed_orbital'];x,y=.29,2.8;h=2e-5
        f=(energy(x+h,y,p)-energy(x-h,y,p))/(2*h)
        np.testing.assert_allclose(f,deriv_x(x,y,p),atol=2e-8)
    def test_covariant_derivative_independent_finite_difference(self):
        k=np.array([.23,2.2]);h=1e-5
        for p in FAMILIES.values():
            E,U,V,W=bands(*k,p)
            for a,b in [(0,0),(0,1),(1,0),(1,1)]:
                vals=[]
                for s in [-1,1]:
                    d=np.eye(2)[a]*h*s;es,us=np.linalg.eigh(h_only(*(k+d),p))
                    us=us*np.exp(-1j*np.angle(np.sum(U.conj()*us,axis=0)))[None,:]
                    v=matrices(*(k+d),p)[1][b]
                    vals.append(-1j*(us[:,1].conj()@v@us[:,0])/(es[1]-es[0]))
                np.testing.assert_allclose((vals[1]-vals[0])/(2*h),derivative(E,V,W,1,0,a,b)[0],atol=2e-8)
    def test_u1_invariance_all_paths(self):
        phases=np.exp(1j*np.array([.38,-.87,1.13]));f=phases.conj()[:,None]*phases[None,:]
        for p in FAMILIES.values():
            E,U,V,W=bands(.32,2.3,p)
            a=kernels_from_bands(E,V,W);b=kernels_from_bands(E,V*f,W*f)
            for j in range(4):np.testing.assert_allclose(a[j],b[j],atol=1e-14)
    def test_response_tensor_mass_parity(self):
        for n in ['constant','bond_y','mixed_orbital']:
            p=FAMILIES[n];pm=replace(p,mass=-p.mass)
            for axes,sign in [((0,0,1),1),((1,0,1),-1)]:
                k=kernels(.35,2.3,p,axes)[1];km=kernels(-.35,2.3,pm,axes)[1]
                np.testing.assert_allclose(k,sign*km,atol=1e-12)
    def test_constant_matches_frozen_round3(self):
        c=critical(FAMILIES['constant'])
        self.assertAlmostEqual(c['mass'],.0993215269145,places=11)
        self.assertAlmostEqual(c['A'],2.31096686245,places=9)
        self.assertAlmostEqual(c['B'],.02293385,delta=1e-6)
        self.assertAlmostEqual(c['kernel'],.003389250067265,places=13)
    def test_hessian_independent_difference(self):
        h=1e-4;k=np.array([.01,np.pi-.03])
        for p in FAMILIES.values():
            g,H,e=grad_hess(*k,p)
            num=np.column_stack([(grad_hess(*(k+np.eye(2)[a]*h),p)[0]-grad_hess(*(k-np.eye(2)[a]*h),p)[0])/(2*h) for a in range(2)])
            np.testing.assert_allclose(H,num,atol=5e-7)
    def test_critical_roots_and_quartic_window(self):
        for p in FAMILIES.values():
            c=critical(p,.04);d=critical(p,.07)
            self.assertGreater(c['A'],0);self.assertGreater(c['B'],0)
            self.assertLess(abs(c['hessian_eigenvalues'][1]),1e-12)
            self.assertLess(abs(c['B']/d['B']-1),2e-5)
    def test_dark_mirror_and_zero_leading_kernel(self):
        p=FAMILIES['dark'];c=critical(p);p=replace(p,mass=c['mass'])
        np.testing.assert_allclose(h_only(.28,-2.2,p),h_only(.28,2.2,p),atol=1e-14)
        np.testing.assert_allclose(kernels(.28,-2.2,p)[1],-kernels(.28,2.2,p)[1],atol=1e-13)
        self.assertLess(abs(c['kernel']),1e-15)
    def test_insulating_sampled_domains(self):
        q=np.linspace(-np.pi,np.pi,71);x,y=np.meshgrid(q,q,indexing='ij')
        for p in FAMILIES.values():
            for m in [.01,.2]:
                E=np.linalg.eigvalsh(h_only(x,y,replace(p,mass=m)))
                self.assertGreater(float(np.min(E[...,1])-np.max(E[...,0])),1e-3)

class ShapeChecks(unittest.TestCase):
    def test_lorentzian_x_integral(self):
        for d in [-4.,0.,3.]:
            v=2*quad(lambda x:1/(np.pi*((x*x-d)**2+1)),0,np.inf,epsabs=1e-11)[0]
            self.assertAlmostEqual(v,float(profile_x(d,'lorentzian')),places=10)
    def test_gaussian_x_integral(self):
        for d in [-4.,0.,3.]:
            v=2*quad(lambda x:np.exp(-.5*(x*x-d)**2)/np.sqrt(2*np.pi),0,np.inf,epsabs=1e-11)[0]
            self.assertAlmostEqual(v,float(profile_x(d,'gaussian')),places=10)
    def test_shape_integrals_against_adaptive(self):
        for s in ['gaussian','lorentzian']:
            for xi,z in [(-1.,0),(-.1,2.8),(.5,4.)]:
                a=phi(xi,z,s,720);b=phi_quad(xi,z,s)[0]
                self.assertLess(abs(a-b),2e-8)
    def test_matched_fwhm_is_not_same_optimum(self):
        r={shape:optimum(shape) for shape in ('gaussian','lorentzian')}
        self.assertGreater(r['lorentzian']['separation_over_FWHM']/r['gaussian']['separation_over_FWHM'],1.8)
        self.assertLess(r['lorentzian']['separation_over_FWHM']/r['gaussian']['separation_over_FWHM'],2.0)
    def test_critical_scaling_change_of_variables(self):
        A,B,K=2.3,.023,.0034
        for g in [1e-4,1e-6]:
            jac=np.sqrt(g/A)*(g/B)**.25/g
            self.assertAlmostEqual(jac,g**(-.25)/(np.sqrt(A)*B**.25),places=12)
    def test_source_constant_optimum_recovered(self):
        r={shape:optimum(shape) for shape in ('gaussian','lorentzian')}['gaussian']
        self.assertAlmostEqual(r['zeta'],2.84126225829,places=7)
        self.assertAlmostEqual(r['separation_over_scale'],2.018192805,places=7)



if __name__=="__main__":unittest.main()
