import unittest
from dataclasses import replace
import numpy as np
from numpy.testing import assert_allclose
from orbital_igmc import CuMnAs, geometry, integrate, mesh, response_si, density_reference
from orbital_igmc.core import TZ, IZ

PROTOCOLS=('printed_D10','density_A9')

class ThreeBand:
    group_sizes=(1,1,1)
    def operators(self,q):
        x,y=np.asarray(q).T; n=len(q)
        h=np.zeros((n,3,3),complex);v=np.zeros((2,n,3,3),complex)
        h[:,0,0]=-1.3+.23*np.cos(x);h[:,1,1]=.3+.3*np.cos(y);h[:,2,2]=2.1+.2*np.cos(x+y)
        v[0,:,0,0]=-.23*np.sin(x);v[1,:,1,1]=-.3*np.sin(y);v[:,:,2,2]=-.2*np.sin(x+y)
        hops=[(0,1,.43*(np.sin(x)+1j*np.sin(y))+.13j*np.cos(x),.43*np.cos(x)-.13j*np.sin(x),.43j*np.cos(y)),
              (0,2,.38*(np.sin(y)+1j*np.sin(x)),.38j*np.cos(x),.38*np.cos(y)),
              (1,2,.49*(np.sin(x+y)+1j*np.cos(x)),.49*(np.cos(x+y)-1j*np.sin(x)),.49*np.cos(x+y))]
        for i,j,a,dx,dy in hops:
            h[:,i,j]=a;h[:,j,i]=a.conj()
            v[0,:,i,j]=dx;v[0,:,j,i]=dx.conj();v[1,:,i,j]=dy;v[1,:,j,i]=dy.conj()
        return h,v,np.diag([-.2,.1,.3])

class CoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=CuMnAs(a_angstrom=3.82);cls.q=mesh(41)
        cls.g=geometry(*cls.m.operators(cls.q),group_sizes=(2,2))
    def evaluate(self,m=None,n=41,mu=0,T=25,p='printed_D10'):
        m=m or self.m
        g=geometry(*m.operators(mesh(n)),group_sizes=m.group_sizes)
        return integrate(g,mu,T,protocol=p)

    def test_01_spectrum_eq38(self):
        x,y=self.q.T;m=self.m
        c=-m.t0*(np.cos(x)+np.cos(y))
        d=np.sqrt((2*m.t1*np.cos(x/2)*np.cos(y/2))**2+(m.exchange*np.cos(m.phi)-m.alpha*np.sin(y))**2+(m.exchange*np.sin(m.phi)+m.alpha*np.sin(x))**2)
        assert_allclose(self.g.energies,np.stack([c-d,c-d,c+d,c+d],-1),atol=3e-15)

    def test_02_derivatives_and_sewing(self):
        q=self.q[:11];h,v,_=self.m.operators(q);step=1e-5
        for a in (0,1):
            dq=np.zeros_like(q);dq[:,a]=step
            fd=(self.m.operators(q+dq)[0]-self.m.operators(q-dq)[0])/(2*step)
            assert_allclose(fd,v[a],atol=5e-11)
            dq[:,a]=2*np.pi
            assert_allclose(self.m.operators(q+dq)[0],TZ@h@TZ,atol=3e-15)

    def test_03_random_U2_gauge(self):
        h,v,m=self.m.operators(self.q[:37]);_,u=np.linalg.eigh(h)
        rng=np.random.default_rng(127)
        for s in (slice(0,2),slice(2,4)):
            z=rng.normal(size=(len(h),2,2))+1j*rng.normal(size=(len(h),2,2))
            w,_=np.linalg.qr(z);u[:,:,s]=u[:,:,s]@w
        a=geometry(h,v,m,group_sizes=(2,2));b=geometry(h,v,m,group_sizes=(2,2),eigenvectors=u)
        for p in PROTOCOLS:
            for key in ('orbital','spin'):
                assert_allclose(integrate(a,0,25,protocol=p)[key],integrate(b,0,25,protocol=p)[key],atol=2e-14)

    def test_04_fixed_orbital_basis_rotation(self):
        h,v,m=self.m.operators(self.q[:31]);rng=np.random.default_rng(19)
        w,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
        a=geometry(h,v,m,group_sizes=(2,2));b=geometry(w.conj().T@h@w,w.conj().T@v@w,w.conj().T@m@w,group_sizes=(2,2))
        for p in PROTOCOLS:
            for key in ('orbital','spin'):
                assert_allclose(integrate(a,0,25,protocol=p)[key],integrate(b,0,25,protocol=p)[key],atol=2e-14)

    def test_05_insulating_plateau(self):
        lo=self.g.energies[:,1].max();hi=self.g.energies[:,2].min()
        for p in PROTOCOLS:
            values=[integrate(self.g,lo+r*(hi-lo),0,protocol=p)['orbital'] for r in (.25,.5,.75)]
            assert_allclose(values,[values[0]]*3,atol=1e-15)
            self.assertGreater(abs(values[0][0]),.01)

    def test_06_neel_reversal(self):
        rev=replace(self.m,phi=self.m.phi+np.pi)
        for p in PROTOCOLS:
            for k in ('orbital','spin','total'):
                assert_allclose(self.evaluate(rev,p=p)[k],-self.evaluate(p=p)[k],atol=2e-14)

    def test_07_mirror_selection(self):
        for phi,zero_axis in ((0,1),(np.pi/2,0)):
            for p in PROTOCOLS:
                r=self.evaluate(replace(self.m,phi=phi),mu=.08,p=p)
                for k in ('orbital','spin'): self.assertLess(abs(r[k][zero_axis]),2e-14)

    def test_08_restore_time_reversal(self):
        for p in PROTOCOLS:
            r=self.evaluate(replace(self.m,exchange=0),n=40,p=p,T=300)
            assert_allclose(r['orbital'],0,atol=2e-14);assert_allclose(r['spin'],0,atol=2e-14)

    def test_09_remove_particle_hole_asymmetry(self):
        for p in PROTOCOLS:
            r=self.evaluate(replace(self.m,t0=0),p=p)
            assert_allclose(r['orbital'],0,atol=2e-14)
            self.assertGreater(abs(r['spin'][0]),.01)

    def test_10_density_projector_reference_insulator(self):
        ref=density_reference(self.m,self.q,0,25)
        fast=integrate(self.g,0,25,protocol='density_A9')
        for k in ('orbital','spin'):assert_allclose(fast[k],ref[k],rtol=3e-7,atol=2e-9)

    def test_11_density_projector_reference_metal(self):
        ref=density_reference(self.m,self.q,-.5,1000)
        fast=integrate(self.g,-.5,1000,protocol='density_A9')
        for k in ('orbital','spin'):assert_allclose(fast[k],ref[k],rtol=3e-7,atol=2e-9)

    def test_12_source_discrepancy_is_not_hidden(self):
        a=integrate(self.g,0,0,protocol='printed_D10');b=integrate(self.g,0,0,protocol='density_A9')
        assert_allclose(a['orbital'],2*b['orbital'],atol=2e-14)
        assert_allclose(a['spin'],b['spin'],atol=2e-14)
        assert_allclose(a['moment'],a['multiband'],atol=2e-14)
        self.assertGreater(np.linalg.norm(a['total']-b['total']),.01)

    def test_13_generic_three_band_reference(self):
        m=ThreeBand();q=mesh(32);g=geometry(*m.operators(q),group_sizes=m.group_sizes)
        a=integrate(g,0,500,protocol='density_A9');ref=density_reference(m,q,0,500)
        for k in ('orbital','spin'):assert_allclose(a[k],ref[k],rtol=2e-6,atol=2e-9)
        p=integrate(g,0,500,protocol='printed_D10')
        self.assertGreater(np.linalg.norm(p['multiband']),1e-5)
        # CuMnAs's factor two is NOT a universal correction for arbitrary bands.
        self.assertGreater(np.linalg.norm(a['orbital']-p['metric']-.5*(p['moment']+p['multiband'])),1e-4)

    def test_21_literal_D10_three_band_sum(self):
        model=ThreeBand();h,v,m=model.operators(np.array([[.37,-.61]]))
        e,u=np.linalg.eigh(h[0]);vv=np.array([u.conj().T@a[0]@u for a in v]);mm=u.conj().T@m@u
        r=np.zeros_like(vv)
        for n in range(3):
            for p in range(3):
                if n!=p:r[:,n,p]=1j*vv[:,n,p]/(e[p]-e[n])
        def triple(m,n,p):
            return -.5*((vv[0,m,n]+(m==n)*vv[0,p,p])*r[1,n,p]-(vv[1,m,n]+(m==n)*vv[1,p,p])*r[0,n,p])
        mom=np.zeros((2,3));multi=np.zeros((2,3));spin=np.zeros((2,3))
        for a in (0,1):
            for n in range(3):
                for m in range(3):
                    if m!=n:
                        ml=sum(triple(m,p,n) for p in range(3) if p!=n)
                        mom[a,n]-=(ml*r[a,n,m]/(e[m]-e[n])).real
                        spin[a,n]+=(mm[m,n]*r[a,n,m]/(e[m]-e[n])).real
                    for p in range(3):
                        if m!=p and p!=n:multi[a,n]+=(triple(m,n,p)*r[a,p,m]/(e[m]-e[p])).real
        g=geometry(h,v,model.operators(np.array([[.37,-.61]]))[2],group_sizes=(1,1,1))
        assert_allclose(mom,g.orbital_moment_term[:,0],atol=1e-15)
        assert_allclose(multi,g.orbital_multiband_term[:,0],atol=1e-15)
        assert_allclose(spin,g.spin_term[:,0],atol=1e-15)

    def test_14_zero_and_linear_frequency(self):
        r=integrate(self.g,0,25,protocol='printed_D10')
        a=response_si(r,self.m,.0001,convention='harmonic')
        b=response_si(r,self.m,.0002,convention='harmonic')
        zero=response_si(r,self.m,0,convention='harmonic')
        neg=response_si(r,self.m,-.0001,convention='harmonic')
        for k in a:
            assert_allclose(b[k],2*a[k],atol=1e-18);assert_allclose(zero[k],0,atol=0)
            assert_allclose(neg[k],a[k].conj(),atol=1e-18)

    def test_15_harmonic_phasor_conversion(self):
        r=integrate(self.g,0,25,protocol='density_A9')
        for k,v in response_si(r,self.m,.0001,convention='harmonic').items():
            assert_allclose(response_si(r,self.m,.0001,convention='phasor')[k],2*v,atol=0)

    def test_16_length_scaling(self):
        m2=replace(self.m,a_angstrom=2*self.m.a_angstrom)
        for p in PROTOCOLS:
            a=response_si(self.evaluate(p=p),self.m,.0001,convention='harmonic')
            b=response_si(self.evaluate(m2,p=p),m2,.0001,convention='harmonic')
            assert_allclose(b['orbital'],2*a['orbital'],atol=2e-17)
            assert_allclose(b['spin'],.5*a['spin'],atol=2e-17)

    def test_17_convergence_in_gap(self):
        a=self.evaluate(n=101);b=self.evaluate(n=161);c=self.evaluate(n=241)
        for k in ('orbital','spin'):
            assert_allclose(a[k],b[k],atol=2e-6)
            assert_allclose(b[k],c[k],atol=2e-9)

    def test_18_zero_T_metal_rejected(self):
        with self.assertRaisesRegex(ValueError,'Zero-T metal'):
            integrate(self.g,-.5,0,protocol='density_A9')

    def test_19_empty_and_full_occupation(self):
        for p in PROTOCOLS:
            for mu in (-10,10):
                r=integrate(self.g,mu,0,protocol=p)
                assert_allclose(r['orbital'],0,atol=2e-15);assert_allclose(r['spin'],0,atol=2e-15)

    def test_20_fail_closed(self):
        h,v,m=self.m.operators(self.q[:4])
        with self.assertRaises(ValueError):geometry(h,v,m,group_sizes=(1,1,1,1))
        with self.assertRaises(ValueError):geometry(h,v,m,group_sizes=(4,))
        with self.assertRaises(ValueError):geometry(h,v,m,group_sizes=(2,2),eigenvectors=np.zeros_like(h))
        vbad=v.copy();vbad[0,0,0,0]+=1
        with self.assertRaisesRegex(ValueError,'splits locally'):geometry(h,vbad,m,group_sizes=(2,2))
        hbad=h.copy();hbad[0,0,1]+=1j
        with self.assertRaises(ValueError):geometry(hbad,v,m,group_sizes=(2,2))
        with self.assertRaises(ValueError):mesh(2)
        with self.assertRaises(ValueError):CuMnAs(a_angstrom=-1)
        with self.assertRaises(ValueError):integrate(self.g,0,25,protocol='auto')
        with self.assertRaises(ValueError):integrate(self.g,0,-1,protocol='density_A9')
        with self.assertRaises(ValueError):response_si(self.evaluate(),self.m,1,convention='harmonic')
        with self.assertRaises(ValueError):response_si(self.evaluate(),self.m,.0001,convention='auto')

if __name__=='__main__':unittest.main()
