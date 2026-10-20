"""Original v18 validation: symmetric four-channel spectral Kubo sector only."""
import unittest
import numpy as np
from qgmap.spectral_kubo_v18 import spectral_bubble,slope
from qgmap.warped_disorder_v17 import Lattice,k_grid,h_and_vel
from qgmap.microscopic_ward_v11 import scba,ladder_vertex,I

class TestSpectralKuboV18(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=Lattice(n=81,alpha=.16,lam=.6)

    def test_no_disorder_ladder_equals_bare(self):
        p=Lattice(n=43,disorder_variance=0.)
        r=spectral_bubble(p,.03,vertex=True)
        s=spectral_bubble(p,.03,vertex=False)
        self.assertAlmostEqual(r['symmetric'],s['symmetric'],delta=1e-12)
        self.assertLess(r['vertex_residual'],1e-13)

    def test_four_channel_identity_bare(self):
        p=Lattice(n=47)
        h,vx,vy=h_and_vel(k_grid(p),p,.023)
        _,gr,ga,_=scba(h,p)
        v1=np.mean(np.trace(vx@(gr-ga)@vy@(gr-ga),axis1=-2,axis2=-1))*(-.25)
        direct=np.mean(np.trace(vx@((gr-ga)/(2j))@vy@((gr-ga)/(2j)),axis1=-2,axis2=-1))
        self.assertLess(abs(v1-direct),1e-12)

    def test_eigenbasis_spectral_bubble_no_disorder(self):
        p=Lattice(n=41,disorder_variance=0.)
        h,vx,vy=h_and_vel(k_grid(p),p,.025)
        energy,U=np.linalg.eigh(h)
        d=p.eta/((p.chemical_potential-energy)**2+p.eta**2)
        ud=np.swapaxes(U.conj(),-1,-2)
        xx=ud@vx@U;yy=ud@vy@U
        value=np.mean(np.einsum('...ij,...ji,...i,...j->...',xx,yy,d,d))
        reference=spectral_bubble(p,.025,vertex=False)['xy']
        self.assertLess(abs(value-reference),5e-13)

    def test_positive_spectral_density(self):
        p=Lattice(n=47)
        h,_,_=h_and_vel(k_grid(p),p,.01)
        _,gr,ga,_=scba(h,p)
        A=(ga-gr)/(2j)
        self.assertGreater(np.min(np.linalg.eigvalsh(A)), -1e-10)

    def test_symmetric_bubble_hermiticity(self):
        v=spectral_bubble(self.p,.02)
        self.assertLess(v['hermiticity_error'],2e-10)
        self.assertLess(v['imaginary_symmetric_part'],2e-10)
        self.assertLess(abs(v['xy']-v['yx']),2e-10)

    def test_scba_and_all_four_channel_residuals(self):
        v=spectral_bubble(self.p,.02)
        self.assertLess(v['scba_residual'],2e-10)
        self.assertLess(v['vertex_residual'],1e-10)
        self.assertEqual(set(v['channels']),{'ra','ar','rr','aa'})

    def test_hermitian_conjugate_channel(self):
        v=spectral_bubble(self.p,.02)
        self.assertAlmostEqual(v['channels']['ra'][0],v['channels']['ar'][1],delta=1e-10)
        self.assertAlmostEqual(v['channels']['rr'][0],v['channels']['aa'][1],delta=1e-10)

    def test_linear_dirac_zero_dense(self):
        p=Lattice(n=121,alpha=0,lam=0)
        self.assertLess(abs(slope(p,.02)['value']),1e-5)

    def test_warp_only_extra_reflection_zero_dense(self):
        p=Lattice(n=121,alpha=0,lam=.6)
        self.assertLess(abs(slope(p,.02)['value']),2e-5)

    def test_nonzero_wilson_branches(self):
        p=Lattice(n=121,alpha=.16,lam=.6)
        self.assertGreater(abs(slope(p,.02)['value']),.005)

    def test_mesh_convergence(self):
        p=Lattice(n=121,alpha=.16,lam=.6)
        q=Lattice(n=161,alpha=.16,lam=.6)
        self.assertLess(abs(slope(p,.02)['value']-slope(q,.02)['value']),1e-5)

    def test_ra_is_not_spectral(self):
        from qgmap.warped_disorder_v17 import slope as ra_slope
        p=Lattice(n=121,alpha=.16,lam=.6)
        self.assertGreater(abs(slope(p,.02)['value']-ra_slope(p,.02)['value']),.003)

    def test_shift_covariance(self):
        p=Lattice(n=161,alpha=.16,lam=.6)
        self.assertLess(abs(slope(p,.02,'shift')['value']),1e-6)

    def test_beta_zero_response(self):
        p=Lattice(n=161,beta=0.,alpha=.16,lam=.6)
        self.assertLess(abs(slope(p,.02)['value']),1e-6)

    def test_beta_reversal(self):
        p=Lattice(n=161,beta=.125,alpha=.16,lam=.6)
        q=Lattice(n=161,beta=-.125,alpha=.16,lam=.6)
        self.assertLess(abs(slope(p,.02)['value']+slope(q,.02)['value']),1e-6)

    def test_streda_sea_integrand_antisymmetric(self):
        p=Lattice(n=19)
        h,vx,vy=h_and_vel(k_grid(p),p,.02)
        _,gr,ga,_=scba(h,p)
        def integrand(v1,v2,G):
            dg=-G@G
            return np.mean(np.trace(v1@G@v2@dg-v1@dg@v2@G,axis1=-2,axis2=-1))
        for G in (gr,ga):
            self.assertLess(abs(integrand(vx,vy,G)+integrand(vy,vx,G)),1e-10)

    def test_unitary_covariance_of_all_channel_ladders(self):
        p=Lattice(n=23)
        h,vx,vy=h_and_vel(k_grid(p),p,.022)
        _,gr,ga,_=scba(h,p)
        theta=.38
        U=np.cos(theta)*I+1j*np.sin(theta)*np.array([[0,1],[1,0]],complex)
        V=U.conj().T
        for left,right in ((gr,ga),(ga,gr),(gr,gr),(ga,ga)):
            g1,C,_=ladder_vertex(left,right,vy,p.disorder_variance)
            g2,C2,_=ladder_vertex(V@left@U,V@right@U,V@vy@U,p.disorder_variance)
            self.assertLess(np.max(abs(C2-V@C@U)),1e-11)

    def test_field_step_behavior(self):
        p=Lattice(n=121,alpha=.16,lam=.6)
        a=slope(p,.02)['value'];b=slope(p,.01)['value'];c=slope(p,.005)['value']
        self.assertLess(abs(c-b),abs(b-a)*.35+1e-10)

if __name__=='__main__':unittest.main()
