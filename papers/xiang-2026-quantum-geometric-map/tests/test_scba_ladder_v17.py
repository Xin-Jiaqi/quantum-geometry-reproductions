"""v17 independent periodic-lattice SCBA/ladder RA regression, NOT full physical PHE."""
import unittest
import numpy as np
from qgmap.warped_disorder_v17 import Lattice,k_grid,h_and_vel,ra_observable,slope,density_ward
from qgmap.microscopic_ward_v11 import I,SY

class V17(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=Lattice(n=101)
        cls.k=k_grid(cls.p)

    def test_parameters(self):
        for bad in ({'n':5},{'eta':0},{'disorder_variance':-1},{'beta':1.01}):
            with self.assertRaises(ValueError):Lattice(**bad)
        with self.assertRaises(ValueError):h_and_vel(self.k,self.p,protocol='bad')

    def test_mirror_y_H_field(self):
        k=self.k[::7,::9].copy(); mirrored=k.copy();mirrored[...,1]*=-1
        U=1j*SY
        for b in (0,.05):
            H,_,_=h_and_vel(k,self.p,b)
            HT,_,_=h_and_vel(mirrored,self.p,-b)
            self.assertLess(np.max(abs(H-U@HT@U.conj().T)),2e-14)

    def test_velocity_mirror_transform(self):
        k=self.k[::5,::11].copy(); mk=k.copy();mk[...,1]*=-1
        U=1j*SY
        _,vx,vy=h_and_vel(k,self.p,.02)
        _,vxm,vym=h_and_vel(mk,self.p,-.02)
        self.assertLess(np.max(abs(vx-U@vxm@U.conj().T)),2e-14)
        self.assertLess(np.max(abs(vy+U@vym@U.conj().T)),2e-14)

    def test_true_velocity_finite_k_difference(self):
        k=self.k[::8,::13];h=1e-5
        H,vx,vy=h_and_vel(k,self.p,.013)
        for a,v in [(0,vx),(1,vy)]:
            kp=k.copy();km=k.copy();kp[...,a]+=h;km[...,a]-=h
            d=(h_and_vel(kp,self.p,.013)[0]-h_and_vel(km,self.p,.013)[0])/(2*h)
            self.assertLess(np.max(abs(d-v)),1e-9)

    def test_onsite_B_derivative_is_SX(self):
        from qgmap.microscopic_ward_v11 import SX
        hp=h_and_vel(self.k,self.p,.02)[0]; hm=h_and_vel(self.k,self.p,-.02)[0]
        self.assertLess(np.max(abs((hp-hm)/.04+SX)),1e-13)

    def test_small_k_continuum_coefficients(self):
        from qgmap.warping_v16 import Surface,matrix_h
        k=np.array([[[.002,.003],[.006,-.004]]])
        H=h_and_vel(k,self.p,0)[0]
        reference=matrix_h(Surface(beta=self.p.beta,alpha=self.p.alpha,lam=self.p.lam),k[...,0],k[...,1])
        self.assertLess(np.max(abs(H-reference)),1e-6)

    def test_accidental_reflection_for_zero_alpha(self):
        p=Lattice(n=43,alpha=0.,lam=.6)
        k=k_grid(p);kk=k.copy();kk[...,1]=np.pi-kk[...,1]
        H,vx,vy=h_and_vel(k,p,.031)
        H2,vx2,vy2=h_and_vel(kk,p,.031)
        self.assertLess(np.max(abs(H-H2)),3e-14)
        self.assertLess(np.max(abs(vx-vx2)),3e-14)
        self.assertLess(np.max(abs(vy+vy2)),3e-14)

    def test_ladder_bethe_salpeter(self):
        r=ra_observable(self.p,.02)
        self.assertLess(r['vertex_residual'],1e-11)
        self.assertLess(r['scba_residual'],1e-10)
        self.assertTrue(np.isfinite(r['symmetric']))

    def test_scba_density_ward_finite_eta(self):
        r=density_ward(self.p,.02)
        self.assertLess(r['abs_error'],5e-10)
        self.assertLess(r['residual'],1e-11)
        self.assertLess(r['causality_max_eigenvalue'],0)

    def test_shift_control_full_BZ(self):
        p=Lattice(n=141,alpha=.16,lam=.6)
        self.assertLess(abs(slope(p,.02,'shift')['value']),1e-5)

    def test_nonzero_bare_and_ladder(self):
        p=Lattice(n=121,alpha=.16,lam=.6)
        bare=slope(p,.02,vertex=False)['value']
        dressed=slope(p,.02,vertex=True)['value']
        self.assertGreater(abs(bare),1e-3)
        self.assertGreater(abs(dressed),1e-3)
        self.assertGreater(abs(bare-dressed),1e-4)

    def test_zero_warp_only_dense_BZ(self):
        p=Lattice(n=161,alpha=0,lam=.6)
        self.assertLess(abs(slope(p,.02)['value']),2e-8)

    def test_beta_zero(self):
        p=Lattice(n=141,alpha=.16,lam=.6,beta=0)
        self.assertLess(abs(slope(p,.02)['value']),1e-5)

    def test_beta_reversal(self):
        p=Lattice(n=141,alpha=.16,lam=.6,beta=.125)
        q=Lattice(n=141,alpha=.16,lam=.6,beta=-.125)
        self.assertLess(abs(slope(p,.02)['value']+slope(q,.02)['value']),2e-5)

    def test_grid_convergence(self):
        coarse=slope(Lattice(n=121),.02)['value']
        dense=slope(Lattice(n=161),.02)['value']
        self.assertLess(abs(coarse-dense),2e-5)

if __name__=='__main__':unittest.main()
