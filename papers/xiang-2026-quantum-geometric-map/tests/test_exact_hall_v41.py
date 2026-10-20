"""Exact finite-torus linear Hall gate; NEVER a mixed EB response gate."""
import unittest
import numpy as np
from qgmap.exact_hall_v41 import torus,hall_lehmann,hall_retarded_crosscheck,twist_chern,audit

class ExactHallV41(unittest.TestCase):
    def test_clean_gap_chern_phases(self):
        for m,c in ((-3,0),(-1,-1),(1,1),(3,0)):
            r=twist_chern(L=3,mass=m,ntheta=5)
            self.assertAlmostEqual(r['chern'],c,delta=1e-10)
    def test_exact_twist_avg_equals_fhs(self):
        r=audit(L=3,mass=-1.,ntheta=5)
        self.assertLess(abs(r['twist_ave_minus_chern']),.02)
    def test_disorder_fixed_sample(self):
        r=audit(L=3,mass=-1.,W=.12,seed=1947,ntheta=5)
        self.assertLess(abs(r['twist_ave_minus_chern']),.02)
        self.assertGreater(r['minimum_fermi_gap_twist_grid'],0)
    def test_finite_frequency_static(self):
        h,vx,vy=torus(3,-1.)
        r=hall_lehmann(h,vx,vy)
        rr=hall_lehmann(h,vx,vy,eta=1e-4)
        self.assertAlmostEqual(r['hall_real'],rr['hall_real'],delta=1e-6)
        self.assertAlmostEqual(rr['hall_imag'],0.,delta=1e-11)
    def test_independent_retarded_kubo(self):
        for W in (0.,.12):
            h,x,y=torus(3,-1.,W,1947,(.24,.43))
            for freq,delta in ((.14,.03),(0.,.05),(.21,.015)):
                a=hall_lehmann(h,x,y,omega=freq,eta=delta)
                b=hall_retarded_crosscheck(h,x,y,omega=freq,eta=delta)
                self.assertAlmostEqual(a['hall_real'],b.real,delta=1e-10)
                self.assertAlmostEqual(a['hall_imag'],b.imag,delta=1e-10)
    def test_cross_contact_zero(self):
        L=3;W=.12;theta=(.17,.22);h=2e-4
        def x(a,b):return torus(L,-1,W,1947,(a,b))[0]
        contact=(x(theta[0]+h,theta[1]+h)-x(theta[0]+h,theta[1]-h)
                -x(theta[0]-h,theta[1]+h)+x(theta[0]-h,theta[1]-h))/(4*h*h)
        self.assertLess(np.max(abs(contact)),1e-8)
    def test_current_derivatives(self):
        L=3;w=.12;s=1234;theta=(.22,.41);h=1e-5
        H,Jx,Jy=torus(L,-1.,w,s,theta)
        for i,J in ((0,Jx),(1,Jy)):
            plus=list(theta);minus=list(theta);plus[i]+=h;minus[i]-=h
            Hp=torus(L,-1.,w,s,plus)[0];Hm=torus(L,-1.,w,s,minus)[0]
            self.assertLess(np.max(abs((Hp-Hm)/(2*h)-J/L)),1e-9)
    def test_fail_closed(self):
        with self.assertRaises(ValueError):torus(2)
        with self.assertRaises(ValueError):torus(3,W=-.1)
        with self.assertRaises(ValueError):torus(3,mass=0.)
        H,Jx,Jy=torus(3)
        with self.assertRaises(ValueError):hall_lehmann(H,Jx,Jy,mu=-100)
        with self.assertRaises(ValueError):twist_chern(3,ntheta=2)

if __name__=='__main__':unittest.main()
