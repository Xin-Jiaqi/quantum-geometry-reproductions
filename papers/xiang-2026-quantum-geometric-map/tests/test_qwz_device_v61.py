"""Independent v61 scientific gates: same-H bulk/edge comparison and negative controls."""
import unittest
import numpy as np
from qgmap.periodic_kubo_v59 import sx,sy,sz,qwz_h,periodic_chern
from qgmap.qwz_finite_marker_v60 import qwz_open
from qgmap.qwz_device_v61 import TX,TY,strip_principal_layer,qwz_cross_leads,qwz_four_terminal
from qgmap.matrix_leads_v60 import matrix_caroli,selfenergy_matrix
from qgmap.lead_transport_v57 import laplacian,four_terminal

class TestV61Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.top=qwz_four_terminal(5,-1,eta=1e-4)
        cls.inv=qwz_four_terminal(5,1,eta=1e-4)
        cls.trivial=qwz_four_terminal(5,3,eta=1e-4)

    def test_01_bloch_fourier_independently(self):
        for x,y in ((.2,.7),(-2.5,1.7),(.0,.0)):
            H=-1*sz+TX*np.exp(1j*x)+TX.conj().T*np.exp(-1j*x)+TY*np.exp(1j*y)+TY.conj().T*np.exp(-1j*y)
            self.assertTrue(np.allclose(H,qwz_h(x,y,-1),atol=1e-14))

    def test_02_exact_open_h(self):
        h,xy=qwz_open(5,4,-1)
        self.assertTrue(np.allclose(h,h.conj().T))
        self.assertTrue(np.allclose(h[:2,2:4],TX))
        self.assertTrue(np.allclose(h[:2,10:12],TY))
        self.assertEqual(h.shape,(40,40))

    def test_03_strips_exact_fourier(self):
        for transverse,forward in ((TY,TX),(TX,TY)):
            h,v=strip_principal_layer(4,-1,transverse,forward)
            self.assertTrue(np.allclose(h,h.conj().T))
            self.assertTrue(np.allclose(v[:2,:2],forward))
            self.assertTrue(np.allclose(h[:2,2:4],transverse))
            self.assertTrue(np.allclose(h[:2,:2],-sz))

    def test_04_lead_orientations(self):
        lead=qwz_cross_leads(5,4,-1,1e-4)
        for j,t in enumerate((TX.conj().T,TY,TX,TY.conj().T)):
            self.assertTrue(np.allclose(lead[j].h01[:2,:2],t))
            self.assertTrue(np.allclose(lead[j].coupling[np.nonzero(np.any(lead[j].coupling,axis=1))[0][0]:][:2,:2],t))

    def test_05_causality(self):
        for lead in qwz_cross_leads(5,5,-1,1e-4):
            _,_,gam,_=selfenergy_matrix(0,lead)
            self.assertGreaterEqual(np.linalg.eigvalsh(gam).min(),-1e-9)

    def test_06_transmission_nonnegative(self):
        self.assertGreaterEqual(self.top['transmissions'].min(),-1e-10)

    def test_07_caroli_probability_bound(self):
        self.assertLessEqual(self.top['transmissions'].max(),1+1e-6)

    def test_08_caroli_onsager_current_balance(self):
        T=self.top['transmissions']
        self.assertLess(np.max(abs(T.sum(axis=0)-T.sum(axis=1))),1e-10)

    def test_09_current_laplacian(self):
        T=self.top['transmissions'];o=four_terminal(T,np.zeros_like(T))
        v=o['voltages']
        self.assertAlmostEqual(v[0],.5,places=12)
        self.assertAlmostEqual(v[2],-.5,places=12)
        self.assertLess(abs((laplacian(T)@v)[1]),1e-12)
        self.assertLess(abs((laplacian(T)@v)[3]),1e-12)

    def test_10_floating_current(self):
        self.assertLess(self.top['floating_error'],1e-12)

    def test_11_total_current(self):
        self.assertLess(self.top['current_error'],1e-12)

    def test_12_bulk_quantization(self):
        self.assertAlmostEqual(periodic_chern(-1,31)['chern_wilson'],1,places=12)
        self.assertAlmostEqual(periodic_chern(1,31)['chern_wilson'],-1,places=12)

    def test_13_positive_chirality_hall_plateau(self):
        self.assertTrue(self.top['hall_resistance_reportable'])
        self.assertLess(abs(self.top['hall_resistance_h_over_q2']-1),3e-6)
        self.assertLess(abs(self.top['drive_current_q2_over_h']-1),.001)

    def test_14_negative_chirality_hall_plateau(self):
        self.assertLess(abs(self.inv['hall_resistance_h_over_q2']+1),3e-6)
        self.assertLess(abs(self.inv['drive_current_q2_over_h']-1),.001)

    def test_15_oriented_hall_sign_inversion(self):
        self.assertAlmostEqual(self.top['hall_ratio'],-self.inv['hall_ratio'],places=10)
        self.assertAlmostEqual(self.top['drive_current_q2_over_h'],self.inv['drive_current_q2_over_h'],places=10)

    def test_16_eta_convergence(self):
        coarse=qwz_four_terminal(5,-1,eta=1e-3)
        fine=qwz_four_terminal(5,-1,eta=1e-5)
        self.assertLess(abs(fine['hall_resistance_h_over_q2']-1),abs(coarse['hall_resistance_h_over_q2']-1))
        self.assertLess(abs(fine['drive_current_q2_over_h']-1),abs(coarse['drive_current_q2_over_h']-1))

    def test_17_size_convergence(self):
        errors=[abs(qwz_four_terminal(L,-1,eta=1e-4)['drive_current_q2_over_h']-1) for L in (3,5,7)]
        self.assertGreater(errors[0],errors[1]);self.assertGreater(errors[1],errors[2])

    def test_18_interior_chemical_potential(self):
        d=qwz_four_terminal(5,-1,mu=.25,eta=1e-4)
        self.assertTrue(d['hall_resistance_reportable'])
        self.assertLess(abs(d['hall_resistance_h_over_q2']-1),.01)

    def test_19_trivial_insulating_no_fictitious_hall_resistance(self):
        self.assertLess(abs(self.trivial['bulk_sigma_xy_q2_over_h']),1e-12)
        self.assertFalse(self.trivial['hall_resistance_reportable'])
        self.assertIsNone(self.trivial['hall_resistance_h_over_q2'])
        self.assertLess(self.trivial['drive_current_q2_over_h'],1e-7)

    def test_20_trivial_eta_artifact_suppression(self):
        a=qwz_four_terminal(5,3,eta=1e-3)
        b=qwz_four_terminal(5,3,eta=1e-4)
        self.assertGreater(a['drive_current_q2_over_h'],b['drive_current_q2_over_h']*90)
        self.assertIsNone(b['hall_resistance_h_over_q2'])

    def test_21_contact_perturbation_conditional(self):
        a=qwz_four_terminal(5,-1,eta=1e-4,contact_scale=.7)
        self.assertGreater(abs(a['drive_current_q2_over_h']-self.top['drive_current_q2_over_h']),1e-4)

    def test_22_local_orbital_phase_invariance(self):
        h,_=qwz_open(3,3,-1)
        leads=qwz_cross_leads(3,3,-1,1e-4)
        phase=np.exp(1j*np.linspace(-.7,.9,len(h)))
        U=np.diag(phase)
        from qgmap.matrix_leads_v60 import MatrixLead
        t0,_=matrix_caroli(0,h,np.zeros_like(h),leads)
        h1=U@h@U.conj().T
        leads1=tuple(MatrixLead(l.h00,l.h01,U@l.coupling,l.eta) for l in leads)
        t1,_=matrix_caroli(0,h1,np.zeros_like(h1),leads1)
        self.assertLess(np.max(abs(t0-t1)),1e-10)
        tbad,_=matrix_caroli(0,h1,np.zeros_like(h1),leads)
        self.assertGreater(np.max(abs(t0-tbad)),1e-5)

    def test_23_mixed_shape(self):
        d=qwz_four_terminal(4,-1,Lx=4,Ly=5,eta=1e-4)
        self.assertTrue(d['hall_resistance_reportable'])
        self.assertLess(abs(d['hall_resistance_h_over_q2']-1),.005)

    def test_24_missing_boundary_refused(self):
        with self.assertRaises(ValueError):qwz_cross_leads(2,5,-1)
        with self.assertRaises(ValueError):strip_principal_layer(1,-1,TY,TX)

    def test_25_bath_regulator_refused(self):
        with self.assertRaises(ValueError):qwz_cross_leads(5,5,-1,0)

    def test_26_threshold_explicit(self):
        with self.assertRaises(ValueError):qwz_four_terminal(5,-1,min_reportable_current=0)
        d=qwz_four_terminal(5,-1,eta=1e-4,min_reportable_current=2)
        self.assertIsNone(d['hall_resistance_h_over_q2'])

    def test_27_boundaries_distinct(self):
        leads=qwz_cross_leads(5,5,-1)
        for i in range(4):
            self.assertEqual(leads[i].coupling.shape,(50,10))
        self.assertEqual(int(np.count_nonzero(np.any(abs(leads[0].coupling)>0,axis=1) & np.any(abs(leads[1].coupling)>0,axis=1))),2)

    def test_28_independent_scattering_matrix_oracle(self):
        # Build a scattering matrix using square roots of independent lead Γ.
        # Checks S†S=I and |S_ab|² against the Caroli trace, not the Hall solver.
        h,_=qwz_open(3,3,-1)
        leads=qwz_cross_leads(3,3,-1,eta=1e-4)
        selfenergies=[selfenergy_matrix(0,l)[0] for l in leads]
        gammas=[1j*(s-s.conj().T) for s in selfenergies]
        gr=np.linalg.inv(-h-sum(selfenergies))
        W=[]
        for gamma in gammas:
            ev,U=np.linalg.eigh(gamma)
            W.append(U@np.diag(np.sqrt(np.maximum(ev,0))))
        weights=np.hstack(W)
        S=np.eye(weights.shape[1])-1j*weights.conj().T@gr@weights
        self.assertLess(np.linalg.norm(S.conj().T@S-np.eye(len(S))),1e-7)
        T,_=matrix_caroli(0,h,np.zeros_like(h),leads)
        n=h.shape[0]
        for a in range(4):
            for b in range(4):
                if a==b:continue
                block=S[a*n:(a+1)*n,b*n:(b+1)*n]
                self.assertAlmostEqual(np.linalg.norm(block)**2,T[a,b],places=9)

if __name__=='__main__':unittest.main()
