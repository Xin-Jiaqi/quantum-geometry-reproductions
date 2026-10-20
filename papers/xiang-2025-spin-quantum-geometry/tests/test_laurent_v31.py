"""V31 multiband Laurent and PT isolated clusters: independent numeric gates."""
import json,sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from methods.finite_eta import second_order
from methods.pt_doublet import hamiltonian,STAGGERED
from methods.pt_synthetic_v31 import PTThreeDoublets
from methods.laurent_v31 import (adiabatic_laurent,independent_population_delta,
    second_order_outer_rate,spectral_clusters,explicit_three_cluster_term,
    reconstruct_response)
F=PTThreeDoublets()
REF=json.loads((ROOT/'results/REFERENCE_V31.json').read_text())

def pt_even_observable():
    rng=np.random.default_rng(14)
    z=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
    h=(z+z.conj().T)/2
    return (h+F.J@h.conj()@F.J.conj().T)/2

class LaurentV31Tests(unittest.TestCase):
    def test_four_band_projector_ranks(self):
        c=spectral_clusters(hamiltonian(.38,.71)[0])
        self.assertEqual(tuple(round(np.trace(P).real) for P in c.projectors),(2,2))
    def test_pt_three_doublet_global(self):
        self.assertLess(np.max(abs(F.J@F.J.conj()+np.eye(6))),1e-14)
        for k in (-.81,.37,1.23):
            self.assertLess(F.check_PT(k),1e-12)
            c=spectral_clusters(F(k))
            self.assertEqual(tuple(round(np.trace(P).real) for P in c.projectors),(2,2,2))
            self.assertGreater(c.gap_min,.5)
    def test_A_B_C_hermitian(self):
        for f,k,mu,T in ((lambda x,y:hamiltonian(x,y)[0],(.38,.71),1.1,.17),
                         (F,(.37,0.),.6,.32)):
            r=adiabatic_laurent(f,*k,mu=mu,temperature=T)
            for key in ('A','B','C_adiabatic','diagonal','mixing','occupation'):
                self.assertLess(np.max(abs(r[key]-r[key].conj().T)),5e-10)
    def test_population_independent_analytic(self):
        for f,k,mu,T in ((lambda x,y:hamiltonian(x,y)[0],(.38,.71),1.1,.17),
                         (F,(.37,0.),.6,.32)):
            r=adiabatic_laurent(f,*k,mu=mu,temperature=T)
            oracle=independent_population_delta(f,*k,mu=mu,temperature=T)
            self.assertLess(np.max(abs(r['delta_printed_candidate']-oracle)),3e-10)
    def test_laurent_adiabatic_limit_fourband(self):
        f=lambda x,y:hamiltonian(x,y)[0]
        coeff=adiabatic_laurent(f,.38,.71,mu=1.1,temperature=.17)
        errors=[]
        for eta in (.08,.04,.02,.01):
            full=second_order(f,.38,.71,STAGGERED,mu=1.1,temperature=.17,
                              eta=eta,step=4e-4)['rho2']
            errors.append(float(np.max(abs(full-reconstruct_response(coeff,eta)))))
        self.assertLess(errors[-1],1e-4)
        for a,b in zip(errors[:-1],errors[1:]):
            self.assertTrue(.42 < b/a < .58)
    def test_laurent_adiabatic_limit_sixband(self):
        coeff=adiabatic_laurent(F,.37,0.,mu=.6,temperature=.32)
        errors=[]
        for eta in (.08,.04,.02,.01):
            full=second_order_outer_rate(F,.37,0.,eta,mu=.6,temperature=.32)
            errors.append(float(np.max(abs(full-reconstruct_response(coeff,eta)))))
        self.assertLess(errors[-1],3e-4)
        for a,b in zip(errors[:-1],errors[1:]):self.assertTrue(.42 < b/a < .59)
    def test_old_v19_local_adiabatic_constant_and_difference(self):
        r=adiabatic_laurent(F,.37,0.,mu=.6,temperature=.32)
        O=pt_even_observable()
        for key in ('C_adiabatic','C_half_population_candidate','delta_printed_candidate'):
            actual=float(np.trace(O@r[key]).real)
            self.assertAlmostEqual(actual,REF['pt_6band_local'][key],delta=2e-8)
        self.assertGreater(abs(REF['pt_6band_local']['delta_printed_candidate']),1e-4)
        self.assertLess(np.max(abs(r['C_adiabatic']-r['C_half_population_candidate']-
                                   r['delta_printed_candidate'])),1e-11)
    def test_real_third_cluster_nonzero(self):
        tri=explicit_three_cluster_term(F,.37,0.,mu=.6,temperature=.32)
        self.assertAlmostEqual(float(np.max(abs(tri))),REF['pt_6band_local']['triple_norm'],delta=3e-8)
        self.assertGreater(np.max(abs(tri)),.001)
        tri2=explicit_three_cluster_term(lambda x,y:hamiltonian(x,y)[0],.38,.71,mu=1.1,temperature=.17)
        self.assertLess(np.max(abs(tri2)),1e-15)
    def test_1d_bz_integrated_population_not_identically_zero(self):
        O=pt_even_observable();results={}
        for N in (24,48):
            ks=-np.pi+(np.arange(N)+.5)*2*np.pi/N
            vals=[]
            for k in ks:
                r=adiabatic_laurent(F,float(k),0.,mu=.6,temperature=.32)
                vals.append([np.trace(O@r[q]).real for q in
                    ('C_adiabatic','C_half_population_candidate','delta_printed_candidate')])
            results[N]=np.mean(vals,axis=0)
        self.assertLess(np.max(abs(results[24]-results[48])),4e-9)
        for value,target in zip(results[48],REF['pt_6band_bz48'].values()):
            self.assertAlmostEqual(value,target,delta=2e-8)
        self.assertGreater(abs(results[48][2]),1e-5)
    def test_outer_rate_diagnostic_changes_constant(self):
        f=lambda x,y:hamiltonian(x,y)[0]
        a=adiabatic_laurent(f,.38,.71,mu=1.1,temperature=.17,outer_rate=1.)
        b=adiabatic_laurent(f,.38,.71,mu=1.1,temperature=.17,outer_rate=2.)
        self.assertGreater(np.max(abs(a['C_adiabatic']-b['C_adiabatic'])),1e-4)
        self.assertIsNone(a['C_half_population_candidate'])
        for rate in (1.,3.):
            coeff=adiabatic_laurent(f,.38,.71,mu=1.1,temperature=.17,outer_rate=rate)
            errors=[]
            for eta in (.06,.03,.015):
                full=second_order_outer_rate(f,.38,.71,eta,rate,mu=1.1,temperature=.17)
                errors.append(float(np.max(abs(full-reconstruct_response(coeff,eta)))))
            self.assertLess(errors[1]/errors[0],.57)
            self.assertLess(errors[2]/errors[1],.57)
    def test_nondegenerate_three_four_five_bands(self):
        for n in (3,4,5):
            rng=np.random.default_rng(2700+n)
            B=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
            B=(B+B.conj().T)/2/np.linalg.norm(B)*.09
            en=np.linspace(-2.5,2.5,n)
            H=lambda k,y:np.diag(en)+np.sin(k)*B
            r=adiabatic_laurent(H,.31,0.,mu=.15,temperature=.42)
            self.assertEqual(r['cluster_ranks'],(1,)*n)
            self.assertLess(np.max(abs(r['delta_printed_candidate']-
                independent_population_delta(H,.31,0.,mu=.15,temperature=.42))),2e-9)
            v1=second_order_outer_rate(H,.31,0.,.04,mu=.15,temperature=.42)
            v2=second_order_outer_rate(H,.31,0.,.02,mu=.15,temperature=.42)
            a=np.max(abs(v1-reconstruct_response(r,.04)))
            b=np.max(abs(v2-reconstruct_response(r,.02)))
            self.assertLess(b/a,.62)
    def test_constant_orbital_gauge_covariance(self):
        rng=np.random.default_rng(238)
        A=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
        Q,_=np.linalg.qr(A)
        rotated=lambda x,y:Q@F(x,y)@Q.conj().T
        normal=adiabatic_laurent(F,.37,0.,mu=.6,temperature=.32)
        moved=adiabatic_laurent(rotated,.37,0.,mu=.6,temperature=.32)
        for key in ('A','B','C_adiabatic','delta_printed_candidate'):
            # Nested 4th-order differences amplify BLAS-dependent roundoff.
            self.assertLess(np.max(abs(Q@normal[key]@Q.conj().T-moved[key])),1e-8)
    def test_near_crossing_fails_closed(self):
        with self.assertRaises(ValueError):spectral_clusters(np.diag([-1.,-1.,-1.+1e-6,2.]))
        with self.assertRaises(ValueError):spectral_clusters(np.eye(2))
        with self.assertRaises(ValueError):spectral_clusters(np.diag([-1.,1.]),atol=1e-5,min_gap=1e-7)
    def test_not_mislabeled_intrinsic(self):
        r=adiabatic_laurent(F,.37,0.,mu=.6,temperature=.32)
        self.assertIn('NOT established physical intrinsic DC',r['scientific_scope'])
        self.assertEqual(r['protocol'],'adiabatic_L1_eta_L2_2eta_zero_frequency')

if __name__=='__main__':unittest.main()
