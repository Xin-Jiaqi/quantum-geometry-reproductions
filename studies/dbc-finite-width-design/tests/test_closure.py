import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy.integrate import quad
from optical_design.closure import (NormalForm,cusp_profile,skew_optimum,
    local_error_bound,optimizer_radius,candidate_winner,geometry_weight_exponent,audit_material_contract)
from optical_design.lineshapes import phi
from optical_design.projectors import projector_kernel
from optical_design.models import FAMILIES,h_only,kernels,matrices

class NormalFormTests(unittest.TestCase):
    def test_cross_term_shear_identity(self):
        n=NormalForm(2.,.08,-.2,.01,C=.4)
        x,y=.31,.28;a=.04;b=.003
        lhs=-n.A*x*x+a*y*y+n.C*x*y*y-n.B*y**4+b*y
        X=x-n.C*y*y/(2*n.A)
        rhs=-n.A*X*X+a*y*y-n.effective_B*y**4+b*y
        self.assertAlmostEqual(lhs,rhs,places=15)
    def test_cross_term_cannot_be_ignored(self):
        n=NormalForm(2.,.08,-.2,.01,C=.4)
        self.assertAlmostEqual(n.effective_B,.06)
        self.assertGreater(abs(n.B/n.effective_B-1),.3)
        with self.assertRaises(ValueError):NormalForm(1.,.1,-1.,1.,C=1.)
    def test_asymmetry_scaling(self):
        n=NormalForm(2.,.03,-.3,.01)
        v1=n.dimensionless(1e-4,-.01,1e-5)['nu']
        v2=n.dimensionless(1e-6,-.01,1e-5)['nu']
        self.assertAlmostEqual(v2/v1,100**.75)
    def test_even_limit_matches_old_shape(self):
        for shape in ('gaussian','lorentzian'):
            for xi,zeta in ((-.07,2.84),(-1.,0.),(.2,4.)):
                self.assertLess(abs(cusp_profile(xi,zeta,0,shape)[0]-phi(xi,zeta,shape,600)),2e-8)
    def test_tilt_sign_symmetry(self):
        for nu in (.1,.5,1.):
            self.assertAlmostEqual(cusp_profile(.1,2.8,nu)[0],cusp_profile(.1,2.8,-nu)[0],places=13)
    def test_exponents_include_kernel_weight(self):
        self.assertEqual(geometry_weight_exponent(),-.25)
        self.assertEqual(geometry_weight_exponent(q=2),.25)
        self.assertEqual(geometry_weight_exponent(p=2),.75)
    def test_invalid_normal_forms(self):
        for args in [(0.,1.,1.,1.),(1.,-1.,1.,1.),(1.,1.,0.,1.),(1.,1.,1.,0.)]:
            with self.assertRaises(ValueError):NormalForm(*args)
    def test_shape_envelope_fails_closed(self):
        with self.assertRaises(ValueError):cusp_profile(20.,1.)
        with self.assertRaises(ValueError):skew_optimum(3.)

class ErrorTests(unittest.TestCase):
    def test_kernel_energy_error_bound_direct_integral(self):
        for shape in ('gaussian','lorentzian'):
            width=.08;dK=.002;dE=.003
            L=(lambda t:np.exp(-.5*t*t)/np.sqrt(2*np.pi)) if shape=='gaussian' else (lambda t:1/(np.pi*(1+t*t)))
            f=lambda x: (.2+dK*x*x)*L((-x*x+dE*x-.04)/width)/width
            f0=lambda x:.2*L((-x*x-.04)/width)/width
            actual=abs(quad(lambda x:f(x)-f0(x),-1,1,epsabs=1e-12)[0])
            bound=local_error_bound(width,2,dK,dE,.2,shape)
            self.assertLess(actual,bound)
    def test_winner_requires_separation(self):
        self.assertTrue(candidate_winner([1.,1.3],.1)['certified_among_candidates'])
        self.assertFalse(candidate_winner([1.,1.1],.1)['certified_among_candidates'])
    def test_strong_concavity_parameter_radius(self):
        mu=2.;eps=.03;shift=eps/mu
        self.assertAlmostEqual(optimizer_radius(eps,mu),shift)
    def test_error_bounds_reject_nan(self):
        with self.assertRaises(ValueError):local_error_bound(np.nan,1,1,1,1)
        with self.assertRaises(ValueError):candidate_winner([1.,np.nan],.1)

class ProjectorTests(unittest.TestCase):
    def test_independent_projector_reference_all_paths(self):
        k=np.array([.31,2.23])
        for p in FAMILIES.values():
            for axes in ((0,0,1),(1,1,1),(0,0,0)):
                r=projector_kernel(lambda z:h_only(*z,p),k,axes=axes)
                expected=kernels(*k,p,axes)[1]
                np.testing.assert_allclose(r['kernels'],expected,atol=2e-9,rtol=3e-6)
    def test_constant_orbital_unitary(self):
        p=FAMILIES['mixed_orbital'];rng=np.random.default_rng(3)
        u=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
        a=projector_kernel(lambda z:h_only(*z,p),[.31,2.23])
        b=projector_kernel(lambda z:u.conj().T@h_only(*z,p)@u,[.31,2.23])
        np.testing.assert_allclose(a['kernels'],b['kernels'],atol=2e-9)
    def test_exact_doublet_trace_equals_two_copies(self):
        p=FAMILIES['constant'];k=[.31,2.23]
        a=projector_kernel(lambda z:h_only(*z,p),k)
        b=projector_kernel(lambda z:np.kron(h_only(*z,p),np.eye(2)),k,nocc=2,groups=[(0,1),(2,3),(4,5)])
        np.testing.assert_allclose(b['kernels'],2*a['kernels'],atol=4e-9,rtol=1e-5)
        np.testing.assert_allclose(b['gaps'],a['gaps'],atol=1e-12)
    def test_unequal_energy_cluster_is_rejected(self):
        p=FAMILIES['constant']
        with self.assertRaisesRegex(ValueError,'Split group'):
            projector_kernel(lambda z:h_only(*z,p),[.3,2.2],groups=[(0,),(1,2)])
    def test_abelian_exact_degeneracy_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'degeneracy'):
            projector_kernel(lambda z:np.diag([-1.,1.,1.]),[.3,.4])
    def test_projector_fd_step_convergence(self):
        p=FAMILIES['mixed_orbital'];k=[.21,2.1];target=kernels(*k,p)[1]
        a=projector_kernel(lambda z:h_only(*z,p),k,step=.003)['kernels']
        b=projector_kernel(lambda z:h_only(*z,p),k,step=.0015)['kernels']
        self.assertLess(np.max(abs(a-target)),1e-8)
        self.assertLess(np.max(abs(b-target)),1e-8)
    def test_projector_frame_rotation_in_degenerate_subspace(self):
        rng=np.random.default_rng(4)
        q=np.linalg.qr(rng.normal(size=(5,2))+1j*rng.normal(size=(5,2)))[0]
        u=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))[0]
        np.testing.assert_allclose(q@q.conj().T,(q@u)@(q@u).conj().T,atol=1e-14)

class ContractTests(unittest.TestCase):
    def fixture(self):
        return dict(energy_unit='eV',momentum_basis='cartesian',response='normal_shift_current',
            position_operator='validated_full',kernel_crosscheck=True,occupation='insulator',
            dimension=2,linewidth_role='spectral_only',tensor_frame='fixed_cartesian',bands_resolved=True,
            local_symmetry_checked=True,competing_features_checked=True,full_bz_validation=True,
            input_hash='synthetic',operator_gauge_id='same-basis',model_version='synthetic-v1',
            insulating_gap_eV=1.,width_eV=.01,line_shape='gaussian',width_convention='standard_deviation')
    def test_complete_contract(self):self.assertTrue(audit_material_contract(self.fixture())['accepted'])
    def test_missing_operator_blocks_material_claim(self):
        d=self.fixture();d['position_operator']='hr_and_centres'
        self.assertFalse(audit_material_contract(d)['accepted'])
    def test_width_mismatch_blocks(self):
        d=self.fixture();d['width_convention']='HWHM'
        self.assertFalse(audit_material_contract(d)['accepted'])
    def test_missing_occupations_blocks(self):
        d=self.fixture();d['occupation']='guessed'
        self.assertFalse(audit_material_contract(d)['accepted'])


class ExtendedScopeTests(unittest.TestCase):
    def test_quadratic_edge_both_line_shapes_and_sign(self):
        from optical_design.closure import quadratic_edge_response
        H=np.diag([2.,3.]);coef=1/(2*np.pi*np.sqrt(6))
        for shape in ('gaussian','lorentzian'):
            at=quadratic_edge_response(1.,1.,H,1.,.01,shape)
            self.assertAlmostEqual(at,coef/2)
            self.assertAlmostEqual(quadratic_edge_response(1.,1.,H,-1.,.01,shape),-at)
            self.assertAlmostEqual(quadratic_edge_response(1.,1.,-H,1.,.01,shape),at)
    def test_quadratic_edge_rejects_saddle(self):
        from optical_design.closure import quadratic_edge_response
        with self.assertRaises(ValueError):quadratic_edge_response(1,1,np.diag([-1,1]),1,.1)
    def test_3d_changes_critical_exponent(self):
        from optical_design.closure import phase_space_exponent
        self.assertEqual(phase_space_exponent([2,4]),-.25)
        self.assertEqual(phase_space_exponent([2,4,2]),.25)
    def test_signed_edges_cancel(self):
        from optical_design.closure import quadratic_edge_response
        for E in (.9,1.,1.1):
            plus=quadratic_edge_response(E,1.,np.eye(2),1.,.01)
            minus=quadratic_edge_response(E,1.,np.eye(2),-1.,.01)
            self.assertAlmostEqual(plus+minus,0.)


class MassiveDiracTests(unittest.TestCase):
    def test_independent_coarea_curve(self):
        from optical_design.closure import massive_dirac_integral
        for u in (.4,.8,.978):
            for m in (.01,.2):
                a=1-u*u
                for w in (1.3*m*np.sqrt(a),2*m):
                    D=np.sqrt(w*w-a*m*m)
                    def f(t):
                        x=D/np.sqrt(a)*np.cos(t);y=(u*w+D*np.sin(t))/a
                        R=np.sqrt(x*x+y*y+m*m)
                        speed=np.hypot(-D/np.sqrt(a)*np.sin(t),D/a*np.cos(t))
                        grad=np.hypot(x/R,y/R-u)
                        return y/R**2*speed/grad
                    direct=quad(f,0,2*np.pi,epsabs=1e-9,epsrel=1e-10)[0]
                    self.assertAlmostEqual(direct,massive_dirac_integral(w,u,m),places=7)
                    self.assertAlmostEqual(massive_dirac_integral(-w,u,m),-direct,places=7)
    def test_threshold_and_difference_sum_rule(self):
        from optical_design.closure import massive_dirac_integral as F
        for u in (.4,.95):
            m=.02;edge=m*np.sqrt(1-u*u)
            self.assertEqual(F(edge/2,u,m),0.)
            delta=quad(lambda w:F(w,u,m)-F(w,u,0),0,edge)[0]
            delta+=quad(lambda w:F(w,u,m)-F(w,u,0),edge,np.inf,epsabs=1e-9)[0]
            self.assertLess(abs(delta),1e-8)

if __name__=='__main__':unittest.main()
