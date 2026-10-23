"""Round27 finite-q gauge identity and continuity checks; 2D finite torus only."""
import json
from pathlib import Path
import unittest
import numpy as np
from flux_lab.finite_q_gauge_ward import (
    soc_torus,site_phase,fermi_density,bond_currents,divergence,
    linear_spectral,linear_liouville,pure_gauge_response,
    continuity_residual,finite_amplitude_gauge)

Z=.42+.09j

def reference():
    t=soc_torus();L=site_phase(t)
    return t,L,pure_gauge_response(t,L,z=Z)

class FiniteQGaugeWard(unittest.TestCase):
    def test_finite_periodic_2d_rashba_hermitian_and_soc_torque(self):
        t,L,a=reference()
        self.assertEqual(t.H.shape,(18,18))
        np.testing.assert_allclose(t.H,t.H.conj().T,atol=1e-14)
        S=np.kron(np.eye(9),np.diag([1.,-1.]))
        self.assertGreater(np.max(abs(t.H@S-S@t.H)),.1)

    def test_microscopic_bond_continuity_operator_exact(self):
        t,_,_=reference()
        raw=1j*(t.H[None]@t.projectors-t.projectors@t.H)
        np.testing.assert_allclose(raw+divergence(t,bond_currents(t)),0,atol=1e-14)

    def test_pure_gauge_spectral_equals_exact_density_matrix_rotation(self):
        _,_,a=reference()
        np.testing.assert_allclose(a['dr_gauge'],a['dr_exact'],atol=2e-13)
        self.assertGreater(np.max(abs(a['dr_exact'])),.1)

    def test_pure_gauge_density_and_full_physical_bond_current_zero(self):
        _,_,a=reference()
        np.testing.assert_allclose(a['density_gauge'],0,atol=2e-13)
        np.testing.assert_allclose(a['current_gauge'],0,atol=2e-13)

    def test_scalar_vector_response_cancel_only_as_pure_gauge(self):
        _,_,a=reference()
        self.assertGreater(np.max(abs(a['density_vector'])),.05)
        self.assertGreater(np.max(abs(a['density_scalar'])),.05)
        np.testing.assert_allclose(a['density_vector']+a['density_scalar'],0,atol=2e-13)

    def test_missing_diamagnetic_bond_contact_fails_negative_control(self):
        _,_,a=reference()
        wrong=a['current_state_vector']+a['current_state_scalar']
        self.assertGreater(np.max(abs(wrong)),.3)
        np.testing.assert_allclose(wrong,-a['current_contact'],atol=1e-13)

    def test_vector_only_incomplete_pure_gauge_fails_density(self):
        _,_,a=reference()
        self.assertGreater(np.max(abs(a['density_vector'])),.05)

    def test_independent_orbital_liouville_matches_spectral(self):
        t,L,a=reference()
        full=linear_liouville(t.H,a['rho'],a['K']+a['Vscalar'],Z)
        np.testing.assert_allclose(full,a['dr_gauge'],atol=8e-13)
        scalar=linear_liouville(t.H,a['rho'],a['Vscalar'],Z)
        np.testing.assert_allclose(scalar,a['dr_scalar'],atol=8e-13)

    def test_each_drive_obeys_finite_q_continuity_with_correct_contact(self):
        t,_,a=reference()
        for drive,dr,contact in (
           (a['Vscalar'],a['dr_scalar'],None),
           (a['K'],a['dr_vector'],a['bond_J_contact']),
           (a['K']+a['Vscalar'],a['dr_gauge'],a['bond_J_contact'])):
            res=continuity_residual(t,dr,drive,a['rho'],Z,contact)
            np.testing.assert_allclose(res,0,atol=2e-13)

    def test_continuity_contact_omission_yields_local_violation(self):
        t,_,a=reference()
        wrong=continuity_residual(t,a['dr_vector'],a['K'],a['rho'],Z)
        self.assertGreater(np.max(abs(wrong)),.1)

    def test_exact_finite_peierls_phase_transforms_bond_current(self):
        t,L,a=reference()
        for eps in (.003,.17,.61):
            U,Hg=finite_amplitude_gauge(t,L,eps)
            np.testing.assert_allclose(Hg,U@t.H@U.conj().T,atol=1e-14)
            np.testing.assert_allclose(bond_currents(t,Hg),U[None]@a['bond_J']@U.conj().T[None],atol=1e-14)

    def test_finite_amplitude_central_difference_independently_recovers_vertices(self):
        t,L,a=reference();eps=2e-6
        _,Hp=finite_amplitude_gauge(t,L,eps)
        _,Hm=finite_amplitude_gauge(t,L,-eps)
        Kfd=(Hp-Hm)/(2*eps)
        Jfd=(bond_currents(t,Hp)-bond_currents(t,Hm))/(2*eps)
        np.testing.assert_allclose(Kfd,a['K'],atol=5e-11)
        np.testing.assert_allclose(Jfd,a['bond_J_contact'],atol=5e-11)

    def test_exact_time_dependent_pure_gauge_schrodinger_equation(self):
        t,L,a=reference(); amp=.23;omega=.59;t0=.41;dt=2e-5
        def state(tval):
            lam=amp*np.sin(omega*tval)
            U,Hg=finite_amplitude_gauge(t,L,lam)
            return U@a['rho']@U.conj().T,Hg-amp*omega*np.cos(omega*tval)*L
        rho,h=state(t0);rp,_=state(t0+dt);rm,_=state(t0-dt)
        residual=1j*(rp-rm)/(2*dt)-(h@rho-rho@h)
        np.testing.assert_allclose(residual,0,atol=5e-10)
        np.testing.assert_allclose(np.einsum('aij,ji->a',t.projectors,rho),
            np.einsum('aij,ji->a',t.projectors,a['rho']),atol=1e-13)

    def test_distinct_nonzero_q_and_frequency_protocols(self):
        t,_,_=reference()
        for qx,qy in ((1,1),(2,1),(1,2)):
            L=site_phase(t,qx=qx,qy=qy)
            for z in (.23+.07j,.74+.11j):
                x=pure_gauge_response(t,L,z=z)
                self.assertLess(np.max(abs(x['density_gauge'])),3e-13)
                self.assertLess(np.max(abs(x['current_gauge'])),3e-13)

    def test_uniform_onsite_and_bad_retarded_gate_fail_closed(self):
        t,L,a=reference()
        for qx,qy in ((0,0),(3,1),(1,3),(-1,1)):
            with self.assertRaises(ValueError):site_phase(t,qx=qx,qy=qy)
        with self.assertRaisesRegex(ValueError,'retarded'):linear_spectral(t.H,a['rho'],a['K'],-.2j)
        with self.assertRaisesRegex(ValueError,'retarded'):pure_gauge_response(t,L,z=.2)
        with self.assertRaises(ValueError):fermi_density(t.H,temperature=0)
        with self.assertRaises(ValueError):soc_torus(nx=2)

    def test_frozen_round27_audit_regression(self):
        t,L,a=reference()
        p=Path(__file__).resolve().parents[1]/'results/frozen/round27_finite_q_gauge_ward.json'
        frozen=json.loads(p.read_text())
        self.assertEqual(frozen['schema'],'cui-round27-finite-q-u1-gauge-ward-v1')
        self.assertEqual(frozen['paper_status'],'PARTIAL')
        fresh={
           'scalar_density_abs_max':float(np.max(abs(a['density_scalar']))),
           'vector_density_abs_max':float(np.max(abs(a['density_vector']))),
           'wrong_missing_contact_current_abs_max':float(np.max(abs(a['current_state_vector']+a['current_state_scalar']))),
           'correct_total_current_abs_max':float(np.max(abs(a['current_gauge']))),
           'correct_total_density_abs_max':float(np.max(abs(a['density_gauge']))),
        }
        for k,v in fresh.items():
            self.assertAlmostEqual(v,frozen[k],delta=3e-11)

if __name__=='__main__':unittest.main()
