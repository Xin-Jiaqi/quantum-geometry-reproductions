"""Independent analytic, finite-field, projector, symmetry and gap controls."""
import unittest
from unittest.mock import patch
import numpy as np
from numpy.testing import assert_allclose as close
from imhe_lab import *
from imhe_lab.core import jet, mul, mm, adj


class CoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.k, cls.w = mesh(24)
        cls.fields = layer_fields(cls.k)
        cls.g = geometry(*cls.fields.jets())
        cls.kfew = np.array([[.73, -.48], [1.21, .32], [-.51, .67]])
        cls.few = layer_fields(cls.kfew)

    def test_eigen_equation_and_normalization_all_slots(self):
        hj, vj = self.fields.jets()
        e, u, gap = eigenjets(hj)
        close(mm(hj, u), mul(u, e[:, :, None, :]), atol=2e-14)
        expected = jet(np.broadcast_to(np.eye(3), self.fields.h.shape))
        close(mm(adj(u), u), expected, atol=5e-15)
        self.assertGreater(gap, 1.)

    def test_base_spectral_curvature_and_velocity(self):
        e, o, v = exact_geometry(*self.fields.at(0, 0))
        close(self.g['energy'][0], e, atol=1e-14)
        close(self.g['omega'][0], o, atol=1e-14)
        close(self.g['velocity'][0], v, atol=1e-14)

    def test_full_band_curvature_and_energy_trace(self):
        close(self.g['omega'].sum(-1), 0, atol=1e-14)
        close(self.g['energy'][3].sum(-1), 0, atol=1e-14)

    def test_mixed_energy_exact_fields(self):
        g = geometry(*self.few.jets())
        numeric = mixed_difference(lambda a,b:exact_geometry(*self.few.at(a,b))[0], .002, .002)
        close(numeric, g['energy'][3], atol=3e-8, rtol=2e-5)

    def test_mixed_curvature_field_step_convergence(self):
        g = geometry(*self.few.jets())
        errors = []
        for h in (.004, .002):
            num = mixed_difference(lambda a,b:exact_geometry(*self.few.at(a,b))[1], h, h)
            errors.append(np.max(abs(num-g['omega'][3])))
        self.assertLess(errors[1], errors[0]/3)
        self.assertLess(errors[1], 4e-8)

    def test_projector_curvature_independent_k_derivatives(self):
        numeric = projector_curvature(lambda k:layer_fields(k).at(.03,.04)[0], self.kfew)
        expected = exact_geometry(*self.few.at(.03,.04))[1]
        close(numeric, expected, atol=3e-9)

    def test_moment_derivatives_from_orbital_matrix_finite_difference(self):
        for a in range(2):
            offset = np.eye(2)[a]*1e-4
            derivative = (layer_fields(self.kfew+offset).m-layer_fields(self.kfew-offset).m)/(2e-4)
            close(derivative, self.few.dm[a], atol=1e-9, rtol=1e-6)

    def test_moment_diagonal_A1_and_three_band_offdiagonals(self):
        e, u = np.linalg.eigh(self.few.h)
        mb = u.conj().swapaxes(-1,-2)@self.few.m@u
        v = np.stack([u.conj().swapaxes(-1,-2)@d@u for d in self.few.dh])
        delta = e[:,None,:]-e[:,:,None]
        delta = np.where(np.eye(3, dtype=bool), np.inf, delta)
        expected = (np.imag(v[0]*v[1].swapaxes(-1,-2))/delta).sum(-1)
        close(np.diagonal(mb,axis1=-2,axis2=-1), expected, atol=2e-14)
        self.assertGreater(np.max(abs(mb-np.eye(3)*expected[:,None,:])), 1e-3)

    def test_two_band_A5_has_no_offdiagonal_intermediate_band(self):
        sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.diag([1,-1])
        h=(.4*sx+.7*sy+1.2*sz)[None]; dh=np.stack((sx[None],sy[None]))
        m,_=orbital_moment(h,dh)
        e,u=np.linalg.eigh(h); mb=u.conj().swapaxes(-1,-2)@m@u
        close(mb[:,0,1],0,atol=1e-14)
        close(mb[:,0,0],mb[:,1,1],atol=1e-14)
        self.assertGreater(abs(mb[0,0,0]),.01)

    def test_intrinsic_A16_exact_conductivity_difference(self):
        r=transport(self.g,self.w,.2,.15)
        num=mixed_difference(lambda a,b:exact_transport(self.fields,self.w,.2,.15,a,b)[0],.002,.002)
        close(r['intrinsic_full'],num,atol=3e-10)
        self.assertGreater(abs(r['intrinsic_full']-r['intrinsic_eq6']),1e-4)
        close(sum(r['intrinsic_terms'].values()),r['intrinsic_full'],atol=1e-15)

    def test_full_drift_and_Eq7_are_distinct_with_independent_check(self):
        f=layer_fields(self.k,magnetic_order=.13);g=geometry(*f.jets());r=transport(g,self.w,.2,.15)
        num=mixed_difference(lambda a,b:exact_transport(f,self.w,.2,.15,a,b)[1],.002,.002)
        close(r['band_drift_full'],num,atol=5e-9,rtol=1e-4)
        self.assertGreater(np.max(abs(r['band_drift_full']-r['extrinsic_eq7'])),1e-4)
        close(r['extrinsic_eq7']+r['band_drift_cross'],r['band_drift_full'],atol=1e-15)

    def test_spinless_T_local_operator_and_geometry_parity(self):
        plus=geometry(*self.few.jets()); fminus=layer_fields(-self.kfew); minus=geometry(*fminus.jets())
        close(fminus.h,self.few.h.conj(),atol=1e-13)
        close(fminus.m,-self.few.m.conj(),atol=1e-13)
        close(minus['energy'][3],-plus['energy'][3],atol=1e-13)
        close(minus['omega'][0],-plus['omega'][0],atol=1e-13)
        close(minus['omega'][3],plus['omega'][3],atol=1e-13)

    def test_time_reversal_preserving_intrinsic_and_extrinsic(self):
        r=transport(self.g,self.w,.2,.15)
        self.assertGreater(abs(r['intrinsic_full']),1e-4)
        close(r['extrinsic_eq7'],0,atol=1e-15)
        close(r['band_drift_full'],0,atol=1e-15)

    def test_magnetic_domain_parity(self):
        a=transport(geometry(*layer_fields(self.k,magnetic_order=.13).jets()),self.w,.2,.15)
        b=transport(geometry(*layer_fields(self.k,magnetic_order=-.13).jets()),self.w,.2,.15)
        close(a['intrinsic_full'],b['intrinsic_full'],atol=2e-14)
        close(a['band_drift_full'],-b['band_drift_full'],atol=2e-14)
        close(a['extrinsic_eq7'],-b['extrinsic_eq7'],atol=2e-14)

    def test_response_exchange_symmetry_and_zero_Hall_work(self):
        r=transport(self.g,self.w,.2,.15); x=r['intrinsic_tensor']
        close(x,-x.T,atol=1e-15)
        bias=np.array([.31,-.7]);self.assertAlmostEqual(float(bias@x@bias),0)
        close(r['band_drift_full'],r['band_drift_full'].T,atol=1e-15)

    def test_zero_couplings(self):
        for attr,deriv in [('p','dp'),('m','dm')]:
            f=layer_fields(self.kfew);setattr(f,attr,np.zeros_like(getattr(f,attr)))
            setattr(f,deriv,np.zeros_like(getattr(f,deriv)))
            g=geometry(*f.jets());close(g['omega'][3],0,atol=1e-15);close(g['energy'][3],0,atol=1e-15)

    def test_constant_orbital_gauge(self):
        rng=np.random.default_rng(913); q,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        a=[q.conj().T@v@q for v in self.few.arrays()]
        rotated=Fields(*a)
        g1=geometry(*self.few.jets());g2=geometry(*rotated.jets())
        close(g1['omega'],g2['omega'],atol=1e-13)
        close(g1['energy'],g2['energy'],atol=1e-13)
        h,dh,ddh,_=layer_hamiltonian(self.kfew)
        mo,_=orbital_moment(q.conj().T@h@q,q.conj().T@dh@q)
        close(mo,q.conj().T@self.few.m@q,atol=1e-13)

    def test_arbitrary_band_phases(self):
        original=np.linalg.eigh
        def phased(h):
            e,u=original(h)
            angles=np.arange(e.size).reshape(e.shape)*.613
            return e,u*np.exp(1j*angles)[:,None,:]
        ref=geometry(*self.few.jets())
        with patch('numpy.linalg.eigh',phased):
            g=geometry(*self.few.jets())
        close(g['omega'],ref['omega'],atol=2e-14)
        close(g['energy'],ref['energy'],atol=2e-14)

    def test_six_band_contact_and_normalization(self):
        rng=np.random.default_rng(37)
        def herm(shape):
            a=rng.normal(size=shape)+1j*rng.normal(size=shape)
            return .03*(a+a.conj().swapaxes(-1,-2))
        h=np.diag(np.arange(6)*1.1)[None]+herm((2,6,6))
        f=Fields(h,herm((2,2,6,6)),herm((2,6,6)),herm((2,2,6,6)),
                 herm((2,6,6)),herm((2,2,6,6)),herm((2,6,6)),herm((2,2,6,6)))
        g=geometry(*f.jets());num=mixed_difference(lambda a,b:exact_geometry(*f.at(a,b))[1],.003,.003)
        close(g['omega'][3],num,atol=2e-9,rtol=1e-4)
        close(g['energy'][3].sum(-1),np.trace(f.c,axis1=-2,axis2=-1).real,atol=1e-13)

    def test_energy_origin_shift(self):
        r=transport(self.g,self.w,.2,.15)
        f=layer_fields(self.k);f.h+=7.3*np.eye(3)
        s=transport(geometry(*f.jets()),self.w,7.5,.15)
        close(r['intrinsic_full'],s['intrinsic_full'],atol=1e-13)

    def test_scalar_gate_is_chemical_potential_shift_not_layer_mixing(self):
        f=layer_fields(self.k);f.p[:]=np.eye(3)
        g=geometry(*f.jets());close(g['omega'][1],0,atol=1e-14)
        close(g['omega'][3],0,atol=1e-14)
        r=transport(g,self.w,.2,.15)
        num=mixed_difference(lambda a,b:exact_transport(f,self.w,.2,.15,a,b)[0],.0005,.0005)
        close(num,r['intrinsic_full'],atol=3e-9)

    def test_k_mesh_convergence(self):
        k,w=mesh(32);r=transport(geometry(*layer_fields(k).jets()),w,.2,.15)
        ref=transport(self.g,self.w,.2,.15)
        close(r['intrinsic_full'],ref['intrinsic_full'],atol=1e-12)

    def test_full_BZ_gapped_filled_band_control(self):
        r=transport(self.g,self.w,-.5,0)
        close(r['intrinsic_full'],0,atol=1e-12)
        close(r['band_drift_full'],0,atol=1e-15)
        close(r['intrinsic_terms']['energy_surface'],0,atol=1e-15)
        self.assertGreater(np.max(abs(self.g['omega'][3,:,0])),1e-3)

    def test_fermi_energy_derivatives(self):
        e=np.array([[.1,.3,.8]]);mu=.25;t=.19;h=1e-5
        f=fermi_derivatives(e,mu,t)
        for n in range(1,4):
            numeric=(fermi_derivatives(e+h,mu,t)[n-1]-fermi_derivatives(e-h,mu,t)[n-1])/(2*h)
            close(numeric,f[n],atol=2e-7)

    def test_printed_A4_inverse_gap_sign_counterexample(self):
        sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.])
        h=1.2*sz[None];dh=np.stack((sx[None],sy[None]));zeros=np.zeros_like(h)
        f=Fields(h,dh,sz[None],np.zeros_like(dh),zeros,np.zeros_like(dh))
        g=geometry(*f.jets());e=g['energy'];v=g['band_velocity']
        delta=e[0,0,0]-e[0,0,1]; dp=e[1,0,0]-e[1,0,1]
        printed=-1j*(v[1,0,0,0,1]*delta+v[0,0,0,0,1]*dp)/delta**2
        checked=g['r'][1,0,0,0,1]
        close(printed,-checked,atol=1e-14)
        self.assertGreater(abs(printed-checked),.1)

    def test_first_order_states_alone_miss_A8_second_order_velocity(self):
        sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.])
        h=1.2*sz[None];dh=np.stack((sz[None],sy[None]))
        f=Fields(h,dh,sz[None],np.zeros_like(dh),sx[None],np.zeros_like(dh))
        g=geometry(*f.jets());u=g['eigenvectors']
        truncated=(u[1].conj().swapaxes(-1,-2)@dh[0]@u[2]
                   +u[2].conj().swapaxes(-1,-2)@dh[0]@u[1])
        close(truncated,0,atol=1e-15)
        self.assertGreater(abs(g['band_velocity'][3,0,0,0,1]),.1)

    def test_recipe_recalculates_and_serializes(self):
        import json, os, subprocess, sys
        from pathlib import Path
        root=Path(__file__).resolve().parents[1]
        process=subprocess.run([sys.executable,str(root/'scripts/reproduce_core.py')],
                               capture_output=True,text=True,check=True,
                               env={**os.environ,'PYTHONPATH':str(root/'src')},timeout=60)
        result=json.loads(process.stdout)
        self.assertEqual(result['status'],'MODEL_READY')
        self.assertLess(result['validation']['field_steps'][-1]['intrinsic_abs_error'],1e-9)
        self.assertTrue(result['source_code_sha256'])

    def test_frozen_nonzero_benchmark(self):
        import json
        from pathlib import Path
        path=Path(__file__).resolve().parents[1]/'results/validation.json'
        frozen=json.loads(path.read_text())
        k,w=mesh(32);r=transport(geometry(*layer_fields(k).jets()),w,.2,.15)
        close(r['intrinsic_full'],frozen['time_reversal_symmetric']['intrinsic_full'],
              rtol=2e-10,atol=1e-13)

    def test_degenerate_and_malformed_inputs_rejected(self):
        f=layer_fields(self.kfew);f.h[:]=0
        with self.assertRaises(ValueError):geometry(*f.jets())
        f=layer_fields(self.kfew);f.p[0,0,1]+=.1
        with self.assertRaises(ValueError):f.jets()
        with self.assertRaises(ValueError):fermi_derivatives(self.g['energy'][0],.2,0)
        with self.assertRaises(ValueError):fermi_derivatives(self.g['energy'][0],.2,-1)
        with self.assertRaises(ValueError):transport(self.g,-self.w,.2,.15)
        with self.assertRaises(ValueError):mixed_difference(lambda a,b:a*b,0,.1)
        with self.assertRaises(ValueError):geometry(*self.fields.jets(),gap_tol=-1)


if __name__=='__main__':unittest.main()
