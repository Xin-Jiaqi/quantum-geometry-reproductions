"""Independent orbital-superoperator / k-finite-difference checks of SM(6).

No production resolvent or Fermi Frechet derivative is used by the oracle.
"""
import unittest
from dataclasses import replace
import numpy as np
from scipy.special import expit
from nmee_lab import (DynamicSpec, evaluate_dynamic, cumnas_dynamic_model,
                      CovariantBloch2D)
from nmee_lab.clean_core.core import KB, MU_B
from nmee_lab.ptbackend import cumnas as cm
from nmee_lab.ptbackend.liouville_doublet import eb_kernel
from test_covariant import toy, physical_model, twist


def orbital_oracle(model, points, spec, step=2e-5):
    def equilibrium(k):
        e, u = np.linalg.eigh(model.hamiltonian(k[None])[0])
        return (u * expit((spec.mu_eV - e) / (KB * spec.temperature_K))) @ u.conj().T
    def comm(a, b):
        return a @ b - b @ a
    def solve(k, z, x):
        h = model.hamiltonian(k[None])[0]; n = len(h)
        # Row-major vectorization, built in orbital space.
        matrix = z * np.eye(n*n) - np.kron(h, np.eye(n)) + np.kron(np.eye(n), h.T)
        return np.linalg.solve(matrix, x.ravel()).reshape(n, n)
    def spin(k):
        return model.operators[spec.magnetic](k[None])[0]
    def magnetic_density(k):
        return solve(k, spec.magnetic_energy_eV, comm(spin(k), equilibrium(k)))
    def derivative(fun, k):
        dk = np.eye(2)[spec.beta] * step
        # Five-point differences, including any k-dependence of spin.
        diff = (-fun(k+2*dk)+8*fun(k+dk)-8*fun(k-dk)+fun(k-2*dk))/(12*step)
        A = model.connection(k[None], spec.beta)[0]
        return diff - 1j * comm(A, fun(k))
    values = []
    for k in points:
        rhoE = solve(k, spec.electric_energy_eV, 1j*derivative(equilibrium, k))
        source = 1j*derivative(magnetic_density, k) + comm(spin(k), rhoE)
        mixed = solve(k, spec.electric_energy_eV+spec.magnetic_energy_eV, source)
        values.append(np.trace(model.operators[spec.observable](k[None])[0] @ mixed))
    return -model.g_factor*MU_B*abs(np.linalg.det(model.reciprocal_nm))/(2*np.pi)**2*np.mean(values)


class DynamicTests(unittest.TestCase):
    def setUp(self):
        self.spec = DynamicSpec(.043+.017j, .071+.023j, .03, 100., 'A_z', 'x', 0)
        self.model = cumnas_dynamic_model(.15)
        self.k = np.array([[.37, -.52], [2.5, -1.1], [3.08, -.9]])

    def test_cumnas_exact_PT_matches_independent_projector(self):
        p = cm.Model(a_nm=.6)
        model = cumnas_dynamic_model(.15, p)
        k = self.k / p.a_nm
        expected, _ = eb_kernel(self.k, .15, .03, p,
            self.spec.electric_energy_eV, self.spec.magnetic_energy_eV, cm.PA@cm.IZ)
        actual = evaluate_dynamic(model, k, self.spec)['response']
        np.testing.assert_allclose(actual, expected.mean(), rtol=2e-10, atol=1e-12)

    def test_cumnas_full_orbital_equation_with_two_orderings(self):
        actual = evaluate_dynamic(self.model, self.k, self.spec)
        expected = orbital_oracle(self.model, self.k, self.spec)
        np.testing.assert_allclose(actual['response'], expected, rtol=2e-7, atol=1e-11)
        self.assertEqual(actual['response'], sum(actual['parts'].values()))
        for part in actual['parts'].values():
            self.assertGreater(abs(part), 1e-8)
            self.assertGreater(abs(part-expected), 1e-8)

    def test_nonzero_position_connection_and_k_dependent_gauge(self):
        base, _, _ = physical_model(toy())
        rotated = twist(base, np.diag([0., .4, 0., .4]))
        s = replace(self.spec, observable='top_z', magnetic='x', mu_eV=0.)
        a = evaluate_dynamic(base, self.k, s)['response']
        b = evaluate_dynamic(rotated, self.k, s)['response']
        np.testing.assert_allclose(a, b, atol=1e-12, rtol=1e-10)
        np.testing.assert_allclose(b, orbital_oracle(rotated, self.k, s), atol=1e-11, rtol=2e-7)
        # Missing derivative of the transformed spin must be numerically visible.
        wrong = replace(rotated, derivative_operators=base.derivative_operators)
        self.assertGreater(abs(evaluate_dynamic(wrong, self.k, s)['response']-a), 1e-8)

    def test_six_band_three_doublets_and_thermal_derivative(self):
        rng = np.random.default_rng(106)
        h0 = np.diag(np.repeat([-.06, .0, .07], 2))
        a = rng.normal(size=(6,6))+1j*rng.normal(size=(6,6)); v = .03*(a+a.conj().T)
        sx = np.kron(np.eye(3), np.array([[0.,1.],[1.,0.]]))
        sz = np.kron(np.eye(3), np.diag([1.,-1.]))
        zero = lambda k,b: np.zeros((len(k),6,6),complex)
        model = CovariantBloch2D(
            lambda k: h0+np.sin(k[:,0,None,None])*v,
            lambda k,b: np.cos(k[:,0,None,None])*v if b==0 else zero(k,b), zero,
            {'x':lambda k:np.broadcast_to(sx,(len(k),6,6)),
             'z':lambda k:np.broadcast_to(sz,(len(k),6,6))},
            {'x':zero,'z':zero}, 2*np.pi*np.eye(2))
        s = replace(self.spec, observable='z', magnetic='x', mu_eV=.004, temperature_K=80.)
        k = np.zeros((1,2))
        actual = evaluate_dynamic(model,k,s)['response']
        np.testing.assert_allclose(actual, orbital_oracle(model,k,s), rtol=2e-7, atol=1e-10)
        self.assertGreater(abs(actual),1e-7)

    def test_real_field_conjugacy_and_sublattice_additivity(self):
        a = evaluate_dynamic(self.model,self.k,self.spec)['response']
        negative = replace(self.spec,electric_energy_eV=-self.spec.electric_energy_eV.conjugate(),
                           magnetic_energy_eV=-self.spec.magnetic_energy_eV.conjugate())
        np.testing.assert_allclose(evaluate_dynamic(self.model,self.k,negative)['response'],a.conjugate(),atol=1e-12)
        b = evaluate_dynamic(self.model,self.k,replace(self.spec,observable='B_z'))['response']
        total = evaluate_dynamic(self.model,self.k,replace(self.spec,observable='z'))['response']
        np.testing.assert_allclose(a+b,total,atol=1e-12)

    def test_chunking_preserves_complex_response(self):
        a = evaluate_dynamic(self.model,self.k,self.spec)
        b = evaluate_dynamic(self.model,self.k,self.spec,chunk_size=1)
        np.testing.assert_allclose(a['response'],b['response'],atol=1e-13,rtol=1e-12)
        self.assertEqual(b['mesh_n'],len(self.k))
        with self.assertRaises(ValueError):
            evaluate_dynamic(self.model,self.k,self.spec,chunk_size=0)

    def test_gapped_diagonal_response_survives_small_drives(self):
        # H(k)=diag(-.5,.7)+sin(kx)*V; O=sigma_z is diagonal at k=0.
        # The printed two/three-band terms vanish there; SM(9)'s diagonal
        # contraction is nonzero. Check the full AC limit independently.
        h = np.diag([-.5,.7]); v = np.array([[.2,.3+.2j],[.3-.2j,-.1]])
        sx = np.array([[0,1],[1,0]]); sz = np.diag([1,-1])
        zero = lambda k,b: np.zeros((len(k),2,2),complex)
        model = CovariantBloch2D(lambda k:h+np.sin(k[:,0,None,None])*v,
            lambda k,b:np.cos(k[:,0,None,None])*v if b==0 else zero(k,b), zero,
            {'x':lambda k:np.broadcast_to(sx,(len(k),2,2)),
             'z':lambda k:np.broadcast_to(sz,(len(k),2,2))}, {'x':zero},2*np.pi*np.eye(2))
        # r_01=-i V_01/(e0-e1); +2 g mu_B f0 (O00-O11) Q01 / gap^2.
        Q01 = (-1j*v[0,1]/(-1.2)).real
        diagonal = 2*model.g_factor*MU_B*2*Q01/1.2**2
        errors=[]
        for eta in (1e-3,1e-4):
            s=DynamicSpec(1j*eta,1.31j*eta,0.,10.,'z','x',0)
            actual=evaluate_dynamic(model,np.zeros((1,2)),s)['response']
            expected=orbital_oracle(model,np.zeros((1,2)),s)
            np.testing.assert_allclose(actual,expected,rtol=1e-6,atol=1e-10)
            errors.append(abs(actual-diagonal))
        self.assertLess(errors[-1],errors[0]*.12)
        self.assertLess(errors[-1]/abs(diagonal),1e-4)

    def test_reject_undefined_drives_and_missing_physical_inputs(self):
        for change in ({'electric_energy_eV':.1}, {'magnetic_energy_eV':-.01j},
                       {'temperature_K':0.}, {'temperature_K':np.nan}, {'beta':2}):
            with self.assertRaises(ValueError):
                evaluate_dynamic(self.model,self.k,replace(self.spec,**change))
        with self.assertRaises(ValueError):
            evaluate_dynamic(replace(self.model,derivative_operators={}),self.k,self.spec)
        with self.assertRaises(ValueError):
            evaluate_dynamic(self.model,[[np.nan,0]],self.spec)


if __name__ == '__main__':
    unittest.main()
