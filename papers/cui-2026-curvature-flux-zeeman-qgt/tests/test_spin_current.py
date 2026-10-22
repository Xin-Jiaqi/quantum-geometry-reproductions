"""Cui et al. Appendix C restricted physical-spin-current acceptance tests."""
import json
import unittest
from pathlib import Path
import numpy as np
from flux_lab.core import PAULI
from flux_lab.spin_current import (spin_orbital_benchmark, spin_current_operator,
     spin_current_derivative,intrinsic_spin_geometry,fermi,integrate_spin_rta_bz)

HERE=Path(__file__).resolve().parents[1]

def numeric_band_derivative(kx,ky,a,b,periodic,step=1e-5):
    result=[]
    for shift in (+step,-step):
        k=[kx,ky];k[b]+=shift
        h,v,_,s=spin_orbital_benchmark(*k,periodic=periodic)
        E,U=np.linalg.eigh(h)
        O=spin_current_operator(v,s)[a]
        result.append([np.vdot(U[:,n],O@U[:,n]).real for n in range(4)])
    return (np.array(result[0])-np.array(result[1]))/(2*step)

def finite_eigenstate_oracle(h,v,s,a,b,step=1e-5):
    """H+E*r_b finite eigenstate derivative with occupation kept fixed."""
    E,U=np.linalg.eigh(h);V=U.conj().T@v[b]@U
    R=np.zeros((4,4),complex)
    for n in range(4):
        for m in range(4):
            if n!=m:R[n,m]=1j*V[n,m]/(E[m]-E[n])
    assert np.max(abs(R-R.conj().T))<1e-12
    O=spin_current_operator(v,s)[a];r=U@R@U.conj().T
    results=[]
    for eps in (+step,-step):
        _,uu=np.linalg.eigh(h+eps*r)
        results.append(np.vdot(uu[:,0],O@uu[:,0]).real)
    return (results[0]-results[1])/(2*step)

class SpinCurrentAppendixC(unittest.TestCase):
    def test_operator_is_conventional_anticommutator(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        np.testing.assert_allclose(spin_current_operator(v,s),np.array([(x@s+s@x)/4 for x in v]))

    def test_intrinsic_complex_qgt_sector_exact_identity(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        d=intrinsic_spin_geometry(h,v,s)
        np.testing.assert_allclose(d['sector_intrinsic'],d['intrinsic'],atol=1e-13)
        self.assertGreater(abs(d['intrinsic'][0,1,0]),.1)

    def test_intrinsic_direct_finite_field_oracle(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        d=intrinsic_spin_geometry(h,v,s)
        for a in range(2):
            for b in range(2):
                self.assertAlmostEqual(finite_eigenstate_oracle(h,v,s,a,b),d['intrinsic'][0,a,b],delta=3e-9)

    def test_intrinsic_spectral_manual_fourband(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        E,U=np.linalg.eigh(h);O=spin_current_operator(v,s)
        d=intrinsic_spin_geometry(h,v,s);raw=np.zeros((4,2,2))
        for n in range(4):
            for m in range(4):
                if n!=m:
                    for a in range(2):
                        for b in range(2):
                            raw[n,a,b]-=2*np.imag(np.vdot(U[:,n],O[a]@U[:,m])*np.vdot(U[:,m],v[b]@U[:,n]))/(E[n]-E[m])**2
        np.testing.assert_allclose(raw,d['intrinsic'],atol=1e-12)

    def test_linear_hessian_zero_and_full_reduces_geometry(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22,periodic=False)
        d=spin_current_derivative(h,v,K,s)
        np.testing.assert_allclose(d['hessian_vertex'],0.,atol=1e-14)
        np.testing.assert_allclose(d['full_RTA'],d['geometry_only'],atol=1e-14)

    def test_linear_finite_difference_matches_geometry_only(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22,periodic=False)
        d=spin_current_derivative(h,v,K,s)
        for a in range(2):
            for b in range(2):
                np.testing.assert_allclose(numeric_band_derivative(.31,-.22,a,b,False),-d['geometry_only'][:,a,b],atol=1e-8)

    def test_periodic_finite_difference_requires_hessian(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        d=spin_current_derivative(h,v,K,s)
        for a in range(2):
            for b in range(2):
                np.testing.assert_allclose(numeric_band_derivative(.31,-.22,a,b,True),-d['full_RTA'][:,a,b],atol=1e-8)
        self.assertGreater(np.max(abs(d['hessian_vertex'])),.1)

    def test_bz_fermi_surface_integration_by_parts(self):
        r=integrate_spin_rta_bz(mesh=40)
        np.testing.assert_allclose(r['direct_fermi_surface'],r['full_fermi_sea'],atol=5e-8)

    def test_bz_nonzero_missing_hessian_not_mesh_error(self):
        r=integrate_spin_rta_bz(mesh=40)
        self.assertGreater(abs(r['missing_hessian'][0,0]),.005)
        self.assertGreater(abs(r['geometry_only'][0,0]),abs(r['full_fermi_sea'][0,0]))

    def test_bz_reference_frozen(self):
        r=integrate_spin_rta_bz(mesh=64)
        target=json.loads((HERE/'results/frozen/round13_spin_appendix_c.json').read_text())
        for key,ref in target['N64_spin_xz_E_x'].items():
            self.assertAlmostEqual(r[key][0,0],ref,delta=2e-11)

    def test_unitary_spin_orbital_basis_covariance(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        rng=np.random.default_rng(15)
        Z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        U,_=np.linalg.qr(Z);rotate=lambda A:U.conj().T@A@U
        d=spin_current_derivative(h,v,K,s)
        e=spin_current_derivative(rotate(h),np.stack([rotate(x) for x in v]),
                                  np.stack([[rotate(x) for x in row] for row in K]),rotate(s))
        for key in ('intrinsic','geometry_only','hessian_vertex','full_RTA'):
            np.testing.assert_allclose(d[key],e[key],atol=1e-12)

    def test_two_band_linear_spin_current_interband_vanishes(self):
        h=.22*PAULI[0]-.3*PAULI[1]+.51*PAULI[2]
        v=np.stack([PAULI[0],PAULI[1]])
        for s in PAULI:
            d=intrinsic_spin_geometry(h,v,s)
            self.assertLess(np.max(abs(d['intrinsic'])),1e-14)

    def test_nonphysical_spin_is_not_implicitly_inferred(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        with self.assertRaisesRegex(ValueError,'physical Pauli spin'):
            intrinsic_spin_geometry(h,v,np.eye(2))

    def test_missing_hessian_fails_closed(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        with self.assertRaisesRegex(ValueError,'Hessian'):
            spin_current_derivative(h,v,None,s)

    def test_degeneracy_fails_closed(self):
        h=np.diag([0.,0.,1.,2.]);v=np.zeros((2,4,4));s=np.eye(4)
        with self.assertRaisesRegex(ValueError,'Degenerate'):
            intrinsic_spin_geometry(h,v,s)

    def test_nonhermitian_and_mixed_hessian_fail_closed(self):
        h,v,K,s=spin_orbital_benchmark(.31,-.22)
        broken=v.copy();broken[0,0,1]+=.12j
        with self.assertRaisesRegex(ValueError,'Hermitian'):
            spin_current_derivative(h,broken,K,s)
        brokenK=K.copy();brokenK[0,1]+=np.eye(4)
        with self.assertRaisesRegex(ValueError,'symmetric'):
            spin_current_derivative(h,v,brokenK,s)

    def test_temperature_zero_fails_for_fermi_surface(self):
        with self.assertRaises(ValueError):fermi(np.arange(4),temperature=0)

if __name__=='__main__':unittest.main()
