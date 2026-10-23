"""Round26: independent Wilson, projector, physical-spin and failure controls."""
import copy
import json
from pathlib import Path
import unittest
import numpy as np
from flux_lab.wannier_operator_contract import synthetic_soc_bundle,global_gauge_transform

Q=np.array([.183,-.271])

def projector_qgt(H,J,S):
    e,u=np.linalg.eigh(H)
    if np.min(np.diff(e))<1e-7:raise ValueError("nondegenerate oracle required")
    ans=np.zeros((len(e),3,3),complex)
    for n in range(len(e)):
        pn=np.outer(u[:,n],u[:,n].conj())
        for m in range(len(e)):
            if n==m:continue
            pm=np.outer(u[:,m],u[:,m].conj())
            for a in range(3):
                for b in range(3):
                    ans[n,a,b]+=1j*np.trace(pn@J[a]@pm@S[b])/(e[m]-e[n])
    return ans

def wilson_oracle(b,q,step=1e-5):
    """Independent finite difference of parallel-transported H and spin."""
    mid=b.fourier(q)
    J=np.zeros_like(mid["covariant_H"])
    DS=np.zeros_like(mid["covariant_spin"])
    for a in range(3):
        dq=b.lattice[:,a]*step/(2*np.pi)
        plus,minus=b.fourier(q+dq),b.fourier(q-dq)
        e,w=np.linalg.eigh(mid["connection"][a])
        U=(w*np.exp(-1j*step*e))@w.conj().T
        forward=lambda X:U@X@U.conj().T
        backward=lambda X:U.conj().T@X@U
        J[a]=(forward(plus["H"])-backward(minus["H"]))/(2*step)
        DS[a]=np.stack([(forward(plus["spin"][j])-backward(minus["spin"][j]))/(2*step)
                         for j in range(3)])
    return J,DS

class SameWannierGauge(unittest.TestCase):
    def setUp(self):
        self.b=synthetic_soc_bundle()

    def test_analytic_fourier_and_cartesian_gradient(self):
        phi=2*np.pi*Q[0];b=self.b
        d=b.fourier(Q)
        expected=b.h_r[1]+np.exp(1j*phi)*b.h_r[2]+np.exp(-1j*phi)*b.h_r[0]
        np.testing.assert_allclose(d["H"],expected,atol=1e-13)
        dx=1j*b.lattice[0,0]*(np.exp(1j*phi)*b.h_r[2]-np.exp(-1j*phi)*b.h_r[0])
        np.testing.assert_allclose(d["partial_H"][0],dx,atol=1e-13)
        np.testing.assert_allclose(d["partial_H"][1:],0,atol=1e-14)

    def test_independent_wilson_parallel_transport(self):
        d=self.b.fourier(Q);J,DS=wilson_oracle(self.b,Q)
        np.testing.assert_allclose(d["covariant_H"],J,atol=5e-9)
        np.testing.assert_allclose(d["covariant_spin"],DS,atol=5e-9)

    def test_position_connection_cannot_be_omitted(self):
        d=self.b.fourier(Q);J,_=wilson_oracle(self.b,Q)
        np.testing.assert_allclose(d["covariant_H"],J,atol=5e-9)
        self.assertGreater(np.max(abs(d["partial_H"]-J)),.005)

    def test_global_unitary_all_operator_and_qgt_covariance(self):
        rng=np.random.default_rng(26)
        z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        U,_=np.linalg.qr(z)
        d=self.b.fourier(Q)
        rot=global_gauge_transform(self.b,U).fourier(Q)
        for key in ("H","spin","connection","partial_H","partial_spin","covariant_H","covariant_spin"):
            np.testing.assert_allclose(rot[key],U.conj().T@d[key]@U,atol=2e-12)
        np.testing.assert_allclose(projector_qgt(d["H"],d["covariant_H"],d["spin"]),
                 projector_qgt(rot["H"],rot["covariant_H"],rot["spin"]),atol=3e-12)

    def test_wrong_spin_gauge_physics_fails_even_if_metadata_falsified(self):
        rng=np.random.default_rng(261)
        U,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
        d=self.b.fourier(Q);r=global_gauge_transform(self.b,U).fourier(Q)
        right=projector_qgt(r["H"],r["covariant_H"],r["spin"])
        wrong=projector_qgt(r["H"],r["covariant_H"],d["spin"])
        self.assertGreater(np.max(abs(wrong-right)),.02)

    def test_soc_real_spin_not_conserved(self):
        d=self.b.fourier(Q)
        torque=np.stack([d["H"]@s-s@d["H"] for s in d["spin"]])
        self.assertGreater(np.max(abs(torque)),.1)
        np.testing.assert_allclose(d["spin"],self.b.spin_r[1],atol=1e-13)

    def test_reciprocal_periodic_hermitian(self):
        d=self.b.fourier(Q)
        for shift in ((1,0),(0,-1),(1,-2)):
            other=self.b.fourier(Q+np.array(shift))
            for key in ("H","spin","connection","covariant_H","covariant_spin"):
                np.testing.assert_allclose(other[key],d[key],atol=2e-13)
        for key in ("H","spin","connection","covariant_H","covariant_spin"):
            x=d[key]
            np.testing.assert_allclose(x,np.swapaxes(x.conj(),-1,-2),atol=2e-13)

    def test_wrong_gauge_id_and_pseudospin_fail_closed(self):
        b=copy.deepcopy(self.b);b.provenance["basis_ids"]["spin"]="foreign"
        with self.assertRaisesRegex(ValueError,"same Wannier"):b.validate()
        b=copy.deepcopy(self.b);b.provenance["spin_kind"]="orbital_pseudospin"
        with self.assertRaisesRegex(ValueError,"physical Pauli"):b.validate()
        b=copy.deepcopy(self.b);b.provenance["spin_evidence"]=""
        with self.assertRaisesRegex(ValueError,"spin evidence"):b.validate()

    def test_units_centers_and_missing_physical_spin_fail_closed(self):
        b=copy.deepcopy(self.b);b.units["spin"]="hbar_over_two"
        with self.assertRaisesRegex(ValueError,"units"):b.validate()
        b=copy.deepcopy(self.b);b.centers[0,0]+=.1
        with self.assertRaisesRegex(ValueError,"centers"):b.validate()
        b=copy.deepcopy(self.b);b.spin_r=None
        with self.assertRaises((ValueError,TypeError)):b.validate()

    def test_hopping_hermiticity_inverse_translation_and_integer_check(self):
        b=copy.deepcopy(self.b);b.h_r[0,0,1]+=.02
        with self.assertRaisesRegex(ValueError,"Hermiticity"):b.validate()
        b=copy.deepcopy(self.b);b.translations[0]=[-2,0]
        with self.assertRaisesRegex(ValueError,"inverse"):b.validate()
        b=copy.deepcopy(self.b);b.translations=b.translations.astype(float);b.translations[0,0]=-.5
        with self.assertRaisesRegex(ValueError,"integer"):b.validate()

    def test_spin_norm_position_basis_and_provenance_guards(self):
        b=copy.deepcopy(self.b);b.spin_r[1]*=1.2
        with self.assertRaisesRegex(ValueError,"norm exceeds"):b.validate()
        b=copy.deepcopy(self.b);b.provenance["basis_ids"]["position"]="other"
        with self.assertRaisesRegex(ValueError,"same Wannier"):b.validate()

    def test_existing_geometry_solver_bridge_and_independent_oracle(self):
        from flux_lab.multiband import zeeman_qgt,covariant_h_derivative
        d=self.b.fourier(Q)
        np.testing.assert_allclose(covariant_h_derivative(d["H"],d["partial_H"],d["connection"]),
                                   d["covariant_H"],atol=1e-13)
        spectral=zeeman_qgt(d["H"],d["covariant_H"],d["spin"],method="spectral")
        project=zeeman_qgt(d["H"],d["covariant_H"],d["spin"],method="projector")
        np.testing.assert_allclose(spectral["tensor"],project["tensor"],atol=1e-12)
        np.testing.assert_allclose(project["tensor"],projector_qgt(d["H"],d["covariant_H"],d["spin"]),atol=1e-12)

    def test_round26_frozen_audit_regression(self):
        f=Path(__file__).resolve().parents[1]/"results/frozen/round26_same_gauge_operator_audit.json"
        ref=json.loads(f.read_text())
        self.assertEqual(ref["schema"],"cui-round26-same-wannier-operator-contract-v1")
        self.assertEqual(ref["paper_status"],"PARTIAL")
        d=self.b.fourier(ref["fractional_k"])
        self.assertAlmostEqual(np.max(abs(d["covariant_H"])),ref["covariant_vertex_abs_max"],delta=1e-12)
        self.assertAlmostEqual(np.max(abs(d["covariant_H"]-d["partial_H"])),
                               ref["missing_connection_error_max"],delta=1e-12)

if __name__=="__main__":unittest.main()
