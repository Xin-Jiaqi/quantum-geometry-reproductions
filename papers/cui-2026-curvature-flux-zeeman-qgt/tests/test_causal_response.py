"""Physical/source frequency, four-sector and independent Liouville acceptance."""
import json
import unittest
from pathlib import Path
import numpy as np
from flux_lab.core import PAULI
from flux_lab.causal_response import (
    _data,interband_kubo,interband_liouville,static_subtracted_gyrotropic,
    reciprocal_residual,source_four_sector_expansion,source_four_sector_geometric)

ROOT=Path(__file__).resolve().parents[1]

def fixture(kx=.37,ky=-.29,tilt=.23):
    # Synthetic |orbital> x |physical spin>; no real-material claim.
    spin=np.stack([np.kron(np.eye(2),s) for s in PAULI])
    tau=np.stack([np.kron(s,np.eye(2)) for s in PAULI])
    mass=.28+.65*(2-np.cos(kx)-np.cos(ky))
    H=(np.sin(kx)*spin[0]+np.sin(ky)*spin[1]+mass*spin[2]
       +2.5*tau[2]+.7*(tau[0]@spin[2])+tilt*np.sin(kx)*np.eye(4))
    V=np.stack([np.cos(kx)*spin[0]+.65*np.sin(kx)*spin[2]+tilt*np.cos(kx)*np.eye(4),
                np.cos(ky)*spin[1]+.65*np.sin(ky)*spin[2]])
    return H,V,spin

class CausalResponse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.H,cls.V,cls.S=fixture()
        cls.kw=dict(mu=2.3,temperature=.025,kappa=1.)

    def test_liouville_against_spectral(self):
        for z in (.035,-.071,.035+.003j,.045+.01j):
            a=interband_kubo(self.H,self.V,self.S,z=z,**self.kw)
            b=interband_liouville(self.H,self.V,self.S,z=z,**self.kw)
            for x,y in zip(a,b):self.assertLess(np.max(abs(x-y)),4e-13)

    def test_exact_reciprocity_with_complex_regulator(self):
        for w,eta in [(.035,0.),(-.047,0.),(.037,.003),(.037,.04)]:
            q=reciprocal_residual(self.H,self.V,self.S,omega_drive=w,eta=eta,**self.kw)
            self.assertLess(q["max_residual"],3e-13)

    def test_omega_cubic_truncation(self):
        errors=[]
        for w in [.08,.04,.02]:
            a=source_four_sector_expansion(self.H,self.V,self.S,omega_source=w,**self.kw)
            b=static_subtracted_gyrotropic(self.H,self.V,self.S,omega_drive=-w,**self.kw)
            errors.append(np.max(abs(a["sigma_source_omega"]-b)))
        self.assertGreater(errors[0]/errors[1],7.5)
        self.assertGreater(errors[1]/errors[2],7.5)

    def test_explicit_four_sectors_match_pairwise_source(self):
        for w in [.007,.025,-.017]:
            a=source_four_sector_expansion(self.H,self.V,self.S,omega_source=w,**self.kw)
            b=source_four_sector_geometric(self.H,self.V,self.S,omega_source=w,**self.kw)
            np.testing.assert_allclose(a["sigma_source_omega"],b["sigma_source_omega"],atol=3e-14)
            for k in ("omegaA","gN","omegaN","gA"):
                self.assertGreater(np.max(abs(b[k])),1e-7)

    def test_source_fourier_sign(self):
        w=.008
        expansion=source_four_sector_expansion(self.H,self.V,self.S,omega_source=w,**self.kw)
        correct=static_subtracted_gyrotropic(self.H,self.V,self.S,omega_drive=-w,**self.kw)
        wrong=static_subtracted_gyrotropic(self.H,self.V,self.S,omega_drive=w,**self.kw)
        self.assertLess(np.max(abs(expansion["sigma_source_omega"]-correct)),2e-6)
        self.assertGreater(np.max(abs(expansion["sigma_source_omega"]-wrong)),1e-3)

    def test_odd_and_even_sector_parity(self):
        a=source_four_sector_geometric(self.H,self.V,self.S,omega_source=.03,**self.kw)
        b=source_four_sector_geometric(self.H,self.V,self.S,omega_source=-.03,**self.kw)
        for k in ("omegaA","gN"):
            np.testing.assert_allclose(a[k],-b[k],atol=3e-14)
            self.assertLess(np.max(abs(a[k].real)),1e-12)
        for k in ("omegaN","gA"):
            np.testing.assert_allclose(a[k],b[k],atol=3e-14)
            self.assertLess(np.max(abs(a[k].imag)),1e-12)

    def test_static_kme_cannot_be_deleted_from_exact_identity(self):
        q=reciprocal_residual(self.H,self.V,self.S,omega_drive=.035,**self.kw)
        self.assertGreater(q["max_wrong_double_subtraction"],1e-2)
        self.assertGreater(np.max(abs(q["alpha_static"])),.5)

    def test_wrong_retarded_kme_branch_detected(self):
        w=.035;eta=.02;z=w+1j*eta;z0=1j*eta
        actual=static_subtracted_gyrotropic(self.H,self.V,self.S,omega_drive=w,eta=eta,**self.kw)
        wrongA=interband_kubo(self.H,self.V,self.S,z=-w+1j*eta,**self.kw)[1]
        wrong0=interband_kubo(self.H,self.V,self.S,z=z0,**self.kw)[1]
        self.assertGreater(np.max(abs(actual-1j*(z*wrongA-z0*wrong0).T)),1e-5)

    def test_four_band_constant_basis_unitary(self):
        rng=np.random.default_rng(26)
        X=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        U,_=np.linalg.qr(X)
        tr=lambda A:U.conj().T@A@U
        a=interband_kubo(self.H,self.V,self.S,z=.037,**self.kw)
        b=interband_kubo(tr(self.H),np.stack([tr(v) for v in self.V]),
                          np.stack([tr(s) for s in self.S]),z=.037,**self.kw)
        for x,y in zip(a,b):np.testing.assert_allclose(x,y,atol=5e-13)

    def test_wrong_spin_basis_is_physically_different(self):
        rng=np.random.default_rng(26)
        X=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        U,_=np.linalg.qr(X)
        tr=lambda A:U.conj().T@A@U
        a=interband_kubo(self.H,self.V,self.S,z=.037,**self.kw)[0]
        wrong=interband_kubo(tr(self.H),np.stack([tr(v) for v in self.V]),
                             self.S,z=.037,**self.kw)[0]
        self.assertGreater(np.max(abs(a-wrong)),1e-4)

    def test_interband_position_hermiticity(self):
        R=_data(self.H,self.V,self.S)["R"]
        for r in R:np.testing.assert_allclose(r,r.conj().T,atol=1e-13)

    def test_insulating_and_metallic_occupations(self):
        for mu in [0.,2.3]:
            q=reciprocal_residual(self.H,self.V,self.S,omega_drive=.045,mu=mu,temperature=0.)
            self.assertLess(q["max_residual"],2e-13)

    def test_zeeman_vertex_scaling(self):
        a=static_subtracted_gyrotropic(self.H,self.V,self.S,omega_drive=.034,kappa=1.,mu=2.3)
        b=static_subtracted_gyrotropic(self.H,self.V,self.S,omega_drive=.034,kappa=.4,mu=2.3)
        np.testing.assert_allclose(b,.4*a,atol=5e-13)

    def test_zero_frequency_dynamic_response(self):
        a=static_subtracted_gyrotropic(self.H,self.V,self.S,omega_drive=0.,eta=.003,**self.kw)
        self.assertLess(np.max(abs(a)),2e-14)

    def test_degeneracy_missing_spin_and_hermiticity_fail_closed(self):
        with self.assertRaisesRegex(ValueError,"Degenerate"):
            interband_kubo(np.diag([0.,0.,1.,2.]),self.V,self.S)
        with self.assertRaises(ValueError):interband_kubo(self.H,self.V,self.S[:2])
        s=self.S.copy();s[0,0,1]+=.23j
        with self.assertRaises(ValueError):interband_kubo(self.H,self.V,s)

    def test_invalid_regulator_and_temperature(self):
        with self.assertRaises(ValueError):interband_kubo(self.H,self.V,self.S,temperature=-1)
        with self.assertRaises(ValueError):static_subtracted_gyrotropic(self.H,self.V,self.S,eta=-.01)

    def test_round8_frozen_reference(self):
        r=json.loads((ROOT/"results/frozen/round8_causal_public_v1.json").read_text())
        chi=static_subtracted_gyrotropic(self.H,self.V,self.S,omega_drive=r["omega_drive"],**self.kw)
        alpha=interband_kubo(self.H,self.V,self.S,z=0.,**self.kw)[1]
        self.assertAlmostEqual(chi[1,2].real,r["delta_chi_yz_real"],places=11)
        self.assertAlmostEqual(chi[1,2].imag,r["delta_chi_yz_imag"],places=11)
        self.assertAlmostEqual(alpha[2,1].real,r["static_kme_zy"],places=11)

    def test_two_band_pauli_control(self):
        H=.31*PAULI[0]-.17*PAULI[1]+.23*PAULI[2]
        V=np.stack([PAULI[0],PAULI[1]])
        a=interband_kubo(H,V,PAULI,mu=0.,temperature=0.,z=.021)
        b=interband_liouville(H,V,PAULI,mu=0.,temperature=0.,z=.021)
        for x,y in zip(a,b):np.testing.assert_allclose(x,y,atol=1e-12)

if __name__=="__main__":unittest.main()
