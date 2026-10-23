"""Round33 actual native Wannier90 dense ASCII H and r parsing (no real source in repo)."""
import shutil
import tempfile
import unittest
from pathlib import Path
import numpy as np
from flux_lab.wannier90_raw_gate import parse_hr,parse_r,staged_h_r_import

def fixture(directory,deg=(2,1,2),expanded=False):
    """Independently construct documented original Wannier90 ASCII grammar."""
    p=Path(directory);p.mkdir(parents=True,exist_ok=True)
    R=[(-1,0,0),(0,0,0),(1,0,0)]
    H0=np.array([[-.28,.055+.02j],[.055-.02j,.42]])
    Hp=np.array([[.11,.038+.013j],[.02-.044j,.08]])
    Hm=Hp.conj().T
    A0=np.zeros((3,2,2),complex)
    A0[0]=np.diag([.14,.39]);A0[1]=np.diag([.02,.23]);A0[2]=np.diag([.0,.0])
    Ap=np.zeros_like(A0)
    Ap[0]=np.array([[.017,.012j],[-.021j,.009]])
    Ap[1]=np.array([[.006,-.003j],[.004j,.01]])
    Ap[2]=np.diag([.001,.002])
    Am=np.swapaxes(Ap.conj(),-1,-2)
    heff=np.stack([Hm,H0,Hp])
    aeff=np.stack([Am,A0,Ap])
    deg=np.asarray(deg,int)
    hraw=heff if expanded else heff*deg[:,None,None]
    araw=aeff if expanded else aeff*deg[:,None,None,None]
    hr=["fixture generated for test only","2","3"," ".join(str(x) for x in deg)]
    rr=["fixture generated for test only","2","3"]
    for i,Rv in enumerate(R):
        for m in range(2):
            for n in range(2):
                z=hraw[i,m,n]
                hr.append(f"{Rv[0]} {Rv[1]} {Rv[2]} {m+1} {n+1} {z.real:.16e} {z.imag:.16e}")
                comps=araw[i,:,m,n]
                vals=" ".join(f"{v:.16e}" for z in comps for v in (z.real,z.imag))
                rr.append(f"{Rv[0]} {Rv[1]} {Rv[2]} {m+1} {n+1} {vals}")
    (p/"source_hr.dat").write_text("\n".join(hr)+"\n")
    (p/"source_r.dat").write_text("\n".join(rr)+"\n")
    return dict(hr=p/"source_hr.dat",r=p/"source_r.dat",heff=heff,aeff=aeff,
                hraw=hraw,araw=araw,centers=np.array([[.14,.02,0],[.39,.23,0]]))

def import_pair(f,*,expanded=False,ws=False,has_wsvec=False,transl=False):
    return staged_h_r_import(f["hr"],f["r"],use_ws_distance=ws,
               write_ndegen_applied=expanded,wsvec_present=has_wsvec,
               transl_inv_full=transl)

def rewrite(p,fn):
    lines=p.read_text().splitlines()
    p.write_text("\n".join(fn(lines))+"\n")

class NativeW90HREvidence(unittest.TestCase):
    def test_two_native_records_and_matched_wigner_seitz_weights(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            h=parse_hr(f["hr"]);r=parse_r(f["r"])
            self.assertEqual(h["nwann"],2)
            self.assertEqual(h["nrpts"],3)
            np.testing.assert_array_equal(h["degeneracy"],[2,1,2])
            np.testing.assert_allclose(h["raw"],f["hraw"],atol=1e-14)
            np.testing.assert_allclose(r["raw"],f["araw"],atol=1e-14)
            result=import_pair(f)
            np.testing.assert_allclose(result.h_r,f["heff"],atol=1e-14)
            np.testing.assert_allclose(result.position_r,f["aeff"],atol=1e-14)
            np.testing.assert_allclose(result.centers,f["centers"],atol=1e-14)

    def test_independent_fourier_vs_analytic_cosine_model(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);bundle=import_pair(f)
            q=.217
            H,A=bundle.fourier([q,.13])
            phase=np.exp(2j*np.pi*q)
            expected=f["heff"][1]+phase*f["heff"][2]+phase.conjugate()*f["heff"][0]
            target=f["aeff"][1]+phase*f["aeff"][2]+phase.conjugate()*f["aeff"][0]
            np.testing.assert_allclose(H,expected,atol=2e-14)
            np.testing.assert_allclose(A,target,atol=2e-14)
            np.testing.assert_allclose(H,H.conj().T,atol=2e-14)
            for ai in A:np.testing.assert_allclose(ai,ai.conj().T,atol=2e-14)

    def test_wrong_double_degeneracy_weight_is_measurably_wrong(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);b=import_pair(f)
            good,_=b.fourier([.17,0])
            R=b.R[:,:2]
            bad=np.einsum("r,rij->ij",np.exp(2j*np.pi*R@np.array([.17,0])),f["hraw"])
            self.assertGreater(np.max(abs(good-bad)),.05)

    def test_expanded_native_output_no_double_division(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td,deg=(1,1,1),expanded=True)
            b=import_pair(f,expanded=True,ws=True,has_wsvec=True,transl=True)
            np.testing.assert_allclose(b.h_r,f["heff"],atol=2e-14)
            np.testing.assert_allclose(b.position_r,f["aeff"],atol=2e-14)
            self.assertEqual(b.provenance_status,"PARSED_H_R_SPIN_NOT_AUTHENTICATED")

    def test_expanded_source_with_nonunit_weights_fail_closed(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            with self.assertRaisesRegex(ValueError,"requires all"):
                import_pair(f,expanded=True)

    def test_unexpanded_ws_distance_requires_explicit_wsvec_adapter(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            with self.assertRaisesRegex(ValueError,"unexpanded ws_distance"):
                import_pair(f,ws=True)
            with self.assertRaisesRegex(ValueError,"wsvec"):
                import_pair(f,has_wsvec=True)

    def test_source_flags_must_be_boolean_not_truthy_strings(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            with self.assertRaisesRegex(ValueError,"explicit boolean"):
                staged_h_r_import(f["hr"],f["r"],
                    use_ws_distance="false",write_ndegen_applied=False)

    def test_negative_weights_and_wrong_degeneracy_row(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            rewrite(f["hr"],lambda s:s[:3]+["2 0 2"]+s[4:])
            with self.assertRaisesRegex(ValueError,"must be positive"):parse_hr(f["hr"])
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            rewrite(f["hr"],lambda s:s[:3]+["2 1 2 4"]+s[4:])
            with self.assertRaisesRegex(ValueError,"degeneracy length"):parse_hr(f["hr"])

    def test_missing_matrix_entry_detected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            rewrite(f["hr"],lambda s:s[:-1])
            with self.assertRaisesRegex(ValueError,"incomplete/extra"):parse_hr(f["hr"])

    def test_duplicate_orbital_entry_detected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            rewrite(f["r"],lambda s:s[:4]+[s[3]]+s[5:])
            with self.assertRaisesRegex(ValueError,"duplicate matrix entry"):parse_r(f["r"])

    def test_mixed_translation_and_nonplanar_R_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            rewrite(f["hr"],lambda s:s[:5]+[s[5].replace("-1 0 0","2 0 0")]+s[6:])
            with self.assertRaisesRegex(ValueError,"mixed translations"):parse_hr(f["hr"])
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            rewrite(f["hr"],lambda s:s[:4]+[s[4].replace("-1 0 0","-1 0 1")]+s[5:])
            with self.assertRaisesRegex(ValueError,"Rz=0"):parse_hr(f["hr"])

    def test_hr_r_translation_order_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            rewrite(f["r"],lambda s:s[:3]+s[7:11]+s[3:7]+s[11:])
            with self.assertRaisesRegex(ValueError,"translations"):import_pair(f)

    def test_hermitian_reverse_translation_failure_detected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            def mutate(lines):
                tokens=lines[4].split();tokens[5]=str(float(tokens[5])+.125)
                lines[4]=" ".join(tokens);return lines
            rewrite(f["hr"],mutate)
            with self.assertRaisesRegex(ValueError,"Hermitian"):import_pair(f)

    def test_offsite_r_position_pair_failure_detected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            def mutate(lines):
                tokens=lines[3].split();tokens[5]=str(float(tokens[5])+.025)
                lines[3]=" ".join(tokens);return lines
            rewrite(f["r"],mutate)
            with self.assertRaisesRegex(ValueError,"Hermitian"):import_pair(f)

    def test_nonfinite_h_r_values_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td)
            rewrite(f["hr"],lambda s:s[:4]+[s[4][:-22]+" nan"]+s[5:])
            with self.assertRaises(ValueError):parse_hr(f["hr"])

    def test_real_material_physical_spin_remains_unverified(self):
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);x=import_pair(f)
            self.assertFalse(x.physical_spin_available)
            self.assertFalse(x.kme_ready)
            self.assertEqual(x.original_paper_status,"PARTIAL")
            self.assertRegex(x.provenance_status,"SPIN_NOT_AUTHENTICATED")
            self.assertEqual(len(x.h_sha256),64)
            self.assertEqual(len(x.r_sha256),64)

    def test_rejects_malformed_fractional_k(self):
        with tempfile.TemporaryDirectory() as td:
            x=import_pair(fixture(td))
            with self.assertRaisesRegex(ValueError,"fractional 2D"):
                x.fourier([.1])
            with self.assertRaisesRegex(ValueError,"fractional 2D"):
                x.fourier([float("nan"),.3])

    def test_frozen_round33_first_remote_raw_file_format_numbers(self):
        import json
        raw=json.loads((Path(__file__).resolve().parents[1]/
            "results/frozen/round33_native_wannier90_hr_r_audit.json").read_text())
        self.assertEqual(raw["schema"],"cui-round33-w90-native-hr-r-degeneracy-spin-blocked-v1")
        self.assertEqual(raw["original_paper_status"],"PARTIAL")
        self.assertFalse(raw["authentic_material_input_found"])
        self.assertFalse(raw["wannier_physical_spin_provenance_authenticated"])
        self.assertFalse(raw["kme_real_material_ready"])
        with tempfile.TemporaryDirectory() as td:
            f=fixture(td);b=import_pair(f);H,A=b.fourier([.217,.13])
            phase=np.exp(2j*np.pi*.217)
            H_analytic=f["heff"][1]+phase*f["heff"][2]+phase.conjugate()*f["heff"][0]
            Re=np.exp(2j*np.pi*(b.R[:,:2]@np.array([.17,0])))
            missing=np.einsum("r,rij->ij",Re,f["hraw"])
            corrected=b.fourier([.17,0])[0]
            calculated={
              "native_weighted_h_reconstruction_max":float(np.max(abs(b.h_r-f["heff"]))),
              "native_weighted_r_reconstruction_max":float(np.max(abs(b.position_r-f["aeff"]))),
              "analytic_complex_fourier_H_max":float(np.max(abs(H-H_analytic))),
              "wrong_omitted_ws_degeneracy_H_max":float(np.max(abs(missing-corrected))),
            }
            for key,val in calculated.items():
                self.assertAlmostEqual(val,raw[key],delta=2e-11,msg=key)
            self.assertEqual(list(b.degeneracy),raw["wigner_seitz_degeneracy"])
            self.assertFalse(b.physical_spin_available)

if __name__=="__main__":unittest.main()
