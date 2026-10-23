"""Round32 fail-closed external scattering data/gauge/provenance eligibility tests.

No real-material EPW files are present; fixtures are explicitly synthetic.
"""
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from flux_lab.epw_material_contract import (
 SCHEMA,UNITS,ROLES,make_synthetic_fixture,validate_normalized,
 scatter_rates,rotate_electronic_gauge,load_bundle)

def setup():
 a,m=make_synthetic_fixture()
 verdict=validate_normalized(a,m)
 return a,m,verdict

def fixture_folder(tmp,a,m):
 p=Path(tmp);np.savez(p/"normalized.npz",**a)
 m=dict(m);m["normalized_sha256"]=hashlib.sha256((p/"normalized.npz").read_bytes()).hexdigest()
 (p/"manifest.json").write_text(json.dumps(m))
 return p

class EpwMaterialEvidenceGate(unittest.TestCase):
 def test_synthetic_truthful_status_and_explicit_array_shapes(self):
  a,m,v=setup()
  self.assertEqual(v['status'],"SYNTHETIC_FIXTURE_ONLY")
  self.assertFalse(v['real_material_verified'])
  self.assertEqual(v['paper_status'],"PARTIAL")
  self.assertEqual((v['nk'],v['nq'],v['nband'],v['nmode']),(3,3,2,1))
  self.assertEqual(v['active_mode_count'],3)

 def test_exact_fractional_k_plus_q_and_inverse_q(self):
  a,m,v=setup()
  for ik in range(3):
   for iq in range(3):
    self.assertEqual(v['kplusq'][ik,iq],(ik+iq)%3)
  np.testing.assert_array_equal(v['q_inverse'],[0,2,1])

 def test_normalized_archive_integrity_and_no_pickle(self):
  a,m,v=setup()
  with tempfile.TemporaryDirectory() as td:
   p=fixture_folder(td,a,m)
   data,meta,result=load_bundle(p)
   self.assertEqual(result['status'],"SYNTHETIC_FIXTURE_ONLY")
   np.testing.assert_allclose(data['g_eV'],a['g_eV'])
   with open(p/'normalized.npz','ab') as f:f.write(b'corrupt')
   with self.assertRaisesRegex(ValueError,"SHA256"):load_bundle(p)

 def test_missing_full_source_never_claims_real(self):
  a,m,_=setup();m['stage']='ab_initio_candidate'
  m['exporter_version']='explicit_external_mapping'
  with self.assertRaisesRegex(ValueError,"source provenance"):
   validate_normalized(a,m)
  with tempfile.TemporaryDirectory() as td:
   with self.assertRaisesRegex(ValueError,"complete original source roles"):
    validate_normalized(a,m,folder=td)

 def test_candidate_hashes_are_integrity_not_ab_initio_authentication(self):
  a,m,_=setup()
  with tempfile.TemporaryDirectory() as td:
   root=Path(td)
   rows=[]
   # THESE ARE ARTIFICIAL BYTES, intended to prove the limit of checksums.
   for role in sorted(ROLES):
    name=role+'.raw'
    data=('FAKE_SYNTHETIC_EVIDENCE_'+role).encode()
    (root/name).write_bytes(data)
    rows.append(dict(role=role,path=name,sha256=hashlib.sha256(data).hexdigest()))
   m['stage']='ab_initio_candidate'
   m['exporter_version']='TEST_ONLY_fake_export'
   m['evidence']=rows
   v=validate_normalized(a,m,folder=root)
   self.assertIn("PHYSICAL_AUTHENTICITY_UNVERIFIED",v['status'])
   self.assertFalse(v['real_material_verified'])
   self.assertEqual(v['paper_status'],"PARTIAL")
   (root/rows[0]['path']).write_text("tampered")
   with self.assertRaisesRegex(ValueError,"SHA256"):validate_normalized(a,m,folder=root)

 def test_source_path_escape_rejected_even_with_valid_hash(self):
  a,m,_=setup()
  with tempfile.TemporaryDirectory() as td:
   m['stage']='ab_initio_candidate';m['exporter_version']='fake'
   m['evidence']=[dict(role=role,path='../outside.dat',sha256='0'*64)
                  for role in sorted(ROLES)]
   with self.assertRaises(ValueError):validate_normalized(a,m,folder=td)

 def test_phonon_negative_active_energy_and_acoustic_mask(self):
  a,m,_=setup()
  bad=copy.deepcopy(a);bad["phonon_eV"][0,0]=0.
  with self.assertRaisesRegex(ValueError,"nonpositive active"):
   validate_normalized(bad,m)
  masked=copy.deepcopy(a)
  masked['phonon_active'][:,0]=False
  masked['phonon_eV'][:,0]=0.
  masked['g_eV'][:]=0.
  val=validate_normalized(masked,m)
  self.assertEqual(val['active_mode_count'],0)
  bad=copy.deepcopy(masked);bad['g_eV'][0,0,0,0,0]=1
  with self.assertRaisesRegex(ValueError,"acoustic"):validate_normalized(bad,m)

 def test_inverse_q_phonon_phase_missing_rejected(self):
  a,m,_=setup()
  bad=copy.deepcopy(a);bad['g_eV'][0,1,0,0,1]+=.05
  with self.assertRaisesRegex(ValueError,"inverse q Hermiticity"):
   validate_normalized(bad,m)
  bad=copy.deepcopy(a);bad['phonon_eV'][1,0]+=.02
  with self.assertRaisesRegex(ValueError,"branch frequencies"):
   validate_normalized(bad,m)

 def test_missing_q_inverse_and_k_plus_q_mesh_rejected(self):
  a,m,_=setup()
  bad=copy.deepcopy(a);bad['q_frac'][2]=[.8,0]
  with self.assertRaisesRegex(ValueError,"inverse phonon"):
   validate_normalized(bad,m)
  bad=copy.deepcopy(a);bad['k_frac'][2]=[.8,0]
  with self.assertRaisesRegex(ValueError,"k\\+q"):
   validate_normalized(bad,m)

 def test_gauge_id_spin_kind_position_embedding_and_units_rejected(self):
  a,m,_=setup()
  bad=copy.deepcopy(m);bad['basis']['gauge_ids']['g']="ANOTHER"
  with self.assertRaisesRegex(ValueError,"shared electronic"):validate_normalized(a,bad)
  bad=copy.deepcopy(m);bad['basis']['spin_kind']='orbital_pseudospin'
  with self.assertRaisesRegex(ValueError,"pseudospin"):validate_normalized(a,bad)
  bad=copy.deepcopy(m);bad['basis']['position_convention']='centers_added_twice'
  with self.assertRaisesRegex(ValueError,"centers exactly once"):validate_normalized(a,bad)
  bad=copy.deepcopy(m);bad['units']['coupling']='Ry'
  with self.assertRaisesRegex(ValueError,"units"):validate_normalized(a,bad)
  bad=copy.deepcopy(m);bad['basis']['phonon_phase']='arbitrary_mode'
  with self.assertRaisesRegex(ValueError,"phase"):validate_normalized(a,bad)

 def test_nonhermitian_physical_spin_nonfinite_input_rejected(self):
  a,m,_=setup()
  bad=copy.deepcopy(a);bad['spin_pauli'][0,2,0,0]=1.3
  with self.assertRaisesRegex(ValueError,"exceeds one"):validate_normalized(bad,m)
  bad=copy.deepcopy(a);bad['h_eV'][0,0,1]=.03j
  with self.assertRaisesRegex(ValueError,"not Hermitian"):validate_normalized(bad,m)
  bad=copy.deepcopy(a);bad['g_eV'][0,0,0,0,0]=np.nan
  with self.assertRaisesRegex(ValueError,"nonfinite"):validate_normalized(bad,m)

 def test_scattering_unit_positivity_and_pauli_blocked_rates(self):
  a,m,v=setup()
  r=scatter_rates(a,v,beta=40,mu=0,sigma_eV=.08)
  self.assertEqual(r['unit'],"s^-1")
  self.assertFalse(r['real_material_verified'])
  self.assertFalse(r['finite_broadening_exact_KMS'])
  self.assertEqual(r['bare'].shape,(2,3,3,1,2,2))
  self.assertGreater(np.max(r['bare']),1e10)
  self.assertTrue(np.all(r['blocked']>=0))
  self.assertTrue(np.all(r['blocked']<=r['bare']*(1+1e-14)))
  self.assertGreater(np.max(r['bare']-r['blocked']),1e8)

 def test_explicit_pauli_blocking_vs_bare_and_temperature(self):
  a,m,v=setup()
  cold=scatter_rates(a,v,beta=140,mu=0,sigma_eV=.08)
  hot=scatter_rates(a,v,beta=10,mu=0,sigma_eV=.08)
  self.assertGreater(np.sum(hot['bare'][0]),np.sum(cold['bare'][0]))
  self.assertGreater(np.max(abs(hot['fermi']-cold['fermi'])),.05)

 def test_independent_k_dependent_unitary_of_h_and_g_rates(self):
  a,m,v=setup()
  sy=np.array([[0,-1j],[1j,0]],complex)
  us=np.array([np.cos(.27*(i+1))*np.eye(2)-1j*np.sin(.27*(i+1))*sy for i in range(3)])
  rotated=rotate_electronic_gauge(a,us,v['kplusq'])
  # Only H/g are rotated; DO NOT claim transformed S or A connection.
  original=scatter_rates(a,v,beta=40,mu=0,sigma_eV=.08)
  after=scatter_rates(rotated,v,beta=40,mu=0,sigma_eV=.08)
  np.testing.assert_allclose(original['bare'],after['bare'],atol=.03,rtol=2e-12)
  np.testing.assert_allclose(original['blocked'],after['blocked'],atol=.03,rtol=2e-12)
  # Position connection was intentionally not transformed by independent k-gauges.

 def test_bad_unmatched_gauge_changes_rates_and_is_rejected(self):
  a,m,v=setup()
  sy=np.array([[0,-1j],[1j,0]],complex)
  us=np.array([np.cos(.3*(i+1))*np.eye(2)-1j*np.sin(.3*(i+1))*sy for i in range(3)])
  corrupt=copy.deepcopy(a)
  matched=rotate_electronic_gauge(a,us,v['kplusq'])
  corrupt['h_eV']=matched['h_eV']
  base=scatter_rates(a,v,beta=40,mu=0,sigma_eV=.08)
  wrong=scatter_rates(corrupt,v,beta=40,mu=0,sigma_eV=.08)
  self.assertGreater(np.max(abs(base['bare']-wrong['bare'])),1e8)

 def test_broadening_and_temperature_fail_closed(self):
  a,m,v=setup()
  for kwargs in (dict(beta=0),dict(sigma_eV=0),dict(mu=float('nan'))):
   with self.assertRaises(ValueError):scatter_rates(a,v,**kwargs)

 def test_source_status_cannot_be_promoted_by_metadata_claim(self):
  a,m,_=setup()
  m['stage']='synthetic_fixture';m['material_id']='Hg3O3'
  verdict=validate_normalized(a,m)
  self.assertEqual(verdict['status'],"SYNTHETIC_FIXTURE_ONLY")
  self.assertFalse(verdict['real_material_verified'])

 def test_frozen_remote_round32_bounded_synthetic_regression(self):
  from pathlib import Path
  raw=json.loads((Path(__file__).resolve().parents[1]/
     "results/frozen/round32_epw_material_scattering_contract.json").read_text())
  self.assertEqual(raw["schema"],"cui-round32-epw-normalized-provenance-kq-gauge-goldenrule-v1")
  self.assertEqual(raw["whole_source_paper"],"PARTIAL")
  self.assertFalse(raw["real_material_verified"])
  self.assertFalse(raw["authentic_real_epw_phonon_matrix_found"])
  self.assertFalse(raw["physical_spin_provenance_authenticated"])
  a,m,v=setup()
  r=scatter_rates(a,v,beta=40,mu=0,sigma_eV=.08)
  sy=np.array([[0,-1j],[1j,0]],complex)
  U=np.array([np.cos(.27*(i+1))*np.eye(2)-1j*np.sin(.27*(i+1))*sy for i in range(3)])
  rotated=rotate_electronic_gauge(a,U,v["kplusq"])
  matched=scatter_rates(rotated,v,beta=40,mu=0,sigma_eV=.08)
  bad=copy.deepcopy(a);bad["h_eV"]=rotated["h_eV"]
  unmatched=scatter_rates(bad,v,beta=40,mu=0,sigma_eV=.08)
  measurements={
   "bare_finite_gaussian_goldenrule_rate_s_inverse_max":
      float(np.max(r["bare"])),
   "pauli_blocking_rate_s_inverse_difference_max":
      float(np.max(r["bare"]-r["blocked"])),
   "matched_independent_k_gauge_goldenrule_rate_s_inverse_residual_max":
      float(np.max(abs(r["bare"]-matched["bare"]))),
   "wrong_H_only_gauge_goldenrule_rate_s_inverse_difference_max":
      float(np.max(abs(r["bare"]-unmatched["bare"]))),
  }
  for key,val in measurements.items():
   self.assertAlmostEqual(val,raw[key],delta=1e-9*max(1,abs(raw[key])),msg=key)
  self.assertEqual(raw["fake_hashed_candidate_status"],
                   "CANDIDATE_BYTES_VALIDATED_PHYSICAL_AUTHENTICITY_UNVERIFIED")

if __name__=='__main__':unittest.main()
