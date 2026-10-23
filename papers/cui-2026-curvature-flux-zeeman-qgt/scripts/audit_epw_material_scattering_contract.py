"""Round32: normalized EPW data eligibility and no-fake-material audit.

Uses ONLY explicitly SYNTHETIC k/q/H/g/S/A arrays plus mock raw-file bytes
to demonstrate candidate-file integrity is NOT ab-initio authentication.
No authentic material e-ph matrix file has been provided in this repository.
"""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from flux_lab.epw_material_contract import (
    SCHEMA,ROLES,make_synthetic_fixture,validate_normalized,
    rotate_electronic_gauge,scatter_rates,load_bundle)

a,m=make_synthetic_fixture()
v=validate_normalized(a,m)
result=scatter_rates(a,v,beta=40.,mu=0,sigma_eV=.08)
sy=np.array([[0,-1j],[1j,0]],complex)
unitaries=np.array([np.cos(.27*(i+1))*np.eye(2)-1j*np.sin(.27*(i+1))*sy for i in range(3)])
rot=rotate_electronic_gauge(a,unitaries,v['kplusq'])
rot_rate=scatter_rates(rot,v,beta=40.,mu=0,sigma_eV=.08)
bad=copy.deepcopy(a);bad['h_eV']=rot['h_eV']
wrong=scatter_rates(bad,v,beta=40.,mu=0,sigma_eV=.08)
with tempfile.TemporaryDirectory() as td:
    folder=Path(td)
    np.savez(folder/'normalized.npz',**a)
    mm=copy.deepcopy(m)
    mm['normalized_sha256']=hashlib.sha256((folder/'normalized.npz').read_bytes()).hexdigest()
    (folder/'manifest.json').write_text(json.dumps(mm))
    _,_,loaded=load_bundle(folder)
    original_status=loaded['status']
    mm['stage']='ab_initio_candidate'
    mm['exporter_version']='MOCK_NOT_REAL'
    mm['evidence']=[]
    for role in sorted(ROLES):
        name=role+'.raw';datum=('EXPLICITLY_FAKE_'+role).encode()
        (folder/name).write_bytes(datum)
        mm['evidence'].append(dict(role=role,path=name,
                        sha256=hashlib.sha256(datum).hexdigest()))
    candidate=validate_normalized(a,mm,folder=folder)
    candidate_status=candidate['status']
    candidate_auth=candidate['real_material_verified']
    (folder/mm['evidence'][0]['path']).write_text('TAMPERED')
    tamper_rejected=False
    try:validate_normalized(a,mm,folder=folder)
    except ValueError as exc:tamper_rejected=('SHA256' in str(exc))
paired_res=max(float(np.max(abs(a['g_eV'][ik,iq]-
    np.swapaxes(a['g_eV'][v['kplusq'][ik,iq],v['q_inverse'][iq]].conj(),-1,-2))))
    for ik in range(3) for iq in range(3))
report={
 'schema':'cui-round32-epw-normalized-provenance-kq-gauge-goldenrule-v1',
 'whole_source_paper':'PARTIAL',
 'real_material_verified':False,
 'authentic_real_epw_phonon_matrix_found':False,
 'synthetic_normalized_archive_status':original_status,
 'fake_hashed_candidate_status':candidate_status,
 'fake_hashed_candidate_real_verified':candidate_auth,
 'raw_source_tamper_checksum_rejected':tamper_rejected,
 'nk':v['nk'],'nq':v['nq'],'nband':v['nband'],'nmode':v['nmode'],
 'inverse_q_mode_hermiticity_residual_max':paired_res,
 'bare_finite_gaussian_goldenrule_rate_s_inverse_max':float(np.max(result['bare'])),
 'pauli_blocking_rate_s_inverse_difference_max':float(np.max(result['bare']-result['blocked'])),
 'matched_independent_k_gauge_goldenrule_rate_s_inverse_residual_max':float(np.max(abs(result['bare']-rot_rate['bare']))),
 'wrong_H_only_gauge_goldenrule_rate_s_inverse_difference_max':float(np.max(abs(result['bare']-wrong['bare']))),
 'no_exact_kms_from_offshell_gaussian_broadening':not result['finite_broadening_exact_KMS'],
 'physical_spin_provenance_authenticated':False,
 'real_material_kme_si_claim':False,
 'raw_binary_epw_reader_claim':False
}
assert report['synthetic_normalized_archive_status']=='SYNTHETIC_FIXTURE_ONLY'
assert "AUTHENTICITY_UNVERIFIED" in report['fake_hashed_candidate_status']
assert not report['fake_hashed_candidate_real_verified'] and report['raw_source_tamper_checksum_rejected']
assert report['inverse_q_mode_hermiticity_residual_max']<1e-12
assert report['matched_independent_k_gauge_goldenrule_rate_s_inverse_residual_max'] / max(report['bare_finite_gaussian_goldenrule_rate_s_inverse_max'],1)<2e-12
assert report['wrong_H_only_gauge_goldenrule_rate_s_inverse_difference_max'] / report['bare_finite_gaussian_goldenrule_rate_s_inverse_max']>1e-3
assert report['pauli_blocking_rate_s_inverse_difference_max']>1e9
assert report['no_exact_kms_from_offshell_gaussian_broadening']
frozen=ROOT/'results/frozen/round32_epw_material_scattering_contract.json'
expected=json.loads(frozen.read_text())
for key,value in report.items():
    if isinstance(value,float):
        assert abs(value-expected[key]) < 1e-9*max(1,abs(expected[key])),(key,value,expected[key])
    else:
        assert value==expected[key],(key,value,expected[key])
print(json.dumps(report,indent=2))
