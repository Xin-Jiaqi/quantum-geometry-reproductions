"""Round33 public-format native Wannier90 ASCII gate: synthetic fixture ONLY."""
from pathlib import Path
import json
import sys
import tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"tests"))
from flux_lab.wannier90_raw_gate import parse_hr,parse_r,staged_h_r_import
from test_wannier90_raw_gate import fixture,import_pair

with tempfile.TemporaryDirectory() as td:
    f=fixture(td)
    hr=parse_hr(f["hr"]);r=parse_r(f["r"])
    bundle=import_pair(f)
    H,A=bundle.fourier([.217,.13])
    z=np.exp(2j*np.pi*.217)
    analyticH=f["heff"][1]+z*f["heff"][2]+z.conjugate()*f["heff"][0]
    analyticA=f["aeff"][1]+z*f["aeff"][2]+z.conjugate()*f["aeff"][0]
    factor=np.exp(2j*np.pi*(bundle.R[:,:2]@np.array([.17,0.])))
    goodH=bundle.fourier([.17,0.])[0]
    wrongH=np.einsum("r,rij->ij",factor,hr["raw"])
    miss=0
    try:import_pair(f,ws=True)
    except ValueError:miss=1
    blochSpinBlocked=not bundle.physical_spin_available and not bundle.kme_ready
    centersMismatch=float(np.max(abs(bundle.centers-f["centers"])))
    normedHres=float(np.max(abs(bundle.h_r-f["heff"])))
    normedAres=float(np.max(abs(bundle.position_r-f["aeff"])))
    g=fixture(Path(td)/"expanded",deg=(1,1,1),expanded=True)
    expand=import_pair(g,expanded=True,ws=True,has_wsvec=True,transl=True)
    expandDiff=float(np.max(abs(expand.h_r-f["heff"])))
report={
 "schema":"cui-round33-w90-native-hr-r-degeneracy-spin-blocked-v1",
 "original_paper_status":"PARTIAL",
 "input_data_kind":"GENERATED_FORMAT_FIXTURE_NO_REAL_WANNIER90_FILE",
 "raw_hr_parser_implemented":True,
 "raw_r_parser_implemented":True,
 "authentic_material_input_found":False,
 "wannier_physical_spin_provenance_authenticated":False,
 "spn_bloch_to_wannier_gauge_implemented":False,
 "kme_real_material_ready":False,
 "nrpts":int(hr["nrpts"]),"num_wann":int(hr["nwann"]),
 "wigner_seitz_degeneracy":list(map(int,hr["degeneracy"])),
 "native_weighted_h_reconstruction_max":normedHres,
 "native_weighted_r_reconstruction_max":normedAres,
 "analytic_complex_fourier_H_max":float(np.max(abs(H-analyticH))),
 "analytic_complex_fourier_A_max":float(np.max(abs(A-analyticA))),
 "kspace_H_hermiticity_max":float(np.max(abs(H-H.conj().T))),
 "kspace_A_hermiticity_max":float(np.max(abs(A-np.swapaxes(A.conj(),-1,-2)))),
 "orbital_centers_r_onsite_max":centersMismatch,
 "wrong_omitted_ws_degeneracy_H_max":float(np.max(abs(goodH-wrongH))),
 "wsvec_dependent_input_rejected":bool(miss),
 "expanded_ndegen_input_max":expandDiff,
 "physical_spin_kinetic_response_blocked":bool(blochSpinBlocked),
 "synthetic_source_hr_sha256_length":len(hr["source_sha256"]),
 "synthetic_source_r_sha256_length":len(r["source_sha256"])
}
for x in ("native_weighted_h_reconstruction_max","native_weighted_r_reconstruction_max",
 "analytic_complex_fourier_H_max","analytic_complex_fourier_A_max",
 "kspace_H_hermiticity_max","kspace_A_hermiticity_max","orbital_centers_r_onsite_max",
 "expanded_ndegen_input_max"):
    assert report[x]<5e-13,(x,report[x])
assert report["wrong_omitted_ws_degeneracy_H_max"]>.05
assert report["physical_spin_kinetic_response_blocked"]
assert report["wsvec_dependent_input_rejected"]
assert not report["authentic_material_input_found"] and report["original_paper_status"]=="PARTIAL"
ref=ROOT/"results/frozen/round33_native_wannier90_hr_r_audit.json"
target=json.loads(ref.read_text())
for k,value in report.items():
    if isinstance(value,float):
        assert abs(value-target[k])<2e-11,(k,value,target[k])
    else:
        assert value==target[k],(k,value,target[k])
print(json.dumps(report,indent=2))
