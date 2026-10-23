"""Round34 native Bloch Pauli .spn.fmt -> W90 spin with independent HR/eig oracle.

Generated official-format fixtures ONLY; no true DFT spin source, no
raw Fortran unformatted .spn and no real-material KME accepted.
"""
from pathlib import Path
import tempfile,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"src"),str(ROOT/"tests")]
from flux_lab.native_spin_chain import (parse_spn_formatted,parse_u,parse_u_dis,
  parse_v,parse_eig,reconstruct_spin,check_against_native_hr)
from test_native_spin_chain import fixture,load,flags
with tempfile.TemporaryDirectory() as td:
    f=fixture(Path(td)/"phase",bloch_phase=True)
    n=fixture(Path(td)/"unphased",bloch_phase=False)
    sp,u,ud,v,eig,out=load(f)
    plain=load(n)[-1]
    gate=check_against_native_hr(out,f["paths"]["hr"],f["paths"]["r"],flags=flags())
    stale=reconstruct_spin(parse_spn_formatted(n["paths"]["spn"]),u,u_dis=ud,eig=eig)
    overlap=np.matmul(ud["matrix"],u["matrix"])
    sol=out["spin_wannier_k"]
    qcheck=np.max(abs(sol-f["Sw"]))
    covariance=np.max(abs(sol-plain["spin_wannier_k"]))
    wrong=np.max(abs(sol-stale["spin_wannier_k"]))
    sr=fixture(Path(td)/"isolated",disentangle=False)
    spin_iso=load(sr)[-1]
    isolated_error=np.max(abs(spin_iso["spin_wannier_k"]-sr["Sw"]))
    with tempfile.TemporaryDirectory() as tmp:
        pass
    # Wrong transpose of column-major U changes both projected spin/H.
    transformed_bad=[]
    for k in range(f["nk"]):
        misread=f["U"][k].T
        T=f["Udis"][k]@misread
        transformed_bad.append(np.stack([T.conj().T@sp["spin_bloch"][k,a]@T for a in range(3)]))
    wrong_u=np.max(abs(sol-np.array(transformed_bad)))
report={
 "schema":"cui-round34-native-formatted-spn-Udis-U-HR-physical-spin-v1",
 "paper_status":"PARTIAL",
 "fixture":"GENERATED_SPINOR_4B_2W_3K_NOT_REAL_SOURCE",
 "spn_grammar_native_ascii":True,
 "u_u_dis_v_column_major_native_ascii":True,
 "disentanglement_gate":True,
 "isolated_band_gate":True,
 "raw_binary_spn_parsed":False,
 "true_source_spn_verified":False,
 "bloch_to_wannier_spin_mathematics_verified":True,
 "actual_material_spin_provenance_authenticated":False,
 "real_material_KME_ready":False,
 "nk":f["nk"],"nbands":f["nb"],"nwann":f["nw"],
 "spn_upper_triangle_storage_roundtrip_max":float(np.max(abs(sp["spin_bloch"]-f["Sb"]))),
 "u_dis_U_product_vs_native_v_max":float(np.max(abs(overlap-v["matrix"]))),
 "parsed_wannier_pauli_vs_independent_projection_max":float(qcheck),
 "bloch_phase_matched_spin_gauge_covariance_max":float(covariance),
 "wrong_stale_bloch_spin_gauge_max":float(wrong),
 "wrong_row_major_u_vs_wannier_spin_max":float(wrong_u),
 "non_disentangled_2_band_spin_max":float(isolated_error),
 "eig_projected_H_vs_independent_hr_fourier_max_eV":float(gate["h_from_eig_vs_hr_max_eV"]),
 "max_physical_projected_spin_norm":float(max(np.linalg.norm(a,2) for a in sol.reshape(-1,2,2))),
 "spn_diagonal_imaginary_raw_max":float(sp["diag_imag_raw_max"])
}
for key in ("spn_upper_triangle_storage_roundtrip_max","u_dis_U_product_vs_native_v_max",
 "parsed_wannier_pauli_vs_independent_projection_max",
 "bloch_phase_matched_spin_gauge_covariance_max","non_disentangled_2_band_spin_max",
 "eig_projected_H_vs_independent_hr_fourier_max_eV"):
    assert report[key]<2e-7,(key,report[key])
for key in ("wrong_stale_bloch_spin_gauge_max","wrong_row_major_u_vs_wannier_spin_max"):
    assert report[key]>.05,(key,report[key])
assert report["max_physical_projected_spin_norm"]<=1+1e-6
assert report["actual_material_spin_provenance_authenticated"] is False
assert report["paper_status"]=="PARTIAL"
frozen=ROOT/"results/frozen/round34_native_spin_chain_audit.json"
if frozen.exists():
    expected=json.loads(frozen.read_text())
    assert expected.keys()==report.keys()
    for k,x in report.items():
        if isinstance(x,float):
            assert abs(x-expected[k])<2e-9,(k,x,expected[k])
        else:
            assert x==expected[k],(k,x,expected[k])
print(json.dumps(report,indent=2))
