"""Round26 independent synthetic operator, Wilson and negative-control audit."""
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"tests"))
from flux_lab.wannier_operator_contract import synthetic_soc_bundle,global_gauge_transform
from test_wannier_operator_contract import projector_qgt,wilson_oracle

q=np.array([.183,-.271])
b=synthetic_soc_bundle()
d=b.fourier(q)
J,DS=wilson_oracle(b,q)
rng=np.random.default_rng(26)
U,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
rot=global_gauge_transform(b,U).fourier(q)
g=projector_qgt(d["H"],d["covariant_H"],d["spin"])
g_rot=projector_qgt(rot["H"],rot["covariant_H"],rot["spin"])
# Use the SAME seeded rotation as the positive joint-gauge control;
# this makes the deliberately unmatched-spin negative control reproducible.
g_bad=projector_qgt(rot["H"],rot["covariant_H"],d["spin"])
g_correct=g_rot
torque=np.stack([d["H"]@s-s@d["H"] for s in d["spin"]])
report={
    "schema":"cui-round26-same-wannier-operator-contract-v1",
    "paper_status":"PARTIAL",
    "fractional_k":q.tolist(),
    "model":"synthetic 4-band spin x orbital SOC; no material-input authentication",
    "units":{"energy":"eV","length":"angstrom","spin":"dimensionless_pauli"},
    "wilson_finite_difference_dH_residual_max":float(np.max(abs(J-d["covariant_H"]))),
    "wilson_finite_difference_dS_residual_max":float(np.max(abs(DS-d["covariant_spin"]))),
    "global_unitary_QGT_residual_max":float(np.max(abs(g-g_rot))),
    "wrong_spin_gauge_QGT_error_max":float(np.max(abs(g_bad-g_correct))),
    "missing_connection_error_max":float(np.max(abs(d["partial_H"]-d["covariant_H"]))),
    "covariant_vertex_abs_max":float(np.max(abs(d["covariant_H"]))),
    "spin_torque_commutator_abs_max":float(np.max(abs(torque))),
    "actual_material_spin_verified":False,"finite_frequency_metallic_KME":False}
assert report["wilson_finite_difference_dH_residual_max"]<5e-9
assert report["wilson_finite_difference_dS_residual_max"]<5e-9
assert report["global_unitary_QGT_residual_max"]<3e-12
assert report["wrong_spin_gauge_QGT_error_max"]>.02
assert report["missing_connection_error_max"]>.005
assert report["spin_torque_commutator_abs_max"]>.1
frozen=json.loads((ROOT/"results/frozen/round26_same_gauge_operator_audit.json").read_text())
assert frozen["schema"]==report["schema"] and frozen["paper_status"]=="PARTIAL"
for key in ("covariant_vertex_abs_max","missing_connection_error_max","spin_torque_commutator_abs_max","wrong_spin_gauge_QGT_error_max"):
    assert abs(report[key]-frozen[key])<1e-11, key
print(json.dumps(report,indent=2))
