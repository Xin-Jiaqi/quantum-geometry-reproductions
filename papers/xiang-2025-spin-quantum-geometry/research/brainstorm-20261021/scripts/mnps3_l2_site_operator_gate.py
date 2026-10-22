"""Read-only four-Mn site/spin operator covariance preflight for MnPS3 L2.

Run:
  python mnps3_l2_site_operator_gate.py --selftest
  python mnps3_l2_site_operator_gate.py path/to/PRIVATE-l2-site-map.json

A source-backed input requires four physical Mn entries for each daughter,
source_index_to_target mapping and SHA256 provenance. No coordinates or
magnetic patterns are guessed from the paper's schematic. The spin-E test
implements NONRELATIVISTIC [E||MzMy] and is not an SOC magnetic symmetry.
The spin-C2x test is a separate AXIAL spinor rotation expected for the
proper spatial product MzMy=C2x, not proof of an SOC symmetry.

This tests ideal onsite projector algebra only, not PAW/Wannier vertex quality.
No VASP, SSH, scheduler calls, remote paths, or secret output.
"""
import argparse
import json
import re
import numpy as np

S0=np.eye(2, dtype=complex)
SX=np.array([[0,1],[1,0]], complex)
SY=np.array([[0,-1j],[1j,0]], complex)
SZ=np.diag([1,-1]).astype(complex)
SPIN={"x":SX/2,"y":SY/2,"z":SZ/2}


def _validate(data, require_source=True):
    for key in ("sites_plus", "sites_minus", "source_index_to_target"):
        if key not in data: raise ValueError("Missing "+key)
    plus=data["sites_plus"];minus=data["sites_minus"]
    mapping=data["source_index_to_target"]
    if len(plus)!=4 or len(minus)!=4 or sorted(mapping)!=list(range(4)):
        raise ValueError("Exactly four Mn sites per state and bijective index mapping required")
    for label,sites in (("plus",plus),("minus",minus)):
        ids=[]
        for site in sites:
            if not all(k in site for k in ("id","layer","sublattice","eta")):
                raise ValueError(label+" site lacks physical labels")
            if site["layer"] not in ("T","B") or site["eta"] not in (-1,1):
                raise ValueError(label+" invalid layer or moment sign")
            ids.append(site["id"])
        if len(set(ids))!=4 or sorted(s["layer"] for s in sites)!=["B","B","T","T"]:
            raise ValueError(label+" duplicated IDs or non-two-plus-two layer counts")
        for layer in ("T","B"):
            local=[site for site in sites if site["layer"]==layer]
            if sorted(s["sublattice"] for s in local)!=["A","B"]:
                raise ValueError(label+" requires A/B Mn labels in each honeycomb layer")
            if sum(s["eta"] for s in local)!=0:
                raise ValueError(label+" requires the source-claimed intralayer Neel ordering")
    for i,j in enumerate(mapping):
        a,b=plus[i],minus[j]
        if a["layer"] == b["layer"]:
            raise ValueError("Eq13 includes Mz; mapped Mn layer must exchange")
        if a.get("element","Mn")!="Mn" or b.get("element","Mn")!="Mn":
            raise ValueError("Only source-indexed Mn site records accepted")
    if require_source:
        for key in ("source_doi","poscar_plus_sha256","poscar_minus_sha256"):
            if not data.get(key): raise ValueError("Missing provenance: "+key)
        for key in ("poscar_plus_sha256","poscar_minus_sha256"):
            if not re.fullmatch(r"[0-9a-fA-F]{64}",data[key]):
                raise ValueError("Full 64-hex SHA256 required for "+key)
        if data["source_doi"]!="10.1038/s41467-026-69635-2":
            raise ValueError("This gate is restricted to the cited MnPS3 L2 paper")


def _operator(weights,axis):
    return np.kron(np.diag(weights),SPIN[axis])


def _rotation_results(data, rotation_name):
    plus,minus=data["sites_plus"],data["sites_minus"]
    perm=np.zeros((4,4),complex)
    for i,j in enumerate(data["source_index_to_target"]): perm[j,i]=1
    rot=S0 if rotation_name=="spin_E_nonrelativistic" else -1j*SX
    U=np.kron(perm,rot)
    weights={
        "layer_plus":np.array([1 if s["layer"]=="T" else -1 for s in plus]),
        "layer_minus":np.array([1 if s["layer"]=="T" else -1 for s in minus]),
        "neel_plus":np.array([s["eta"] for s in plus]),
        "neel_minus":np.array([s["eta"] for s in minus])
    }
    out={}
    for obj in ("layer","neel"):
        for axis in ("x","y","z"):
            O=_operator(weights[obj+"_plus"],axis)
            P=_operator(weights[obj+"_minus"],axis)
            transported=U@O@U.conj().T
            unit_res={s:float(np.linalg.norm(transported-s*P)) for s in (1,-1)}
            parity=min(unit_res,key=unit_res.get)
            passed=unit_res[parity] < 1e-10
            out[obj+"_"+axis]={"parity":int(parity) if passed else None,
                                 "residual":unit_res[parity],
                                 "status":"ALGEBRA_PASS" if passed else "NOT_SIMPLE_PARITY"}
    return out


def audit(data, require_source=True):
    _validate(data,require_source=require_source)
    plus=data["sites_plus"]
    layer=np.array([1 if s["layer"]=="T" else -1 for s in plus])
    neel=np.array([s["eta"] for s in plus])
    return {
       "status":"SOURCE_MAP_ALGEBRA_CHECKED" if require_source else "SYNTHETIC_SELFTEST_ONLY",
       "source_provenance_verified_by_code":False,
       "warning":"Validating labels and hashes is not independently verifying their source or physical SOC Hamiltonian",
       "layer_neel_weight_overlap":int(np.dot(layer,neel)),
       "spin_E_nonrelativistic":_rotation_results(data,"spin_E_nonrelativistic"),
       "spin_C2x_axial_SOC_candidate":_rotation_results(data,"spin_C2x_axial_SOC_candidate"),
       "actual_SOC_magnetic_group":"UNDETERMINED (needs genuine spin+orbital symmetry of both states)",
       "physical_layer_spin_exchange_matrices":"NOT_AVAILABLE (need accepted same-gauge vertex exports)",
    }


def selftest():
    eta=[1,-1,-1,1]
    plus=[dict(id="synthetic_plus_"+str(i),layer=l,sublattice=s,eta=eta[i])
          for i,(l,s) in enumerate((("T","A"),("T","B"),("B","A"),("B","B")))]
    mapping=[2,3,0,1]
    minus=[dict(id="synthetic_minus_"+str(i),layer=l,sublattice=s,eta=eta[i])
           for i,(l,s) in enumerate((("T","A"),("T","B"),("B","A"),("B","B")))]
    d=dict(sites_plus=plus,sites_minus=minus,source_index_to_target=mapping)
    r=audit(d,require_source=False)
    assert r["layer_neel_weight_overlap"]==0
    for key in ("spin_E_nonrelativistic","spin_C2x_axial_SOC_candidate"):
        assert r[key]["layer_x"]["parity"] == -1
    assert r["spin_E_nonrelativistic"]["layer_y"]["parity"] == -1
    assert r["spin_C2x_axial_SOC_candidate"]["layer_y"]["parity"] == 1
    assert r["actual_SOC_magnetic_group"].startswith("UNDETERMINED")
    for invalid in (dict(d,source_index_to_target=[2,3,0,0]),dict(d,source_index_to_target=[0,1,2,3])):
        try:
            audit(invalid,require_source=False)
        except ValueError: pass
        else: raise AssertionError("Invalid mapping accepted")
    print("SELFTEST PASS: bijection, layer swap, orthogonality, spin-E != axial-C2x, rejects broken inputs")
    print(json.dumps(r,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("input",nargs="?",help="Private source-backed four-Mn site-map JSON")
    parser.add_argument("--selftest",action="store_true")
    args=parser.parse_args()
    if args.selftest: selftest()
    elif args.input:
        with open(args.input,encoding="utf8") as f: data=json.load(f)
        print(json.dumps(audit(data),indent=2))
    else: parser.error("Provide a source-backed private JSON or --selftest")
