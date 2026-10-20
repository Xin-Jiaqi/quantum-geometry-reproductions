"""Stage-6 AgF2 DFT receipt auditor. Does NOT execute DFT or certify physics.

No external dependencies. Only verifies supplied assertions/metadata consistency,
optional local POSCAR hash, and a conditional domain-sign diagnostic. A real
physics validation requires inspection of raw run outputs and magnetic operations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

STATES = ("B1", "B2")
SCHEMA = "bsf_igmc_agf2_stage6_v1"


def _sha256(s: Any) -> bool:
    return isinstance(s, str) and len(s) == 64 and all(c in "0123456789abcdef" for c in s.lower())


def _finite(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _present(value: Any) -> bool:
    return value is not None and value != "" and value != "UNSET"


def parity_diagnostic(a: float, b: float, expected_sign: int, ea: float, eb: float) -> dict[str, Any]:
    """Compare x(B2)=sign*x(B1); do not mistake a noise-level zero for switching."""
    if not all(_finite(x) for x in (a, b, ea, eb)) or ea < 0 or eb < 0:
        raise ValueError("IGMC values and nonnegative absolute errors must be finite")
    sigma = math.hypot(ea, eb)
    scale = max(abs(a), abs(b))
    if scale <= 5 * max(ea, eb, 1e-14):
        return {"status": "UNRESOLVED_ZERO_SIGNAL", "expected_sign": expected_sign}
    residual = b - expected_sign * a
    rel = abs(residual) / scale
    return {"status": "CONSISTENT_NUMERICALLY" if rel <= 0.05 and (sigma == 0 or abs(residual) <= 5 * sigma) else "REQUIRES_PHYSICS_REVIEW",
            "expected_sign": expected_sign, "relative_residual": rel,
            "absolute_residual": residual, "combined_error": sigma}


def audit(payload: dict[str, Any], *, root: Path | None = None) -> dict[str, Any]:
    issues: list[str] = []
    notes: list[str] = []
    if payload.get("schema_version") != SCHEMA:
        issues.append("SCHEMA_NOT_MATCHED")
    camp = payload.get("campaign") or {}
    if camp.get("material") != "AgF2":
        issues.append("MATERIAL_NOT_AGF2")
    branch = camp.get("magnetic_branch")
    strain = camp.get("strain_percent")
    if branch not in ("type_I", "type_II"):
        issues.append("MAGNETIC_BRANCH_NOT_FROZEN")
    if not _finite(strain) or (branch == "type_I" and strain != 0) or (branch == "type_II" and strain != 4):
        issues.append("BRANCH_STRAIN_MISMATCH")
    if camp.get("coordinate_system") != "same Cartesian laboratory axes x=[100],y=[010],z=[001]":
        issues.append("LAB_AXES_NOT_FROZEN")
    if camp.get("soc") is not True:
        issues.append("SOC_NOT_ENABLED_IN_PROTOCOL")

    for k in ("protocol_sha256", "source_structure_sha256", "potcar_lock_sha256"):
        if not _sha256(camp.get(k)):
            issues.append("SOURCE_OR_PROTOCOL_LOCK_MISSING_" + k)
    for k in ("source_structure_url", "method_source_url", "functional", "vdw_method"):
        if not _present(camp.get(k)):
            issues.append("SOURCE_OR_PROTOCOL_MISSING_" + k)
    if not _finite(camp.get("ag_hubbard_u_ev")) or camp["ag_hubbard_u_ev"] < 0:
        issues.append("HUBBARD_U_NOT_EXPLICIT")
    if camp.get("physics_approved_by_user") is not True:
        issues.append("SCIENCE_INPUT_NOT_USER_APPROVED")

    states = payload.get("states") or {}
    state_results: dict[str, dict[str, Any]] = {}
    for name in STATES:
        state = states.get(name) or {}
        problems: list[str] = []
        for field in ("poscar_sha256", "relaxed_outcar_sha256"):
            if not _sha256(state.get(field)):
                problems.append("MISSING_" + field.upper())
        if not _present(state.get("poscar_path")):
            problems.append("POSCAR_PATH_NOT_SUPPLIED")
        elif root is not None:
            relative = Path(str(state["poscar_path"]))
            if relative.is_absolute() or ".." in relative.parts:
                problems.append("POSCAR_PATH_UNSAFE")
            else:
                file = (root / relative)
                if not file.is_file():
                    problems.append("POSCAR_FILE_MISSING")
                elif _sha256(state.get("poscar_sha256")) and hashlib.sha256(file.read_bytes()).hexdigest() != state["poscar_sha256"]:
                    problems.append("POSCAR_HASH_MISMATCH")
        if not isinstance(state.get("num_fu"), int) or isinstance(state.get("num_fu"), bool) or state["num_fu"] <= 0:
            problems.append("FU_COUNT_NOT_VERIFIED")
        for field in ("total_energy_ev", "max_force_ev_per_ang"):
            if not _finite(state.get(field)):
                problems.append("MISSING_" + field.upper())
        if _finite(state.get("max_force_ev_per_ang")) and state["max_force_ev_per_ang"] < 0:
            problems.append("NEGATIVE_MAX_FORCE")
        if state.get("scf_converged") is not True or state.get("ionic_converged") is not True:
            problems.append("DFT_CONVERGENCE_NOT_CONFIRMED")
        if state.get("magnetic_group_verified") is not True:
            problems.append("MAGNETIC_GROUP_UNVERIFIED")
        if state.get("magnetic_axis") not in ("x", "y", "z"):
            problems.append("MAGNETIC_AXIS_UNVERIFIED")
        if not _finite(state.get("direct_gap_ev")) or state["direct_gap_ev"] <= 0:
            problems.append("INSULATING_DIRECT_GAP_NOT_VERIFIED")
        if not _finite(state.get("indirect_gap_ev")) or state["indirect_gap_ev"] <= 0:
            problems.append("INSULATING_INDIRECT_GAP_NOT_VERIFIED")
        state_results[name] = {"missing_or_unverified": problems}

    evidence = payload.get("calculation_evidence") or {}
    for item in ("software_and_version", "kmesh_and_convergence_record", "energy_cutoff_convergence_record",
                 "vacuum_convergence_record", "magnetic_order_comparison_record"):
        if not _present(evidence.get(item)):
            issues.append("DFT_EVIDENCE_MISSING_" + item)
    if "SOURCE" in " ".join(issues) or "SCIENCE_INPUT_NOT_USER_APPROVED" in issues:
        status = "BLOCKED_SOURCE_PROTOCOL"
    elif issues or any(v["missing_or_unverified"] for v in state_results.values()):
        status = "BLOCKED_DFT_MAGNETIC"
    else:
        status = "BLOCKED_OPERATOR_AUDIT"

    if all(isinstance((states.get(s) or {}).get("num_fu"), int) and not isinstance((states.get(s) or {}).get("num_fu"), bool) for s in STATES):
        if states["B1"]["num_fu"] != states["B2"]["num_fu"]:
            issues.append("DOMAIN_FORMULA_UNIT_COUNTS_DIFFER")
    energies = [(states.get(s) or {}).get("total_energy_ev") for s in STATES]
    nus = [(states.get(s) or {}).get("num_fu") for s in STATES]
    delta_mev_fu = None
    if all(_finite(x) for x in energies) and all(isinstance(x, int) and not isinstance(x, bool) and x>0 for x in nus):
        delta_mev_fu = 1000 * (energies[1]/nus[1] - energies[0]/nus[0])
        notes.append("B2-B1 total energy normalized to formula units; NOT the inter-magnetic-order Type-I/II splitting.")
    parity = {}
    if status == "BLOCKED_OPERATOR_AUDIT":
        if all((states.get(s) or {}).get("operators_gauge_validated") is True for s in STATES) and _present(evidence.get("wannier_operator_audit_record")):
            status = "AWAITING_IGMC_RESPONSE"
            freq = evidence.get("field_frequency_thz")
            g = evidence.get("zeeman_g_factor")
            if _finite(freq) and freq>0 and _finite(g) and g>0:
                vals = [states[s].get("igmc") for s in STATES]
                if all(isinstance(v, dict) for v in vals):
                    expected = -1 if branch == "type_I" else 1
                    for component in ("D_xx", "D_yy"):
                        a = vals[0].get(component); b = vals[1].get(component)
                        ea = vals[0].get(component + "_abs_error"); eb = vals[1].get(component + "_abs_error")
                        units_ok = (vals[0].get("units") == vals[1].get("units") == "sheet_alpha_P_per_B")
                        if not units_ok:
                            parity[component] = {"status":"UNITS_NOT_SHEET_ALPHA_P_PER_B"}
                        elif all(_finite(t) for t in (a,b,ea,eb)):
                            parity[component] = parity_diagnostic(a,b,expected,ea,eb)
                        else:
                            parity[component] = {"status":"MISSING_RESPONSE_OR_ERROR_ESTIMATE"}
                    if all(v.get("status") == "CONSISTENT_NUMERICALLY" for v in parity.values()):
                        status = "REVIEW_REQUIRED_NUMERICALLY_CONSISTENT"
                    else:
                        status = "REVIEW_REQUIRED_RESPONSE_PARTIAL"
            else:
                notes.append("Frequency and Zeeman g must be recorded before judging finite-frequency response.")

    if issues and status == "BLOCKED_OPERATOR_AUDIT":
        status = "BLOCKED_DFT_MAGNETIC"
    if payload.get("data_origin") != "real_dft":
        notes.append("No real DFT provenance asserted. Never count this as a material calculation.")
        if status.startswith("REVIEW_REQUIRED"):
            status = "BLOCKED_REAL_DFT_PROVENANCE"
    return {"status":status, "data_origin":payload.get("data_origin"),
            "issues":issues, "state_checks":state_results, "delta_E_B2_minus_B1_meV_per_fu":delta_mev_fu,
            "symmetry_parity_diagnostic":parity,
            "notes":notes + ["No response result here is independently authenticated: scientific review of real OUTCAR, magnetic symmetry and Wannier gauge is mandatory."]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receipt", type=Path, required=True)
    ap.add_argument("--root", type=Path, help="Optional PRIVATE root for hashing relative POSCAR files")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    report = audit(json.loads(args.receipt.read_text(encoding="utf-8")), root=args.root)
    message = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(message, encoding="utf-8")
    print(message, end="")


if __name__ == "__main__":
    main()
