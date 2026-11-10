#!/usr/bin/env python3
# ruff: noqa: E402
"""Frozen same-orientation eight-state KM check; no parameter scan."""

from __future__ import annotations

import hashlib
import json
import os
import signal
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

for variable in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[variable] = "1"

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = next(
    parent
    for parent in HERE.parents
    if (parent / "studies/bsf-nmee-connection/bsf_nmee").is_dir()
)
STUDY = REPO / "studies/bsf-nmee-connection"
sys.path[:0] = [str(STUDY), str(REPO / "papers/jia-2026-nonlin-magnetoelectric-edelstein")]

from bsf_nmee.bilayer_audit import AuditBilayer, AuditSpec, integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift
from bsf_nmee.sliding_tight_binding import P_BOTTOM, P_TOP
from nmee_lab.clean_core.core import MU_B


GRID = (61, 60, 96)
B_STEP_T = 0.03125
DU_MEV = 0.05
MASS_MEV = 45.0
T_VALUES_MEV = (0.0, 10.0)
# Partner registries are centered fractions delta=lambda-0.5 = -/+0.25.
REGISTRIES = (
    {"label": "centered_minus_0.25", "centered_fraction": -0.25, "ab_ba_lambda": 0.25},
    {"label": "centered_plus_0.25", "centered_fraction": 0.25, "ab_ba_lambda": 0.75},
)
DIRECTION = np.array([0.5, np.sqrt(3.0) / 2.0])
W = (P_TOP - P_BOTTOM) / 2.0
TIME_CAP_SECONDS = 175
BATCH = 512

# Predeclared numerical-equality tolerances, not material uncertainty bars.
C_ZERO_TOL = 1.0e-10
G_ZERO_TOL = 1.0e-9
C_EQUALITY_TOL = 1.0e-10
C_RANK_TOL = 1.0e-10
G_RANK_TOL = 1.0e-9
SUM_TOL = 1.0e-10
PARTNER_C_TOL = 1.0e-9
PARTNER_G_TOL = 1.0e-8


class TimeCapExceeded(RuntimeError):
    """Raised when the frozen local run exceeds its wall-time cap."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def adjoint(array: np.ndarray) -> np.ndarray:
    return array.conj().swapaxes(-1, -2)


def signed(value: float, tolerance: float) -> int:
    if value > tolerance:
        return 1
    if value < -tolerance:
        return -1
    return 0


def rank_label(lower: float, upper: float, tolerance: float) -> str:
    difference = abs(upper) - abs(lower)
    if difference > tolerance:
        return "upper"
    if difference < -tolerance:
        return "lower"
    return "unresolved"


def response_kernel(
    hamiltonian: np.ndarray,
    velocity: np.ndarray,
    spins: list[np.ndarray],
) -> tuple[np.ndarray, float, float]:
    energies, states = np.linalg.eigh(hamiltonian)
    occupied = states[:, :, :4]
    empty = states[:, :, 4:]
    velocity_eo = adjoint(empty) @ velocity @ occupied
    differences = energies[:, :4, None] - energies[:, None, 4:]
    values = np.array(
        [
            (
                -2.0
                * np.imag((adjoint(occupied) @ spin @ empty) * velocity_eo.swapaxes(1, 2))
                / differences**2
            ).sum(axis=(1, 2))
            for spin in spins
        ]
    )
    return values, float(energies[:, 3].max()), float(energies[:, 4].min())


def blank_accumulator() -> dict[str, Any]:
    return {"C": np.zeros(3), "valence_max": -np.inf, "conduction_min": np.inf}


def set_time_cap() -> None:
    def handler(_signum: int, _frame: Any) -> None:
        raise TimeCapExceeded(f"wall-time cap of {TIME_CAP_SECONDS} seconds exceeded")

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(TIME_CAP_SECONDS)


def main() -> None:
    start = time.monotonic()
    set_time_cap()
    result_path = HERE / "result.json"
    acceptance_path = HERE / "acceptance.json"
    result: dict[str, Any] = {
        "complete": False,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scope": "Frozen KM_DIRECT_1 only: unchanged eight-state direct stacking, t=0/10 meV, two layer-exchange partner registries, no tuning or optional stages.",
        "registry_clarification_before_result_acquisition": {
            "recorded": True,
            "wording": "The frozen +/-0.25 fractions are centered coordinates delta=lambda-0.5. Existing ab_ba_shift therefore receives lambda=0.25 and 0.75; literal lambda=-0.25 and 0.25 are not layer-exchange partners.",
            "scientific_parameters_changed": False,
        },
        "parameters": {
            "stacking": "direct",
            "mass_meV": MASS_MEV,
            "t_meV": list(T_VALUES_MEV),
            "U0_meV": 0.0,
            "du_meV": DU_MEV,
            "B_step_T": B_STEP_T,
            "grid": list(GRID),
            "direction_xy": DIRECTION.tolist(),
            "registries": list(REGISTRIES),
            "time_cap_seconds": TIME_CAP_SECONDS,
            "batch": BATCH,
            "blas_threads": 1,
        },
        "units": {"C": "(nm V T)^-1", "G": "(nm V T eV)^-1", "gap": "eV"},
        "tolerances": {
            "C_zero": C_ZERO_TOL,
            "G_zero": G_ZERO_TOL,
            "C_layer_equality": C_EQUALITY_TOL,
            "C_rank": C_RANK_TOL,
            "G_rank": G_RANK_TOL,
            "layer_sum": SUM_TOL,
            "partner_C": PARTNER_C_TOL,
            "partner_G": PARTNER_G_TOL,
            "classification": "numerical equality only; not a physical uncertainty estimate",
        },
        "provenance": {
            "hypothesis_queue_frozen_file": "hypothesis_queue_frozen.json",
            "hypothesis_queue_sha256": sha256(HERE / "hypothesis_queue_frozen.json"),
            "bilayer_audit_sha256": sha256(STUDY / "bsf_nmee/bilayer_audit.py"),
            "bilayer_protocol_sha256": sha256(STUDY / "bsf_nmee/bilayer_protocol.py"),
        },
        "rows": [],
    }
    acceptance: dict[str, Any] = {
        "complete": False,
        "status": "running",
        "test": "KM_DIRECT_1",
        "primary_rule": "pass only if exactly one t=10 layer reverses relative to the same-sign t=0 direct baseline at a fixed registry; fail if resolved and absent; unresolved is inconclusive",
        "reporting_clarification": "After acquisition, material sign accounting was added: 7268 and 7269 bilayers share lower-negative/upper-positive signs, so their net-sign difference is compensation. This did not alter this model test or its acceptance rule.",
    }
    try:
        template = AuditBilayer(
            AuditSpec(stacking="direct", staggered_eV=MASS_MEV / 1000.0)
        )
        kpoints, weights, grid_meta = integration_grid(template, *GRID)
        result["nk"] = len(kpoints)
        result["grid_meta"] = grid_meta
        models = {
            registry["label"]: AuditBilayer(
                AuditSpec(
                    shift=ab_ba_shift(registry["ab_ba_lambda"]),
                    stacking="direct",
                    staggered_eV=MASS_MEV / 1000.0,
                    tunneling_eV=0.010,
                )
            )
            for registry in REGISTRIES
        }
        accumulators = {
            (registry["label"], t_mev, u_index): blank_accumulator()
            for registry in REGISTRIES
            for t_mev in T_VALUES_MEV
            for u_index in (-1, 0, 1)
        }
        for start_index in range(0, len(kpoints), BATCH):
            chunk = slice(start_index, start_index + BATCH)
            ks = kpoints[chunk]
            ws = weights[chunk]
            for registry in REGISTRIES:
                label = registry["label"]
                model = models[label]
                h10 = model.matrices(ks)
                v10 = sum(
                    DIRECTION[axis] * model.matrices(ks, axis=axis)
                    for axis in range(2)
                )
                spins = [
                    sum(
                        DIRECTION[axis] * model.spins[f"{layer}_{component}"]
                        for axis, component in enumerate(("x", "y"))
                    )
                    for layer in ("bottom", "top", "total")
                ]
                for t_mev in T_VALUES_MEV:
                    h = h10.copy()
                    velocity = v10.copy()
                    if t_mev == 0.0:
                        for matrix in (h, velocity):
                            matrix[:, :4, 4:] = 0.0
                            matrix[:, 4:, :4] = 0.0
                    for u_index in (-1, 0, 1):
                        shifted = h + u_index * DU_MEV / 1000.0 * W
                        target = accumulators[(label, t_mev, u_index)]
                        for field_sign in (1, -1):
                            values, valence_max, conduction_min = response_kernel(
                                shifted
                                - field_sign
                                * B_STEP_T
                                * model.mono.g_factor
                                * MU_B
                                * model.spins["total_z"],
                                velocity,
                                spins,
                            )
                            target["C"] += (
                                field_sign * (values @ ws) / (2.0 * B_STEP_T)
                            )
                            target["valence_max"] = max(
                                target["valence_max"], valence_max
                            )
                            target["conduction_min"] = min(
                                target["conduction_min"], conduction_min
                            )
        rows: list[dict[str, Any]] = []
        for registry in REGISTRIES:
            label = registry["label"]
            for t_mev in T_VALUES_MEV:
                records = [accumulators[(label, t_mev, index)] for index in (-1, 0, 1)]
                cminus, c0, cplus = (record["C"] for record in records)
                gate_slope = (cplus - cminus) / (2.0 * DU_MEV / 1000.0)
                common_gap = min(record["conduction_min"] for record in records) - max(
                    record["valence_max"] for record in records
                )
                row = {
                    "registry": label,
                    "centered_fraction": registry["centered_fraction"],
                    "ab_ba_lambda": registry["ab_ba_lambda"],
                    "shift_60degree_basis": list(ab_ba_shift(registry["ab_ba_lambda"])),
                    "t_meV": t_mev,
                    "order": ["lower", "upper", "net"],
                    "Cminus": cminus.tolist(),
                    "C0": c0.tolist(),
                    "Cplus": cplus.tolist(),
                    "G": gate_slope.tolist(),
                    "C_signs": [signed(float(value), C_ZERO_TOL) for value in c0[:2]],
                    "G_signs": [signed(float(value), G_ZERO_TOL) for value in gate_slope[:2]],
                    "static_largest_layer": rank_label(c0[0], c0[1], C_RANK_TOL),
                    "gate_sensitive_layer": rank_label(
                        gate_slope[0], gate_slope[1], G_RANK_TOL
                    ),
                    "common_gap_eV": common_gap,
                    "insulating": bool(common_gap > 0.0),
                    "layer_sum_error_C0": float(c0[0] + c0[1] - c0[2]),
                    "layer_sum_error_G": float(
                        gate_slope[0] + gate_slope[1] - gate_slope[2]
                    ),
                }
                rows.append(row)
        result["rows"] = rows

        row_map = {(row["registry"], row["t_meV"]): row for row in rows}
        isolated_checks = []
        reversal_checks = []
        for registry in REGISTRIES:
            label = registry["label"]
            isolated = row_map[(label, 0.0)]
            coupled = row_map[(label, 10.0)]
            isolated_signs = isolated["C_signs"]
            coupled_signs = coupled["C_signs"]
            isolated_checks.append(
                {
                    "registry": label,
                    "same_resolved_sign": bool(
                        isolated_signs[0] != 0 and isolated_signs[0] == isolated_signs[1]
                    ),
                    "equal_within_tolerance": bool(
                        abs(isolated["C0"][0] - isolated["C0"][1])
                        <= C_EQUALITY_TOL
                    ),
                    "C_lower_upper": isolated["C0"][:2],
                    "signs": isolated_signs,
                }
            )
            changed = [
                isolated_signs[index] != 0
                and coupled_signs[index] != 0
                and isolated_signs[index] != coupled_signs[index]
                for index in range(2)
            ]
            unresolved = any(sign_value == 0 for sign_value in isolated_signs + coupled_signs)
            reversal_checks.append(
                {
                    "registry": label,
                    "isolated_signs": isolated_signs,
                    "coupled_signs": coupled_signs,
                    "layer_reversed": {"lower": changed[0], "upper": changed[1]},
                    "exactly_one_layer_reversed": bool(sum(changed) == 1),
                    "unresolved": unresolved,
                }
            )
        first = row_map[(REGISTRIES[0]["label"], 10.0)]
        second = row_map[(REGISTRIES[1]["label"], 10.0)]
        partner = {
            "C_swap_max_abs": float(
                max(
                    abs(first["C0"][0] - second["C0"][1]),
                    abs(first["C0"][1] - second["C0"][0]),
                )
            ),
            "G_antiswap_max_abs": float(
                max(
                    abs(first["G"][0] + second["G"][1]),
                    abs(first["G"][1] + second["G"][0]),
                )
            ),
        }
        partner["passed"] = bool(
            partner["C_swap_max_abs"] <= PARTNER_C_TOL
            and partner["G_antiswap_max_abs"] <= PARTNER_G_TOL
        )
        numerical_checks = {
            "all_insulating": all(row["insulating"] for row in rows),
            "layer_sums_close": all(
                abs(row["layer_sum_error_C0"]) <= SUM_TOL
                and abs(row["layer_sum_error_G"]) <= SUM_TOL
                for row in rows
            ),
            "same_orientation_t0": all(
                check["same_resolved_sign"] and check["equal_within_tolerance"]
                for check in isolated_checks
            ),
            "partner_relation": partner["passed"],
        }
        numerical_pass = all(numerical_checks.values())
        any_unresolved = any(check["unresolved"] for check in reversal_checks)
        any_one_layer_reversal = any(
            check["exactly_one_layer_reversed"] for check in reversal_checks
        )
        if not numerical_pass or any_unresolved:
            status = "inconclusive"
        elif any_one_layer_reversal:
            status = "passed"
        else:
            status = "failed_as_observed"
        acceptance.update(
            {
                "complete": True,
                "status": status,
                "numerical_checks": numerical_checks,
                "isolated_checks": isolated_checks,
                "partner_check": partner,
                "reversal_checks": reversal_checks,
                "interpretation": {
                    "passed": "The unchanged direct eight-state model produces a one-layer reversal from identical same-sign isolated blocks at a frozen registry.",
                    "failed_as_observed": "The unchanged direct eight-state model retains both isolated-layer signs at t=10 meV; it does not reproduce the provisional 7268 one-layer reversal under the frozen test.",
                    "inconclusive": "A numerical, gap, symmetry, or sign-resolution check prevents scoring the frozen reversal test.",
                }[status],
                "limits": [
                    "Internal U is not a calibrated external gate field.",
                    "The KM sign and magnitude are not fitted to MoSSe 7268.",
                    "One refined quadrature is a bounded model check, not a material uncertainty analysis.",
                    "The partner allocation swap is not a layer-direction reversal: both t=10 layer C signs remain positive.",
                    "The provisional 7268 and archived 7269 bilayers both have lower-negative/upper-positive signs; their different net signs reflect compensation, while the within-7268 isolated-to-bilayer reversal is a separate target.",
                    "No channel removal, displacement extension, parameter scan, or extra model state was run.",
                ],
            }
        )
        result["complete"] = True
        result["seconds"] = time.monotonic() - start
        acceptance["seconds"] = result["seconds"]
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {error}"
        result["seconds"] = time.monotonic() - start
        acceptance.update(
            {
                "complete": False,
                "status": "time_cap_or_error",
                "error": result["error"],
                "seconds": result["seconds"],
            }
        )
    finally:
        signal.alarm(0)
        result_path.write_text(json.dumps(result, indent=2) + "\n")
        acceptance_path.write_text(json.dumps(acceptance, indent=2) + "\n")
    print(json.dumps({"result": str(result_path), "acceptance": acceptance}, indent=2))


if __name__ == "__main__":
    main()
