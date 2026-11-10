#!/usr/bin/env python3
"""Bounded local synthesis of layer direction, gate control, and cancellation.

This script reads only small archived CSV/JSON records.  It does not run an
electronic-structure calculation, access a remote machine, or fit a new model.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent


def find_repo_root(path: Path) -> Path:
    for parent in (path, *path.parents):
        if (parent / ".git").exists():
            return parent
    raise RuntimeError("repository root not found")


ROOT = find_repo_root(HERE)
SOURCES = {
    "layer_controls": ROOT
    / "studies/bsf-nmee-connection/results/materials/layer-channel-sensitivity-20261110/all_controls.csv",
    "layer_representative": ROOT
    / "studies/bsf-nmee-connection/results/materials/layer-channel-sensitivity-20261110/representative.csv",
    "criterion_controls": ROOT
    / "studies/bsf-nmee-connection/results/materials/material-criterion-20261110/all_controls.csv",
    "criterion_tests": ROOT
    / "studies/bsf-nmee-connection/results/materials/material-criterion-20261110/prediction_tests.csv",
    "criterion_checks": ROOT
    / "studies/bsf-nmee-connection/results/materials/material-criterion-20261110/final_checks.json",
    "km_controls": ROOT
    / "studies/bsf-nmee-connection/results/mechanism/km-channel-control-20261110/all_results.csv",
    "inse_mirror": ROOT
    / "studies/bsf-nmee-connection/results/materials/inse-mirror-gate-20261109/source_data.csv",
    "mosse_7266": Path(
        "/Users/jiaqi/Desktop/worksofJiaqi/project-x/outputs/MoSSe-NMEE/7266-blind/RESULTS.md"
    ),
    "prediction_freeze": HERE / "prediction_freeze.json",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = []
        for line, row in enumerate(csv.DictReader(handle), start=2):
            row["_source_line"] = str(line)
            rows.append(row)
        return rows


def f(value: object) -> float:
    if isinstance(value, (str, int, float)):
        return float(value)
    raise TypeError(f"expected a numeric scalar, got {type(value).__name__}")


def sign_label(value: float, tolerance: float = 0.0) -> str:
    if value > tolerance:
        return "positive"
    if value < -tolerance:
        return "negative"
    return "zero_or_unresolved"


def amplitude_trend(c_value: float, slope: float) -> str:
    product = c_value * slope
    if product > 0:
        return "strengthens"
    if product < 0:
        return "weakens"
    return "zero_or_unresolved"


def dominant(lower: float, upper: float) -> str:
    al, au = abs(lower), abs(upper)
    if al > au:
        return "lower"
    if au > al:
        return "upper"
    return "tie"


def metrics(c_lower: float, c_upper: float, g_lower: float, g_upper: float) -> dict[str, object]:
    c_sum = abs(c_lower) + abs(c_upper)
    g_sum = abs(g_lower) + abs(g_upper)
    c_net = c_lower + c_upper
    g_net = g_lower + g_upper
    c_retention = abs(c_net) / c_sum if c_sum else None
    g_retention = abs(g_net) / g_sum if g_sum else None
    static_largest = dominant(c_lower, c_upper)
    gate_sensitive = dominant(g_lower, g_upper)
    return {
        "C_lower": c_lower,
        "C_upper": c_upper,
        "C_net_recomputed": c_net,
        "G_lower": g_lower,
        "G_upper": g_upper,
        "G_net_recomputed": g_net,
        "direction_lower": sign_label(c_lower),
        "direction_upper": sign_label(c_upper),
        "positive_control_amplitude_trend_lower": amplitude_trend(c_lower, g_lower),
        "positive_control_amplitude_trend_upper": amplitude_trend(c_upper, g_upper),
        "static_largest_layer": static_largest,
        "gate_sensitive_layer": gate_sensitive,
        "static_gate_rank_mismatch": static_largest != gate_sensitive,
        "static_lower_absolute_share": abs(c_lower) / c_sum if c_sum else None,
        "gate_lower_absolute_share": abs(g_lower) / g_sum if g_sum else None,
        "static_rank_margin": abs(abs(c_lower) - abs(c_upper)) / c_sum if c_sum else None,
        "gate_rank_margin": abs(abs(g_lower) - abs(g_upper)) / g_sum if g_sum else None,
        "static_net_retention": c_retention,
        "gate_net_retention": g_retention,
        "static_cancellation_fraction": 1 - c_retention if c_retention is not None else None,
        "gate_cancellation_fraction": 1 - g_retention if g_retention is not None else None,
    }


def make_case(
    *,
    case_id: str,
    family: str,
    material: str,
    configuration: str,
    field_definition: str,
    evidence_class: str,
    source_key: str,
    source_line: str,
    numerical_setting: str,
    c_lower: float,
    c_upper: float,
    g_lower: float,
    g_upper: float,
    c_unit: str,
    g_unit: str,
    endpoint_sign_change: str,
    notes: str,
) -> dict[str, object]:
    row: dict[str, object] = {
        "case_id": case_id,
        "family": family,
        "material": material,
        "configuration": configuration,
        "field_definition": field_definition,
        "evidence_class": evidence_class,
        "source_path": relative(SOURCES[source_key]),
        "source_line": source_line,
        "source_sha256": sha256(SOURCES[source_key]),
        "numerical_setting": numerical_setting,
        "C_unit": c_unit,
        "G_unit": g_unit,
        "endpoint_layer_sign_change": endpoint_sign_change,
        "notes": notes,
    }
    row.update(metrics(c_lower, c_upper, g_lower, g_upper))
    return row


cases: list[dict[str, object]] = []

# Physical full-interface response under independently screened external fields.
seen: set[tuple[object, ...]] = set()
for row in read_csv(SOURCES["layer_controls"]):
    if row["channel"] != "full" or f(row["alpha"]) != 1.0:
        continue
    key = (
        row["material"],
        row["nb"],
        row["C0_lower"],
        row["C0_upper"],
        row["G_lower"],
        row["G_upper"],
    )
    if key in seen:
        continue
    seen.add(key)
    material = row["material"]
    configuration = "MoSSe_7269_A" if material == "MoSSe" else "InSe_864_A"
    cases.append(
        make_case(
            case_id=f"{configuration}_external_nb{row['nb']}",
            family="full_material",
            material=material,
            configuration=configuration,
            field_definition="screened_external_pm0.10_V_nm",
            evidence_class="existing_fixed_ion_independent_SCF",
            source_key="layer_controls",
            source_line=row["_source_line"],
            numerical_setting=f"6x6 archive; nb={row['nb']}; full response duplicated across frozen selected-space windows",
            c_lower=f(row["C0_lower"]),
            c_upper=f(row["C0_upper"]),
            g_lower=f(row["G_lower"]),
            g_upper=f(row["G_upper"]),
            c_unit="1e-9 (nm V T)^-1",
            g_unit="1e-9 (V^2 T)^-1",
            endpoint_sign_change=str(
                row["endpoint_sign_change_lower"] == "True"
                or row["endpoint_sign_change_upper"] == "True"
            ),
            notes="Physical full interface; absolute mesh convergence remains exploratory.",
        )
    )

# Physical full interfaces under a common diagnostic internal layer-energy step.
seen.clear()
for row in read_csv(SOURCES["criterion_controls"]):
    if row["source"] != "validate" or row["mode"] != "internal_step" or f(row["alpha"]) != 1.0:
        continue
    key = (
        row["material"],
        row["nb"],
        row["C0_lower"],
        row["C0_upper"],
        row["G_lower"],
        row["G_upper"],
    )
    if key in seen:
        continue
    seen.add(key)
    material = row["material"]
    ids = {"MoSSe": "MoSSe_7269_A", "InSe": "InSe_864_A", "MoS2": "MoS2_1930"}
    configuration = ids[material]
    cases.append(
        make_case(
            case_id=f"{configuration}_internal_nb{row['nb']}",
            family="full_material",
            material=material,
            configuration=configuration,
            field_definition="diagnostic_internal_layer_step_pm0.001_eV",
            evidence_class="existing_zero_gate_cache_hypothetical_internal_step",
            source_key="criterion_controls",
            source_line=row["_source_line"],
            numerical_setting=f"6x6 archive; nb={row['nb']}; alpha=1 physical interface",
            c_lower=f(row["C0_lower"]),
            c_upper=f(row["C0_upper"]),
            g_lower=f(row["G_lower"]),
            g_upper=f(row["G_upper"]),
            c_unit="1e-9 (nm V T)^-1",
            g_unit="1e-9 (nm V T eV)^-1",
            endpoint_sign_change=str(
                row["gate_endpoint_sign_change_lower"] == "True"
                or row["gate_endpoint_sign_change_upper"] == "True"
            ),
            notes="Common diagnostic, not a self-consistent external gate; MoS2 lacks matching symmetric small-gate SCFs.",
        )
    )

# The m001-related InSe B domain, reconstructed from archived -0.10, 0, +0.10 states.
mirror_rows = read_csv(SOURCES["inse_mirror"])
mirror_groups: dict[str, dict[float, dict[str, str]]] = defaultdict(dict)
for row in mirror_rows:
    if row["domain"] == "B" and row["mesh"] == "6" and row["step"] == "standard":
        mirror_groups[row["cutoff"]][f(row["F_V_nm"])] = row
for cutoff, by_field in sorted(mirror_groups.items(), key=lambda item: int(item[0])):
    if not {-0.1, 0.0, 0.1}.issubset(by_field):
        continue
    zero, minus, plus = by_field[0.0], by_field[-0.1], by_field[0.1]
    g_lower = (f(plus["C_lower"]) - f(minus["C_lower"])) / 0.2
    g_upper = (f(plus["C_upper"]) - f(minus["C_upper"])) / 0.2
    endpoint_change = (
        sign_label(f(minus["C_lower"])) != sign_label(f(plus["C_lower"]))
        or sign_label(f(minus["C_upper"])) != sign_label(f(plus["C_upper"]))
    )
    cases.append(
        make_case(
            case_id=f"InSe_864_B_external_cutoff{cutoff}",
            family="full_material",
            material="InSe",
            configuration="InSe_864_B_m001_partner",
            field_definition="screened_external_pm0.10_V_nm",
            evidence_class="existing_fixed_ion_independent_SCF",
            source_key="inse_mirror",
            source_line=f"{minus['_source_line']},{zero['_source_line']},{plus['_source_line']}",
            numerical_setting=f"6x6 archive; retained bands={cutoff}; standard magnetic step",
            c_lower=f(zero["C_lower"]),
            c_upper=f(zero["C_upper"]),
            g_lower=g_lower,
            g_upper=g_upper,
            c_unit="1e-9 (nm V T)^-1",
            g_unit="1e-9 (V^2 T)^-1",
            endpoint_sign_change=str(endpoint_change),
            notes="Symmetry-related domain, not an independent material family; net response is mesh-sensitive.",
        )
    )

# Selected eight-state KM controls: held-out phase cases and fine strong-coupling contrasts.
for row in read_csv(SOURCES["km_controls"]):
    t = f(row["t_meV"])
    u = f(row["U_meV"])
    alpha = f(row["retained_alpha"])
    source = row["source"]
    selected = (
        source == "weak_validation.json" and t == 2 and u == 0 and row["channel"] == "Vb_Ct" and alpha == 0.5
    ) or (
        source == "phase_validation.json" and t == 2 and u == 0 and row["channel"] == "Vb_Ct" and alpha == -0.5
    ) or (
        source == "contrast_validation.json"
        and t in (10, 20)
        and u == -4
        and ((row["channel"] == "full" and alpha == 1) or (row["channel"] == "Vb_Ct" and alpha == 0))
    )
    if not selected:
        continue
    cases.append(
        make_case(
            case_id=f"KM_t{t:g}_U{u:g}_{row['channel']}_alpha{alpha:g}",
            family="eight_state_KM",
            material="KM8",
            configuration="fraction_1_4_mass45meV_angle60deg",
            field_definition="internal_upper_minus_lower_onsite_U",
            evidence_class="bounded_local_model_control",
            source_key="km_controls",
            source_line=row["_source_line"],
            numerical_setting=f"source={source}; nk={row['nk']}; t={t:g}meV; U={u:g}meV; alpha={alpha:g}",
            c_lower=f(row["C_lower"]),
            c_upper=f(row["C_upper"]),
            g_lower=f(row["G_lower"]),
            g_upper=f(row["G_upper"]),
            c_unit="(nm V T)^-1",
            g_unit="(nm V T eV)^-1",
            endpoint_sign_change="not_evaluated_here",
            notes="Diagnostic path scaling; not a proposed material operation or material fit.",
        )
    )


interventions: list[dict[str, object]] = []
for row in read_csv(SOURCES["layer_representative"]):
    record: dict[str, object] = {
        "case_id": f"{row['material']}_{row['channel']}_alpha{row['alpha']}",
        "material": row["material"],
        "channel_control": row["channel"],
        "retained_alpha": f(row["alpha"]),
        "source_path": relative(SOURCES["layer_representative"]),
        "source_line": row["_source_line"],
        "source_sha256": sha256(SOURCES["layer_representative"]),
        "notes": "Hamiltonian intervention on pre-existing screened gate states; not rescreened and not a direct laboratory operation.",
    }
    record.update(
        metrics(
            f(row["C0_lower"]),
            f(row["C0_upper"]),
            f(row["G_lower"]),
            f(row["G_upper"]),
        )
    )
    interventions.append(record)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"no rows for {path}")
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def aggregate(rows: list[dict[str, object]]) -> dict[str, object]:
    def bounds(field: str) -> list[float]:
        values = [f(row[field]) for row in rows]
        return [min(values), max(values)]

    return {
        "numerical_settings": len(rows),
        "static_largest_layers": sorted({str(row["static_largest_layer"]) for row in rows}),
        "gate_sensitive_layers": sorted({str(row["gate_sensitive_layer"]) for row in rows}),
        "rank_mismatch_in_every_setting": all(bool(row["static_gate_rank_mismatch"]) for row in rows),
        "static_lower_absolute_share_range": bounds("static_lower_absolute_share"),
        "gate_lower_absolute_share_range": bounds("gate_lower_absolute_share"),
        "static_rank_margin_range": bounds("static_rank_margin"),
        "gate_rank_margin_range": bounds("gate_rank_margin"),
        "static_net_retention_range": bounds("static_net_retention"),
        "gate_net_retention_range": bounds("gate_net_retention"),
    }


full_cases = [row for row in cases if row["family"] == "full_material"]
groups: dict[str, list[dict[str, object]]] = defaultdict(list)
for full_case in full_cases:
    groups[f"{full_case['configuration']}::{full_case['field_definition']}"].append(full_case)

external_cases = [
    row for row in full_cases if row["field_definition"] == "screened_external_pm0.10_V_nm"
]
internal_cases = [
    row for row in full_cases if row["field_definition"] == "diagnostic_internal_layer_step_pm0.001_eV"
]


def select(case_id: str) -> dict[str, object]:
    return next(row for row in cases if row["case_id"] == case_id)


mosse_external = select("MoSSe_7269_A_external_nb112")
inse_a_external = select("InSe_864_A_external_nb160")
inse_b_external = select("InSe_864_B_external_cutoff160")
km_plus = select("KM_t2_U0_Vb_Ct_alpha0.5")
km_minus = select("KM_t2_U0_Vb_Ct_alpha-0.5")

criterion_tests = read_csv(SOURCES["criterion_tests"])
criterion_boundary: dict[str, dict[str, dict[str, int]]] = {}
for material in ("MoSSe", "MoS2", "InSe"):
    modes = sorted({row["mode"] for row in criterion_tests if row["material"] == material})
    criterion_boundary[material] = {}
    for mode in modes:
        blind = [
            row
            for row in criterion_tests
            if row["material"] == material and row["mode"] == mode and row["blind_new_control"] == "True"
        ]
        full = [
            row
            for row in criterion_tests
            if row["material"] == material and row["mode"] == mode and f(row["alpha"]) == 1.0
        ]
        criterion_boundary[material][mode] = {
            "held_out_half_amplitude_groups": len(blind),
            "held_out_layer_slope_signs_correct": sum(
                row[side] == "True"
                for row in blind
                for side in ("sign_correct_lower", "sign_correct_upper")
            ),
            "held_out_layer_slope_signs_total": 2 * len(blind),
            "physical_endpoint_groups": len(full),
            "physical_endpoint_layer_gain_directions_correct": sum(
                row[side] == "True"
                for row in full
                for side in ("gain_direction_correct_lower", "gain_direction_correct_upper")
            ),
            "physical_endpoint_layer_gain_directions_total": 2 * len(full),
        }

criterion_checks = json.loads(SOURCES["criterion_checks"].read_text(encoding="utf-8"))
mixing = {
    item["material"]: {
        "background_direct_gap_eV": item["background_direct_gap_eV"],
        "background_eta_max": item["background_eta_max"],
        "classification": "post_hoc_applicability_descriptor_not_a_threshold",
    }
    for item in criterion_checks["mixing_diagnostics"]
}

strong_km = [
    row
    for row in cases
    if row["family"] == "eight_state_KM" and "contrast_validation.json" in str(row["numerical_setting"])
]

summary = {
    "scope": {
        "description": "bounded local synthesis of existing eight-state KM and small archived DFT-derived response tables",
        "new_DFT": False,
        "remote_access_or_jobs": False,
        "twelve_state_fit": False,
        "new_material_result": False,
    },
    "observable_separation": {
        "layer_response_direction": "sign(C_l)",
        "positive_control_magnitude_change": "sign(C_l*G_l), provided the layer does not cross zero",
        "absolute_gate_sensitive_layer": "argmax_l |G_l|; this is distinct from relative sensitivity |G_l|/|C_l|",
        "static_net_retention": "|C_lower+C_upper|/(|C_lower|+|C_upper|)",
        "gate_net_retention": "|G_lower+G_upper|/(|G_lower|+|G_upper|)",
        "warning": "A net sign or net slope alone cannot identify either local amplitude trend.",
    },
    "H2_static_largest_vs_gate_sensitive": {
        "result": "supported as a possibility, not established as a universal law",
        "screened_external_settings": len(external_cases),
        "screened_external_configuration_count": len(
            {str(row["configuration"]) for row in external_cases}
        ),
        "screened_external_rank_mismatches": sum(
            bool(row["static_gate_rank_mismatch"]) for row in external_cases
        ),
        "diagnostic_internal_step_settings": len(internal_cases),
        "diagnostic_internal_configuration_count": len(
            {str(row["configuration"]) for row in internal_cases}
        ),
        "diagnostic_internal_rank_mismatches": sum(
            bool(row["static_gate_rank_mismatch"]) for row in internal_cases
        ),
        "by_configuration_and_field": {name: aggregate(rows) for name, rows in sorted(groups.items())},
        "independence_limits": [
            "Band/projection settings are numerical repeats, not independent samples.",
            "InSe A and B are m001-related domains, not independent material families.",
            "MoS2 is tested only with a diagnostic internal energy step, not matching self-consistent +/-0.10 V/nm SCFs.",
            "All full-interface static ranks occur near layer compensation; InSe has an especially small static rank margin. Repetition within this archive is not a convergence proof or a robust material ordering.",
        ],
        "counterexamples": {
            "strong_KM_rank_matches": sum(
                not bool(row["static_gate_rank_mismatch"]) for row in strong_km
            ),
            "strong_KM_cases": len(strong_km),
            "representative_material_intervention_rank_matches": sum(
                not bool(row["static_gate_rank_mismatch"]) for row in interventions
            ),
            "representative_material_intervention_cases": len(interventions),
            "meaning": "The observed physical-interface mismatch is not an algebraic identity and need not survive channel interventions or parameter changes.",
        },
    },
    "quantitative_examples": {
        "MoSSe_7269_external": {
            "static_largest": mosse_external["static_largest_layer"],
            "gate_sensitive": mosse_external["gate_sensitive_layer"],
            "upper_share_of_absolute_gate_slope": 1 - f(mosse_external["gate_lower_absolute_share"]),
            "static_net_retention": mosse_external["static_net_retention"],
            "gate_net_retention": mosse_external["gate_net_retention"],
        },
        "InSe_864_A_external": {
            "static_largest": inse_a_external["static_largest_layer"],
            "gate_sensitive": inse_a_external["gate_sensitive_layer"],
            "upper_share_of_absolute_gate_slope": 1 - f(inse_a_external["gate_lower_absolute_share"]),
            "static_net_retention": inse_a_external["static_net_retention"],
            "gate_net_retention": inse_a_external["gate_net_retention"],
        },
        "InSe_864_B_external": {
            "static_largest": inse_b_external["static_largest_layer"],
            "gate_sensitive": inse_b_external["gate_sensitive_layer"],
            "lower_share_of_absolute_gate_slope": inse_b_external["gate_lower_absolute_share"],
            "static_net_retention": inse_b_external["static_net_retention"],
            "gate_net_retention": inse_b_external["gate_net_retention"],
        },
    },
    "testable_rules_and_status": [
        {
            "rule": "The zero-gate sign does not fix whether positive control strengthens or weakens a layer.",
            "status": "supported by two controlled pairs",
            "evidence": {
                "InSe_domains": {
                    "A_directions": [inse_a_external["direction_lower"], inse_a_external["direction_upper"]],
                    "A_trends": [
                        inse_a_external["positive_control_amplitude_trend_lower"],
                        inse_a_external["positive_control_amplitude_trend_upper"],
                    ],
                    "B_directions": [inse_b_external["direction_lower"], inse_b_external["direction_upper"]],
                    "B_trends": [
                        inse_b_external["positive_control_amplitude_trend_lower"],
                        inse_b_external["positive_control_amplitude_trend_upper"],
                    ],
                    "A_and_B_net_G": [inse_a_external["G_net_recomputed"], inse_b_external["G_net_recomputed"]],
                },
                "KM_relative_channel_phase": {
                    "alpha_plus_half_directions": [km_plus["direction_lower"], km_plus["direction_upper"]],
                    "alpha_plus_half_trends": [
                        km_plus["positive_control_amplitude_trend_lower"],
                        km_plus["positive_control_amplitude_trend_upper"],
                    ],
                    "alpha_minus_half_directions": [km_minus["direction_lower"], km_minus["direction_upper"]],
                    "alpha_minus_half_trends": [
                        km_minus["positive_control_amplitude_trend_lower"],
                        km_minus["positive_control_amplitude_trend_upper"],
                    ],
                },
            },
        },
        {
            "rule": "Static cancellation and gate-slope cancellation are separate diagnostics.",
            "status": "supported; no universal conversion inferred",
            "evidence": "MoSSe and InSe both have strongly cancelling static layer responses, but their retained net gate slopes differ sharply; see quantitative_examples and case_metrics.csv.",
        },
        {
            "rule": "A selected interlayer channel may change response direction, gate-sensitive layer, or gain direction separately.",
            "status": "conditional mechanism, not universal",
            "evidence": "MoSSe and InSe representative interventions in intervention_metrics.csv; the same directional channel class has opposite sensitivity effects.",
        },
        {
            "rule": "A quadratic signed-channel interpolation is not a universal material criterion.",
            "status": "falsified by retained InSe tests",
            "criterion_test_counts": criterion_boundary,
            "mixing_diagnostic": mixing,
            "boundary": "At alpha=1, InSe has 0/8 correct layer gain-direction classifications in each field definition, whereas MoSSe and MoS2 retain all available classifications. InSe's larger eta is post hoc and is not a threshold.",
        },
    ],
    "next_falsifiable_material_prediction": {
        "file": relative(SOURCES["prediction_freeze"]),
        "sha256": sha256(SOURCES["prediction_freeze"]),
        "target": "MoSSe BiDB 7268",
        "prediction": "C_lower<0, C_upper>0; upper is statically larger, but lower has larger absolute screened +/-0.10 V/nm gate slope; no layer sign reversal in that interval.",
        "confidence": "low heuristic; 7266 supplies only a protocol-mismatched static sign/rank cue, and existing compensated static rank margins can be comparable to numerical sensitivity",
        "important_nonprediction": "Whether positive gate strengthens or weakens the two amplitudes is intentionally left unpredicted.",
    },
    "limitations": [
        "The material values use archived 6x6 response grids and finite band windows; they are exploratory, not final material constants.",
        "InSe net values are cancellation-sensitive and changed materially with mesh even though layer directions and control trends survived.",
        "Internal layer-energy steps cannot be converted into external V/nm by a universal scalar; real self-consistent gate states change more than the mean layer offset.",
        "Hamiltonian channel removal/scaling is diagnostic and was not rescreened; it is not a directly available experimental knob.",
        "Absolute gate-sensitive layer and relative fractional sensitivity answer different questions; near-zero C can make fractional sensitivities misleading.",
        "No experimental detectability, orbital magnetic response, ionic relaxation, or universal threshold is established here.",
        "The historical 7266 direct finite-field response is not quantitatively interchangeable with the frozen 7268 independently screened +/-0.10 V/nm protocol.",
    ],
}


case_path = HERE / "case_metrics.csv"
intervention_path = HERE / "intervention_metrics.csv"
summary_path = HERE / "analysis_summary.json"
write_csv(case_path, cases)
write_csv(intervention_path, interventions)
summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

manifest = {
    "inputs": {
        key: {"path": relative(path), "sha256": sha256(path)} for key, path in SOURCES.items()
    },
    "outputs": {
        path.name: sha256(path)
        for path in (Path(__file__), case_path, intervention_path, summary_path)
    },
    "reproduction": "python3 studies/bsf-nmee-connection/results/screening/parallel-idea-loop-20261110/physics/analyze.py",
}
(HERE / "source_manifest.json").write_text(
    json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
)

print(json.dumps({"cases": len(cases), "interventions": len(interventions), "summary": str(summary_path)}))
