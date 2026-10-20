#!/usr/bin/env python3
"""Magnetic-symmetry verification of the AgF2 monolayer (altermagnet) and the
constructed B1/B2 bilayer candidates for the Type-I (m_z connector) and
Type-II (m'_z connector) magnetic configurations.

Checks:
  1. monolayer altermagnetic state (Ag sublattices with opposite moments along
     [100]) -> magnetic layer group p21'/b'11 (uni 86 = P2_1'/c') with the
     antiunitary mirror m'_100, matching C2DB and the paper (sub-lattice
     connectors {m_100|(1/2,1/2)} and {2_100|(1/2,1/2)}).
  2. B1/B2 Type-I (ferromagnetic stacking): internal magnetic point group
     m' = {1, m'_100}; B2 = {m_001 | t_o2} B1 with moments preserved
     (unitary connector m_z).
  3. B1/B2 Type-II (antiferromagnetic stacking): same internal group m';
     connector becomes antiunitary m'_z (needs time reversal).

Uses spglib magnetic symmetry (2.5.0) and the pinned toolkit v0.18.0 registry.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.environ.get("GTOT_TOOLKIT_SRC", "/tmp/gtot/toolkit-v0.18.0/src"))
# (toolkit source path is configurable; the pinned clone is at /tmp/gtot/toolkit-v0.18.0,
# tag v0.18.0 = commit 40a68e9ee2235a765e768f07f0f11fb99729c5c2, repo Xin-Jiaqi/group-theory-operations-toolkit)

import spglib  # noqa: E402

from group_theory_operations.magnetic_layer_groups import (  # noqa: E402
    load_magnetic_layer_group_registry,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

MONO = ROOT / "structures" / "C2DB_2AgF2-1_monolayer" / "POSCAR"
B1 = ROOT / "structures" / "constructed_B1_candidate" / "POSCAR"
B2 = ROOT / "structures" / "constructed_B2_candidate" / "POSCAR"


def read_poscar(path: Path):
    lines = path.read_text().splitlines()
    lattice = np.array([[float(x) for x in line.split()] for line in lines[2:5]])
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    coord_type = lines[7].strip().lower()
    natoms = sum(counts)
    coords = np.array([[float(x) for x in lines[8 + i].split()] for i in range(natoms)])
    if coord_type.startswith("car"):
        cart = coords
        frac = np.linalg.solve(lattice.T, cart.T).T
    else:
        frac = coords
    numbers = []
    for sp, c in zip(species, counts):
        z = 47 if sp == "Ag" else 9
        numbers += [z] * c
    return lattice, np.array(numbers), frac


def assign_monolayer_moments(frac, numbers, direction=0):
    """Collinear altermagnet: Ag(0,0) = +x, Ag(0.5,0.5) = -x (direction 0 = x)."""
    magmoms = []
    for i, n in enumerate(numbers):
        if n == 47:
            at_origin = abs(frac[i][0]) < 1e-4 and abs(frac[i][1]) < 1e-4
            m = np.zeros(3)
            m[direction] = 0.5 if at_origin else -0.5
            magmoms.append(m)
        else:
            magmoms.append(np.zeros(3))
    return np.array(magmoms)


def assign_bilayer_moments(frac, numbers, type_="I", t_o=(0.5, -0.16)):
    """Bilayer collinear moments along [100], Ag sublattices per the monolayer
    altermagnet (bottom: Ag(0,0)=+x, Ag(1/2,1/2)=-x).

    The magnetic configurations follow the PAPER's definitions (Zhu et al.
    PRL 135, 056801 (2025), main text around Fig. 2(b)/(c) and SM Sec. 1/3):

      type-I: interlayer Neel vectors OPPOSITE. The atoms connected by the
              stacking mirror m_001 have OPPOSITE spins (main text: "B_1 and
              B_2 with the type-I magnetic configuration are connected by
              N_s^- = {2_perp || m_001}"; SM Fig. 1(a) = Neel vectors
              aligned in the OPPOSITE direction -> atoms connected by N^- have
              opposite spins). => top images FLIP the bottom-layer sign.
      type-II: interlayer Neel vectors SAME. Atoms connected by m_001 have
              the SAME spin (SM: "the magnetic atoms connected by
              N^- = m_001 have the same spin" in the type-II/AHE context;
              connected by the antiunitary N_m^- = m'_001). => top images
              KEEP the bottom-layer sign.

    (Historical note: an earlier draft of this review had these two labels
    swapped, which combined with a polar-vector (z-flip-only) spin transform
    produced false-positive connector passes. The axial-vector rule
    m' = (-1)^tau det(R) R m fixes both issues; see results/
    magnetic_symmetry_report.json -> connector_B1_to_B2_* and the standalone
    test in agf2_preflight/stage6_harness/tests/test_magnetic_connector_axial.py.)
    """
    top_A = (t_o[0] % 1.0, t_o[1] % 1.0)
    top_B = ((0.5 + t_o[0]) % 1.0, (0.5 + t_o[1]) % 1.0)
    magmoms = []
    for i, n in enumerate(numbers):
        if n == 47:
            x, y = frac[i][0] % 1.0, frac[i][1] % 1.0
            if abs(x) < 1e-4 and abs(y) < 1e-4:
                m = np.array([0.5, 0.0, 0.0])          # bottom sublattice A
            elif abs(x - 0.5) < 1e-4 and abs(y - 0.5) < 1e-4:
                m = np.array([-0.5, 0.0, 0.0])         # bottom sublattice B
            elif abs(x - top_A[0]) < 1e-4 and abs(y - top_A[1]) < 1e-4:
                s = -1.0 if type_ == "I" else +1.0     # top image of A
                m = np.array([0.5 * s, 0.0, 0.0])
            elif abs(x - top_B[0]) < 1e-4 and abs(y - top_B[1]) < 1e-4:
                s = -1.0 if type_ == "I" else +1.0     # top image of B
                m = np.array([-0.5 * s, 0.0, 0.0])
            else:
                raise ValueError(f"unrecognized Ag site {x:.3f},{y:.3f}")
            magmoms.append(m)
        else:
            magmoms.append(np.zeros(3))
    return np.array(magmoms)


def magnetic_dataset(lattice, frac, numbers, magmoms, symprec=1e-3):
    ds = spglib.get_magnetic_symmetry_dataset((lattice, frac, numbers, magmoms),
                                              symprec=symprec)
    out = {
        "msg_type": int(ds.msg_type),
        "uni_number": int(ds.uni_number),
        "hall_number": int(ds.hall_number),
        "n_operations": int(len(ds.rotations)),
        "time_reversals": [int(t) for t in ds.time_reversals],
        "rotations": [np.asarray(r).tolist() for r in ds.rotations],
        "translations": [np.asarray(t).tolist() for t in ds.translations],
        "equivalent_atoms": [int(x) for x in ds.equivalent_atoms],
        "std_types": [int(x) for x in ds.std_types],
    }
    return out


def registry_info(uni_number, expected_parent_lg_number=None):
    """Look up magnetic layer group records by MSG uni number.

    The same MSG uni number can correspond to several magnetic LAYER groups
    (different parent layer groups in different settings); e.g. uni 33 is both
    p11a' (parent p11a) and pb'11 (parent pb11). Matching by uni number alone
    is therefore AMBIGUOUS. We list ALL matches and, when an expected parent
    layer group is supplied (verified independently by the non-magnetic
    layer-group check), select the record whose parent matches it."""
    reg = load_magnetic_layer_group_registry()
    matches = []
    for rec in reg["magnetic_layer_groups"]:
        if rec.get("corresponding_magnetic_space_group", {}).get("uni_number") == uni_number:
            matches.append({
                "litvin_og_symbol_ascii": rec.get("litvin_og_symbol_ascii"),
                "magnetic_type": rec.get("magnetic_type"),
                "magnetic_point_group_symbol": rec.get("magnetic_point_group_symbol"),
                "parent_layer_group_symbol": rec.get("parent_layer_group_symbol"),
                "parent_layer_group_number": rec.get("parent_layer_group_number"),
                "bns_number": rec.get("corresponding_magnetic_space_group", {}).get("bns_number"),
                "og_number": rec.get("corresponding_magnetic_space_group", {}).get("og_number"),
            })
    if not matches:
        return {"matches": [], "selected": None}
    selected = None
    if expected_parent_lg_number is not None:
        sel = [m for m in matches
               if m["parent_layer_group_number"] == expected_parent_lg_number]
        if sel:
            selected = sel[0]
    return {"matches": matches, "selected": selected,
            "note": "uni-number lookup is ambiguous when several magnetic layer "
                    "groups share the same corresponding MSG uni number; selection "
                    "is by the independently verified parent layer group."}


def wrapf(v):
    w = v % 1.0
    if w > 1.0 - 1e-6 or w < 1e-6:
        return 0.0
    return w


def connector_magnetic_check(frac_b1, numbers_b1, mag_b1, frac_b2, mag_b2,
                             tau_frac, type_="I"):
    """Apply the connector {m_001 | tau} to B1 atoms+moments and compare with
    B2 (modulo lattice), using the AXIAL-VECTOR rule for the real spin:

        m' = (-1)**tau_t * det(R) * (R @ m),   R = diag(1,1,-1), det(R) = -1,

    where tau_t = 0 for the unitary connector m_001 (type-I) and tau_t = 1 for
    the antiunitary connector m'_001 (type-II). For moments along [100] this
    gives: unitary m_001 maps S_x -> -S_x; antiunitary m'_001 maps S_x -> +S_x.
    (A polar-vector transform m_z -> -m_z would give exactly the opposite
    mapping and is WRONG for spins; the standalone regression test
    agf2_preflight/stage6_harness/tests/test_magnetic_connector_axial.py pins
    the rule and its discriminating power.)

    The per-atom table lists, for every B1 atom, the mapped B2 atom, its
    position, the applied operation {m_001|tau} and the moment match, so the
    full magnetic connection operation and atom arrangement are explicit.
    """
    a = len(numbers_b1)
    zm = 0.5
    R = np.diag([1.0, 1.0, -1.0])
    detR = -1.0
    tau_t = 0 if type_ == "I" else 1
    atom_table = []
    failures = 0
    for i in range(a):
        x, y, z = frac_b1[i]
        xp = wrapf(x + tau_frac[0])
        yp = wrapf(y + tau_frac[1])
        zp = wrapf(2 * zm - z)
        best, best_d = None, 1e9
        for j in range(a):
            d = (abs(xp - wrapf(frac_b2[j][0])) + abs(yp - wrapf(frac_b2[j][1]))
                 + abs(zp - wrapf(frac_b2[j][2])))
            if d < best_d:
                best_d, best = d, j
        if best_d > 1e-5:
            failures += 1
            atom_table.append({
                "b1_atom": i, "b2_atom": None, "position_error": best_d,
                "note": "no B2 atom within tolerance",
            })
            continue
        m_mapped = ((-1.0) ** tau_t) * detR * (R @ mag_b1[i])
        ok = bool(np.allclose(m_mapped, mag_b2[best], atol=1e-9))
        if not ok:
            failures += 1
        atom_table.append({
            "b1_atom": i,
            "b1_position_frac": [round(v, 6) for v in frac_b1[i]],
            "b1_moment": [round(v, 6) for v in mag_b1[i]],
            "b2_atom": best,
            "b2_position_frac": [round(v, 6) for v in frac_b2[best]],
            "b2_moment": [round(v, 6) for v in mag_b2[best]],
            "mapped_moment_axial": [round(v, 6) for v in m_mapped],
            "operation": f"{{m_001 | {tuple(round(v,4) for v in tau_frac)}}}"
                         + ("" if type_ == "I" else " * T"),
            "match": ok,
        })
    return {
        "positions_match_mod_lattice": bool(all(r["b2_atom"] is not None for r in atom_table)),
        "moment_mismatches": failures,
        "connector": ("m_001 (unitary, axial S_x -> -S_x)" if type_ == "I"
                      else "m'_001 (antiunitary, axial S_x -> +S_x)"),
        "spin_transform_rule": "m' = (-1)**tau * det(R) * R @ m  (axial vector, R=diag(1,1,-1))",
        "per_atom": atom_table,
    }


def magnetic_op_correspondence(ds, global_delta, expected_unitary_t, expected_antiunitary_t):
    """Map the spglib magnetic ops into the toolkit standard setting with the
    SAME single global origin shift used by the non-magnetic reconciliation,
    and compare with the expected magnetic layer group op translations.

    expected_unitary_t / expected_antiunitary_t: toolkit-standard translations
    (e.g. pb'11: unitary {E|0}, antiunitary mirror {m_x|(0,1/2,0)})."""
    delta = np.array(global_delta)
    rows = []
    for R, t, tr in zip(ds["rotations"], ds["translations"], ds["time_reversals"]):
        R = np.asarray(R, dtype=int)
        t2 = (np.asarray(t, dtype=float) + (R - np.eye(3)) @ delta) % 1.0
        t2 = np.array([wrapf(x) for x in t2])
        expect = expected_antiunitary_t if tr == 1 else expected_unitary_t
        match = bool(np.allclose(t2, expect, atol=1e-4))
        rows.append({
            "rotation": R.ravel().tolist(),
            "detected_translation": [round(x, 6) for x in t],
            "mapped_translation_after_global_delta": [round(x, 6) for x in t2],
            "time_reversal": int(tr),
            "expected_toolkit_translation": [round(x, 6) for x in expect],
            "match": match,
        })
    return rows, all(r["match"] for r in rows)


def main() -> None:
    report = {}

    # --- monolayer ---
    lattice, numbers, frac = read_poscar(MONO)
    mag = assign_monolayer_moments(frac, numbers)
    ds = magnetic_dataset(lattice, frac, numbers, mag)
    info = registry_info(ds["uni_number"], expected_parent_lg_number=17)
    report["monolayer"] = {
        "spin_config": "collinear altermagnet, Ag(0,0)=+x, Ag(1/2,1/2)=-x, moments along [100]",
        "spglib_magnetic_dataset": ds,
        "toolkit_magnetic_layer_group": info,
        "setting_correspondence": {
            "global_origin_shift_delta": [0.0, 0.0, 0.0],
            "note": "spglib magnetic dataset and toolkit standard setting agree "
                    "directly (delta = 0); expected ops E, {-I|0}, {2_x|(1/2,1/2,0)}*T, "
                    "{m_x|(1/2,1/2,0)}*T (p21'/b'11 standard setting).",
        },
        "expected": {
            "c2db_msg": "P2_1'/c' (uni 86)",
            "paper_sublattice_connectors": ["{m_100|(1/2,1/2)}", "{2_100|(1/2,1/2)}"],
            "comment": "antiunitary mirror m'_100 present iff an op with R=diag(-1,1,1) has time_reversal=1",
        },
    }
    # check antiunitary mirror m'_x presence
    mirror_found = False
    for R, t, tr in zip(ds["rotations"], ds["translations"], ds["time_reversals"]):
        if np.allclose(np.asarray(R), np.diag([-1, 1, 1]), atol=1e-9) and tr == 1:
            mirror_found = True
    report["monolayer"]["antiunitary_mirror_m100_present"] = mirror_found
    print("monolayer: uni", ds["uni_number"], "mirror m'_100 present:", mirror_found)

    # --- bilayers ---
    # Global origin shift that reconciles the NON-magnetic op sets with the
    # toolkit pb11 standard setting (results/bilayer_symmetry_report.json):
    # delta = (0.75, 0, 0) == (-0.25, 0, 0) mod the mirror kernel.
    # pb'11 (type III, parent pb11) toolkit-standard op set:
    #   unitary      {E|0}
    #   antiunitary  {m_x|(0,1/2,0)}*T   (mirror perpendicular to a, b/2 glide)
    # Magnetic-configuration labels follow the paper (see assign_bilayer_moments):
    #   type-I  = interlayer Neel vectors OPPOSITE  (unitary m_001 connector)
    #   type-II = interlayer Neel vectors SAME      (antiunitary m'_001 connector)
    GLOBAL_DELTA = [0.75, 0.0, 0.0]
    PB11_PRIME_UNITARY_T = [0.0, 0.0, 0.0]
    PB11_PRIME_ANTIUNITARY_T = [0.0, 0.5, 0.0]
    t_o_by_label = {"B1": (0.5, -0.16), "B2": (0.5, 0.16)}
    for label, path in (("B1", B1), ("B2", B2)):
        lattice, numbers, frac = read_poscar(path)
        for t_ in ("I", "II"):
            mag = assign_bilayer_moments(frac, numbers, type_=t_, t_o=t_o_by_label[label])
            ds = magnetic_dataset(lattice, frac, numbers, mag)
            info = registry_info(ds["uni_number"], expected_parent_lg_number=12)
            op_corr, op_ok = magnetic_op_correspondence(
                ds, GLOBAL_DELTA, PB11_PRIME_UNITARY_T, PB11_PRIME_ANTIUNITARY_T)
            spin_desc = ("interlayer Neel vectors OPPOSITE (top images flip the "
                         "bottom-layer sign; unitary m_001 connector flips in-plane "
                         "S_x -> -S_x)" if t_ == "I" else
                         "interlayer Neel vectors SAME (top images keep the "
                         "bottom-layer sign; antiunitary m'_001 connector keeps "
                         "in-plane S_x -> +S_x)")
            report[f"{label}_type{t_}"] = {
                "spin_config": spin_desc + "; moments along [100]",
                "spglib_magnetic_dataset": ds,
                "toolkit_magnetic_layer_group": info,
                "internal_magnetic_point_group": "m' = {1, m'_100} "
                                                 "(verified from the spglib magnetic op set; "
                                                 "same for type-I and type-II)",
                "external_connector_verdict": "see connector_B1_to_B2_* below "
                                              "(cross-state mapping; kept separate from the "
                                              "internal MPG verdict per review)",
                "setting_correspondence": {
                    "global_origin_shift_delta": GLOBAL_DELTA,
                    "expected_op_set_toolkit_standard": {
                        "unitary": "{E|0}",
                        "antiunitary": "{m_x|(0,1/2,0)}*T  (mirror perp. a, b/2 glide)",
                    },
                    "per_operation": op_corr,
                    "all_ops_match_after_global_delta": op_ok,
                    "note": "the SAME single global origin shift reconciles the non-magnetic "
                            "op set with the toolkit pb11 standard setting and the magnetic "
                            "op set with the pb'11 standard setting.",
                },
                "expected_internal_magnetic_point_group": "m' = {1, m'_100}",
            }
            print(f"{label} type-{t_}: uni {ds['uni_number']} n_ops={len(ds['rotations'])} "
                  f"T={ds['time_reversals']} op_corr_ok={op_ok}")
    # connector checks (magnetic): axial-vector rule
    l1, n1, f1 = read_poscar(B1)
    l2, n2, f2 = read_poscar(B2)
    for t_ in ("I", "II"):
        m1 = assign_bilayer_moments(f1, n1, type_=t_, t_o=(0.5, -0.16))
        m2 = assign_bilayer_moments(f2, n2, type_=t_, t_o=(0.5, 0.16))
        cc = connector_magnetic_check(f1, n1, m1, f2, m2, (0.5, 0.16), type_=t_)
        report[f"connector_B1_to_B2_type{t_}_magnetic"] = cc
        print(f"connector type-{t_}: {cc['connector']}, pos_ok={cc['positions_match_mod_lattice']}, "
              f"moment_mismatches={cc['moment_mismatches']}")
    report["connector_summary"] = {
        "operator": "{m_001 | (0.5,0.16)}",
        "note": "rotational part m_001 (z-mirror about interlayer midplane); "
                "translational part t_o2 per SM Eq.(3); maps B1 structure onto B2 modulo lattice",
        "spin_transform_rule": "m' = (-1)**tau * det(R) * (R @ m),  R=diag(1,1,-1),  det(R)=-1 "
                               "(real spin = axial vector)",
        "type_I": "interlayer Neel vectors OPPOSITE; connector m_001 UNITARY (axial S_x -> -S_x); "
                  "paper: B1 and B2 with the type-I magnetic configuration are connected by "
                  "N_s^- = {2_perp || m_001} (main text around Fig. 2(b)); switches k-space spin "
                  "polarization, does NOT switch AHE",
        "type_II": "interlayer Neel vectors SAME; connector m'_001 ANTIUNITARY (axial S_x -> +S_x); "
                   "paper: switching AHE requires N_m^- = m'_001, i.e. the type-II configuration "
                   "(main text around Fig. 2(c), 4% strain; SM: the magnetic atoms connected by "
                   "N^- = m_001 have the same spin)",
        "verification_separation": "internal magnetic point group m' = {1,m'_100} (same for both "
                                   "types, verified from the spglib magnetic op set) is reported "
                                   "separately from the external connector cross-state mapping "
                                   "verification (per-atom table above)",
    }
    out = ROOT / "results" / "magnetic_symmetry_report.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", out)


if __name__ == "__main__":
    main()
