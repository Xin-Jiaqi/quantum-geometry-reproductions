#!/usr/bin/env python3
"""Construct AgF2 bilayer B1/B2 CANDIDATES from the C2DB monolayer by the
published stacking rule of Zhu et al. PRL 135, 056801 (2025).

Rule (paper main text / SM Sec. 1):
    B = S + Ohat S,  Ohat = { m_001 | t_o },
    m_001 = mirror perpendicular to [001] (z -> 2*z_m - z about the interlayer
            mid-plane z_m), t_o = in-plane fractional translation of the top
            layer relative to the bottom layer:
              B1: t_o1 = (0.5, -0.16)
              B2: t_o2 = (0.5, +0.16)

Coordinate conventions (documented, per HARNESS_TASK A.3):
  * lateral cell = monolayer cell a=(3.8270,0,0) A, b=(0,5.8797,0) A
    (C2DB 2AgF2-1, LG17 p2_1/b11), same for both layers (paper: same lateral
    cell as the single layer).
  * z axis out of plane; mirror plane at z_m = Lz/2 (cell centre); bottom
    layer on the +z side of the mirror plane, top layer its mirror image.
  * t_o applied to the top layer in fractional in-plane coordinates:
      Cartesian shift = (t_o.x * |a|, t_o.y * |b|, 0).
  * INTERLAYER DISTANCE IS NOT PUBLISHED: nominal F-F interlayer gap
    g = 3.0 A (vdW F-F contact estimate) is ASSUMED for the unrelaxed
    candidate; d(Ag-Ag) = g + 2*delta, delta = F-Ag offset in monolayer.
  * Vacuum per paper protocol: 20 A between periodic images along z.

These are CONSTRUCTED, UNRELAXED CANDIDATES, NOT the authors' relaxed B1/B2
structures (which are not publicly available; PRL data-availability statement).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

A_VEC = (3.827024970065565, 0.0, 0.0)
B_VEC = (0.0, 5.879665745808328, 0.0)
MONO_POSCAR = ROOT / "structures" / "C2DB_2AgF2-1_monolayer" / "POSCAR"

GAP_F_F = 3.0  # Assumed nominal interlayer F-F gap (Angstrom), unrelaxed candidate
VACUUM = 20.0  # Author protocol: 20 A vacuum


def read_poscar(path: Path):
    lines = path.read_text().splitlines()
    lattice = [[float(x) for x in line.split()] for line in lines[2:5]]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    coord_type = lines[7].strip().lower()
    natoms = sum(counts)
    coords = []
    for i in range(8, 8 + natoms):
        coords.append([float(x) for x in lines[i].split()])
    numbers = []
    for sp, c in zip(species, counts):
        z = 47 if sp == "Ag" else 9
        numbers += [z] * c
    return lattice, numbers, coords, coord_type


def main() -> None:
    lattice, numbers, cart, coord_type = read_poscar(MONO_POSCAR)
    assert coord_type.startswith("car")
    assert numbers == [47, 47, 9, 9, 9, 9]  # Ag Ag F F F F
    ag_z = cart[0][2]  # Ag z (both Ag at same z for this monolayer)
    assert abs(cart[1][2] - ag_z) < 1e-6
    delta = abs(cart[2][2] - ag_z)
    print(f"monolayer Ag z={ag_z:.6f}, delta={delta:.6f}")

    d_ag = GAP_F_F + 2 * delta
    Lz = d_ag + 2 * delta + VACUUM
    zm = Lz / 2.0
    print(f"d(Ag-Ag)={d_ag:.6f} A, Lz={Lz:.6f} A, mirror z_m={zm:.6f} A")

    summary = {}
    for label, t_frac in (("B1", (0.5, -0.16)), ("B2", (0.5, +0.16))):
        t_cart = (t_frac[0] * A_VEC[0], t_frac[1] * B_VEC[1], 0.0)
        # bottom layer: translate monolayer so Ag plane sits at zm + d_ag/2
        z_ag_bottom = zm + d_ag / 2.0
        zshift = z_ag_bottom - ag_z
        bottom = [(x, y, z + zshift) for (x, y, z) in cart]
        # top layer = {m_001|t_o} acting on the bottom layer
        top = [(x + t_cart[0], y + t_cart[1], 2 * zm - z) for (x, y, z) in bottom]

        out_cell = [A_VEC, B_VEC, (0.0, 0.0, Lz)]
        species = ["Ag", "F"]
        out_numbers = []
        out_coords = []
        # VASP species-blocked ordering: ALL Ag first (bottom then top),
        # then ALL F (bottom then top), matching counts "4 8"
        for sp in species:
            znum = 47 if sp == "Ag" else 9
            idx = [i for i, n in enumerate(numbers) if n == znum]
            for i in idx:
                out_numbers.append(numbers[i])
                out_coords.append(bottom[i])
            for i in idx:
                out_numbers.append(numbers[i])
                out_coords.append(top[i])

        lines = [
            f"AgF2 bilayer {label} CANDIDATE (constructed from C2DB 2AgF2-1 monolayer via "
            f"Ohat={{m_001|({t_frac[0]},{t_frac[1]:+.2f})}}); UNRELAXED - NOT author-relaxed structure",
            "1.0",
            f"  {out_cell[0][0]:.16f}  {out_cell[0][1]:.16f}  {out_cell[0][2]:.16f}",
            f"  {out_cell[1][0]:.16f}  {out_cell[1][1]:.16f}  {out_cell[1][2]:.16f}",
            f"  {out_cell[2][0]:.16f}  {out_cell[2][1]:.16f}  {out_cell[2][2]:.16f}",
            " ".join(species),
            " ".join(str(out_numbers.count(47)) for _ in ["x"]) + " " + str(out_numbers.count(9)),
            "Cartesian",
        ]
        for x, y, z in out_coords:
            lines.append(f"  {x:.8f}  {y:.8f}  {z:.8f}")
        poscar_text = "\n".join(lines) + "\n"

        outdir = ROOT / "structures" / f"constructed_{label}_candidate"
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "POSCAR").write_text(poscar_text, encoding="utf-8")

        prov = {
            "schema_version": "stage1_structure_provenance_v1",
            "structure_label": f"AgF2 bilayer {label} CONSTRUCTED CANDIDATE (unrelaxed)",
            "construction": {
                "rule": "B = S + Ohat S; Ohat = {m_001 | t_o} (Zhu et al. PRL 135, 056801 (2025), "
                         "main text and SM Sec. 1)",
                "t_o_fractional": list(t_frac),
                "t_o_cartesian_ang": list(t_cart),
                "mirror_plane_z": zm,
                "layer_exchange": "top = {m_001|t_o} bottom (mirror about z_m then in-plane translation)",
                "same_lateral_cell": True,
                "pbc": [True, True, False],
                "assumptions": {
                    "interlayer_gap_f_f_ang": GAP_F_F,
                    "d_ag_ag_ang": d_ag,
                    "note": "interlayer distance not published; nominal vdW gap assumed; "
                            "must be relaxed in DFT and compared with source energetics",
                },
            },
            "monolayer_source": {
                "uid": "2AgF2-1",
                "layer_group": "p2_1/b11 (LG17)",
                "database": "C2DB (local dump project-cy/c2db.db, system 4195)",
                "license": "CC BY-SA 4.0",
                "disclaimer": "NOT the author-relaxed structure of Zhu et al. PRL 135, 056801 (2025); "
                              "that structure is not publicly available (PRL data-availability statement).",
            },
            "cell": {"a_ang": A_VEC[0], "b_ang": B_VEC[1], "c_ang": Lz,
                     "alpha_deg": 90.0, "beta_deg": 90.0, "gamma_deg": 90.0,
                     "vacuum_ang": VACUUM},
            "composition": {"Ag": 4, "F": 8, "formula": "Ag4F8", "formula_units": 4, "per_layer": 2},
            "hashes": {"POSCAR_sha256": hashlib.sha256(poscar_text.encode("utf-8")).hexdigest()},
        }
        (outdir / "provenance.json").write_text(json.dumps(prov, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        summary[label] = prov["hashes"]
        print(f"{label}: POSCAR written, sha256={prov['hashes']['POSCAR_sha256']}")
        print(f"  z-range bottom: {min(c[2] for c in bottom):.3f}..{max(c[2] for c in bottom):.3f}; "
              f"top: {min(c[2] for c in top):.3f}..{max(c[2] for c in top):.3f}")

    # connection check: {m_001 | tau} B1 == B2 with tau = t_o2 (=-t_o1 mod lattice)
    # (SM Eq. 3: t_o2 = -N^- t_o1 + t_0 ; the connection's translational part is t_o2)
    b1 = _load(ROOT / "structures" / "constructed_B1_candidate" / "POSCAR")
    b2 = _load(ROOT / "structures" / "constructed_B2_candidate" / "POSCAR")
    zm1 = b1["lattice"][2][2] / 2.0
    tau = (0.5, 0.16, 0.0)
    mapped = []
    for (x, y, z) in b1["coords"]:
        mapped.append((x + tau[0] * b1["lattice"][0][0],
                       y + tau[1] * b1["lattice"][1][1],
                       2 * zm1 - z))
    tol = 1e-5  # Angstrom; absorbs POSCAR 8-dp rounding
    a_len = b1["lattice"][0][0]
    b_len = b1["lattice"][1][1]
    ok = True
    for m_ in mapped:
        best = min(
            (abs(m_[0] - c[0] - ia * a_len) + abs(m_[1] - c[1] - ib * b_len) + abs(m_[2] - c[2])
             for c in b2["coords"] for ia in (-1, 0, 1) for ib in (-1, 0, 1))
        )
        if best > tol:
            ok = False
            break
    print("connection check B2 == {m_001 | (0.5,0.16)} B1 (mod lattice):",
          "PASS" if ok else "FAIL")


def _load(path: Path):
    lines = path.read_text().splitlines()
    lattice = [[float(x) for x in line.split()] for line in lines[2:5]]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    natoms = sum(counts)
    coords = [[float(x) for x in lines[8 + i].split()] for i in range(natoms)]
    numbers = []
    for sp, c in zip(species, counts):
        numbers += [47 if sp == "Ag" else 9] * c
    return {"lattice": lattice, "coords": coords, "numbers": numbers}


if __name__ == "__main__":
    main()
