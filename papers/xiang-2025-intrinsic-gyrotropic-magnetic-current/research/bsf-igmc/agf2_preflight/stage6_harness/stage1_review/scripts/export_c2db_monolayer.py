#!/usr/bin/env python3
"""Export the C2DB AgF2 monolayer (2AgF2-1, layer group p2_1/b11, LG17) from
the locally stored C2DB SQLite dump into POSCAR + CIF + provenance JSON.

Input : local C2DB row dump JSON (path configurable via the environment
        variable AGF2_C2DB_ROW_JSON; default /tmp/agf2_2AgF2-1_structure.json)
        extracted from the local C2DB dump (system id 4195, UID 2AgF2-1).
Output: ../structures/C2DB_2AgF2-1_monolayer/{POSCAR,AgF2_monolayer.cif,provenance.json}

LOCAL_DATA_REQUIRED: this script cannot run from a blank checkout - it needs
the C2DB database row that is not part of this repository (C2DB itself is a
third-party database). The frozen POSCAR/CIF/provenance in the repository were
produced by this script and are hash-pinned (results/STRUCTURES_MANIFEST.json);
re-running it requires the local dump row.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "structures" / "C2DB_2AgF2-1_monolayer"

# Configurable input path (operator/environment choice, not fabricated data).
C2DB_ROW_JSON = Path(os.environ.get(
    "AGF2_C2DB_ROW_JSON", "/tmp/agf2_2AgF2-1_structure.json"))

NUM2SYM = {47: "Ag", 9: "F"}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    src = json.loads(C2DB_ROW_JSON.read_text())
    numbers = src["numbers"]
    cell = src["cell"]  # row-major 3x3
    pos = src["positions_cartesian"]  # n*3 flat
    natoms = len(numbers)
    coords = [tuple(pos[3 * i: 3 * i + 3]) for i in range(natoms)]

    # sanity
    assert src["label"] == "2AgF2-1 newdb"
    assert numbers.count(47) == 2 and numbers.count(9) == 4
    assert abs(cell[2]) < 1e-6 and abs(cell[5]) < 1e-6 and abs(cell[6]) < 1e-6 and abs(cell[7]) < 1e-6

    species = ["Ag", "F"]
    counts = [numbers.count(47), numbers.count(9)]
    order = []
    for sp in species:
        z = 47 if sp == "Ag" else 9
        for i, n in enumerate(numbers):
            if n == z:
                order.append(i)

    lines = [
        "AgF2 monolayer (C2DB 2AgF2-1, exfoliated02-21; LG17 p2_1/b11) -- constructed POSCAR from database row; NOT author-relaxed structure",
        "1.0",
        f"  {cell[0]:.16f}  {cell[1]:.16f}  {cell[2]:.16f}",
        f"  {cell[3]:.16f}  {cell[4]:.16f}  {cell[5]:.16f}",
        f"  {cell[6]:.16f}  {cell[7]:.16f}  {cell[8]:.16f}",
        " ".join(species),
        " ".join(str(c) for c in counts),
        "Cartesian",
    ]
    for i in order:
        x, y, z = coords[i]
        lines.append(f"  {x:.8f}  {y:.8f}  {z:.8f}")
    poscar_text = "\n".join(lines) + "\n"
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "POSCAR").write_text(poscar_text, encoding="utf-8")

    # CIF (for visualization only)
    cif_lines = [
        "data_AgF2_monolayer_C2DB_2AgF2-1",
        "_cell_length_a  %.8f" % cell[0],
        "_cell_length_b  %.8f" % cell[4],
        "_cell_length_c  %.8f" % cell[8],
        "_cell_angle_alpha  90.0",
        "_cell_angle_beta   90.0",
        "_cell_angle_gamma  90.0",
        "_symmetry_space_group_name_H-M  'P 1'",
        "_symmetry_Int_Tables_number  1",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
    ]
    for i in order:
        x, y, z = coords[i]
        fx, fy, fz = x / cell[0], y / cell[4], z / cell[8]
        cif_lines.append(
            f"{NUM2SYM[numbers[i]]}  {NUM2SYM[numbers[i]]}  {fx:.8f}  {fy:.8f}  {fz:.8f}"
        )
    (OUT / "AgF2_monolayer.cif").write_text("\n".join(cif_lines) + "\n", encoding="utf-8")

    provenance = {
        "schema_version": "stage1_structure_provenance_v1",
        "structure_label": "C2DB 2AgF2-1 monolayer",
        "source": {
            "database": "Computational 2D Materials Database (C2DB)",
            "uid": "2AgF2-1",
            "registry": "C2DB all-data registry (C2DB_all_data.csv row 4194)",
            "folder": "/home/niflheim2/cmr/C2DB-ASR/ICSD-COD/2el/Ag4F8",
            "protocol_label": "exfoliated02-21",
            "layer_group": "p2_1/b11",
            "layer_group_number": 17,
            "local_dump": "project-cy/c2db.db (ASE sqlite format), system id 4195",
            "local_dump_file_mtime": "2024-09-26",
            "license": "CC BY-SA 4.0 (per local DATA_LICENSES.md registry: "
                       "https://cmr.fysik.dtu.dk/c2db/c2db.html ; https://creativecommons.org/licenses/by-sa/4.0/)",
            "license_text": "CC BY-SA 4.0 (Creative Commons Attribution-ShareAlike 4.0 "
                            "International); full legal text at "
                            "https://creativecommons.org/licenses/by-sa/4.0/legalcode",
            "disclaimer": "Third-party PBE/GPAW-relaxed monolayer from the C2DB; "
                          "NOT the author-relaxed structure of Zhu et al. PRL 135, 056801 (2025).",
            "conversion_algorithm": "scripts/export_c2db_monolayer.py: reads the local C2DB row dump "
                                    "(project-cy/c2db.db, system id 4195, UID 2AgF2-1; input path "
                                    "configurable via env AGF2_C2DB_ROW_JSON, default "
                                    "/tmp/agf2_2AgF2-1_structure.json; numbers, 3x3 row-major cell, "
                                    "flat Cartesian positions); asserts composition Ag2F4 and an "
                                    "orthorhombic cell (a,b in-plane, z = vacuum slab height); "
                                    "reorders atoms into VASP species blocks [Ag, F]; writes "
                                    "Cartesian POSCAR (8 dp) and a visualization CIF (P1, no "
                                    "symmetry applied); hashes POSCAR, CIF and the raw database-row "
                                    "dump. Cell vectors and positions are copied verbatim from the "
                                    "database row - nothing is fabricated or symmetrized.",
            "original_row_verification": "source_row_uid_sha256 = SHA-256 of the raw JSON dump of "
                                         "the database row (system 4195); cross-checked against the "
                                         "2022-11-30 dump (system 10217, same cell, "
                                         "E=-19.3594 eV)",
            "retrieval_route": "independent retrieval: https://cmr.fysik.dtu.dk/c2db/ (search UID "
                               "2AgF2-1) or the registry C2DB_all_data.csv row 4194",
        },
        "cell": {
            "a_ang": cell[0], "b_ang": cell[4], "c_ang": cell[8],
            "alpha_deg": 90.0, "beta_deg": 90.0, "gamma_deg": 90.0,
            "vacuum_z_ang": cell[8] - (max(c[2] for c in coords) - min(c[2] for c in coords)),
            "pbc": [True, True, False],
        },
        "composition": {"Ag": 2, "F": 4, "formula": "Ag2F4", "formula_units": 2},
        "coordinates": {
            "mode": "Cartesian",
            "species_order_in_poscar": ["Ag", "F"],
            "numbers": numbers,
            "positions_cartesian": coords,
        },
        "hashes": {
            "POSCAR_sha256": sha256_bytes(poscar_text.encode("utf-8")),
            "CIF_sha256": sha256_bytes((OUT / "AgF2_monolayer.cif").read_bytes()),
            "source_row_uid_sha256": hashlib.sha256(json.dumps(src, sort_keys=True).encode()).hexdigest(),
        },
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("Wrote POSCAR, CIF, provenance under", OUT)
    print(json.dumps(provenance["hashes"], indent=2))


if __name__ == "__main__":
    main()
