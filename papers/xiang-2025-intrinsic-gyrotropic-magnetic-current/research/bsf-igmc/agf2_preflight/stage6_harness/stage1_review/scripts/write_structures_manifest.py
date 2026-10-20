#!/usr/bin/env python3
"""Write results/STRUCTURES_MANIFEST.json: SHA-256 of every structure file and
analysis result produced by the Stage-1 review, with role annotations."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

ROLE = {
    "structures/C2DB_2AgF2-1_monolayer/POSCAR": "monolayer source (C2DB 2AgF2-1, CC BY-SA 4.0)",
    "structures/C2DB_2AgF2-1_monolayer/AgF2_monolayer.cif": "monolayer CIF (visualization)",
    "structures/C2DB_2AgF2-1_monolayer/provenance.json": "monolayer provenance",
    "structures/constructed_B1_candidate/POSCAR": "B1 CONSTRUCTED CANDIDATE (unrelaxed)",
    "structures/constructed_B1_candidate/provenance.json": "B1 provenance/construction record",
    "structures/constructed_B2_candidate/POSCAR": "B2 CONSTRUCTED CANDIDATE (unrelaxed)",
    "structures/constructed_B2_candidate/provenance.json": "B2 provenance/construction record",
}


def main() -> None:
    out = ROOT / "results" / "STRUCTURES_MANIFEST.json"
    entries = {}
    for p in sorted(list((ROOT / "structures").rglob("*")) + list((ROOT / "results").glob("*.json"))):
        if not p.is_file() or p.resolve() == out.resolve():
            # the manifest itself is excluded from its own entries: a
            # self-referential hash can never be a fixed point, and excluding
            # it makes regeneration byte-identical on a blank checkout
            continue
        rel = str(p.relative_to(ROOT))
        entries[rel] = {
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "size_bytes": p.stat().st_size,
            "role": ROLE.get(rel, "analysis artifact"),
        }
    manifest = {
        "schema_version": "stage1_structures_manifest_v1",
        "generated": "2026-10-20",
        "note": "Structure files are either the public C2DB monolayer (CC BY-SA 4.0) or "
                "clearly-labeled constructed unrelaxed candidates; none is the authors' "
                "relaxed B1/B2 (not public per PRL data-availability statement). This "
                "manifest itself is excluded from its own entries (self-referential "
                "hashes are not a fixed point); regeneration is byte-identical.",
        "files": entries,
    }
    out.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {out} ({len(entries)} files)")


if __name__ == "__main__":
    main()
