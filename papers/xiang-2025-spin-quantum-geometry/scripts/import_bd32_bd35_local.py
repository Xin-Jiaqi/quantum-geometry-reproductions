#!/usr/bin/env python3
"""Safely import verified BD32-BD35 exploratory source snapshots into this repo.

No network, HPC, privileged operations, or production computing.
Usage (from any directory):
    python papers/xiang-2025-spin-quantum-geometry/scripts/import_bd32_bd35_local.py \
      --archive /path/to/spinqg_handoff_sources_bd32_bd35_20261023.zip

The ZIP is distributed separately through the user's ChatGPT conversation;
it is NOT in GitHub. The importer fails closed on an incorrect hash, unsafe
paths, duplicates, or an existing target; it never changes the paper
reproduction core, its reference data, or closed MnPS3 experiments.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

EXPECTED_ZIP_SHA256 = "605843ccf73d64d29f05da7181b915429167337a852f2b025f6adfc1c66d1602"
EXPECTED_SIZE_BYTES = 42998
EXPECTED_SCRIPTS = {
    "scripts/BD32_ptse2_inversion_registry_gate.py",
    "scripts/BD32_PT_trace_model_gate.py",
    "scripts/BD32_torus_index_gate.py",
    "scripts/BD33_periodic_PT_spin_response_gate.py",
    "scripts/BD34_interlayer_Se_sublayer_PT_gate.py",
    "scripts/BD35_independent_audit.py",
}
# scripts/BD34... and BD35... are historically *synthetic models*, not fitted material calculations.
ALLOWED_PREFIXES = ("scripts/", "reports/", "reference/")
PAPER_ROOT = Path(__file__).resolve().parents[1]
TARGET = PAPER_ROOT / "research" / "lee-ptse2-bd32-bd35"


def check(archive: Path) -> tuple[list[tuple[str, bytes]], str]:
    if not archive.is_file():
        raise ValueError("SOURCE_ARCHIVE_MISSING: " + str(archive))
    raw = archive.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_ZIP_SHA256 or len(raw) != EXPECTED_SIZE_BYTES:
        raise ValueError(
            "SOURCE_ARCHIVE_INTEGRITY_FAIL: size="
            + str(len(raw)) + ", sha256=" + digest
        )
    items: list[tuple[str, bytes]] = []
    with ZipFile(archive) as z:
        names = [x.filename for x in z.infolist() if not x.is_dir()]
        if len(names) != len(set(names)):
            raise ValueError("DUPLICATE_ZIP_MEMBERS")
        if not EXPECTED_SCRIPTS.issubset(set(names)):
            raise ValueError("SOURCE_ARCHIVE_MISSING_CRITICAL_SCRIPTS")
        for info in z.infolist():
            if info.is_dir():
                continue
            rel = PurePosixPath(info.filename)
            if (
                rel.is_absolute() or ".." in rel.parts
                or not any(info.filename.startswith(p) for p in ALLOWED_PREFIXES)
                or rel.suffix not in (".py", ".md", ".json", ".csv")
                or info.file_size > 1_000_000
            ):
                raise ValueError("UNSAFE_ZIP_PATH: " + info.filename)
            items.append((info.filename, z.read(info)))
    return items, digest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--verify-only", action="store_true")
    args = ap.parse_args()
    try:
        items, digest = check(args.archive)
        manifest = {
            "source_zip_sha256": digest,
            "source_size_bytes": EXPECTED_SIZE_BYTES,
            "science_status": "SYNTHETIC_MODEL_AND_SYMMETRY_ONLY",
            "original_file_sha256": {
                p: hashlib.sha256(b).hexdigest() for p, b in items
            },
        }
        print(f"SOURCE VERIFIED: {len(items)} files, sha256={digest}")
        if args.verify_only:
            return 0
        if TARGET.exists():
            raise ValueError("TARGET_ALREADY_EXISTS_NO_OVERWRITE: " + str(TARGET))
        # Archive was validated in memory before creating any output directories.
        TARGET.mkdir(parents=True, exist_ok=False)
        for filename, data in items:
            out = TARGET / filename
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(data)
        (TARGET / "IMPORT_MANIFEST.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print("IMPORTED into " + str(TARGET))
        print("Next: run BD32 site and torus scripts, then inspect BD35 caveats.")
        return 0
    except (ValueError, OSError, EOFError) as exc:
        print("IMPORT FAILED: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
