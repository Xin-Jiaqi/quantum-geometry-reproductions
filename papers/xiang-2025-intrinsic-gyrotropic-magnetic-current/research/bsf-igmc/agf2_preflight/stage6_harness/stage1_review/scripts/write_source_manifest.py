#!/usr/bin/env python3
"""Write SOURCES_MANIFEST.json as a DECLARED registry of Zhu et al. (2025)
sources used by the Stage-1 review.

Registry only: does not download anything. Copyrighted artifacts (publisher/
arXiv PDFs, TeX bundle, extracted figures, full-text SM dump) are PRIVATE LOCAL
ARCHIVES: they are hashed here for provenance but are never committed to the
public repository (see sources/README.md and the PR #8 review, item 1).
Only URLs, SHA-256, short protocol excerpts and archive notes are published.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
# Private archives live OUTSIDE the repository tree (sibling of the repo root,
# inside the session workspace) so that release verifiers (verify_round9.py)
# never see copyrighted attachments in the public staging area. The manifest
# records their hashes and URLs but no absolute local paths.
# Configurable via env AGF2_PRIVATE_ARCHIVE_DIR (LOCAL_DATA_REQUIRED: the
# archives themselves are not part of the repository).
ARCHIVE_DIR = Path(os.environ.get(
    "AGF2_PRIVATE_ARCHIVE_DIR",
    str(Path(__file__).resolve().parents[9] / "private-agf2-stage1-archives")))

# (path_under_archive, url, role, availability)
DECLARED = [
    ("zhu2025_prl_published.pdf",
     "https://journals.aps.org/prl/abstract/10.1103/dmzg-ck2t",
     "PRL published paper (Zhu et al., PRL 135, 056801 (2025))",
     "private_local_only"),
    ("zhu2025_prl_sm.pdf",
     "https://journals.aps.org/prl/abstract/10.1103/dmzg-ck2t (Supplemental Material)",
     "PRL published Supplemental Material (SM Sec. 1-4)",
     "private_local_only"),
    ("zhu2025_arxiv_v2.pdf",
     "https://arxiv.org/abs/2502.17095",
     "arXiv preprint v2 (full derivation text)",
     "private_local_only"),
    ("zhu2025_arxiv_v2.tar.gz",
     "https://arxiv.org/abs/2502.17095",
     "arXiv TeX source bundle (redistributed nowhere)",
     "private_local_only"),
    ("zhu2025_prl_sm.txt",
     None,
     "Full-text pdftotext extraction of the SM (private working copy only; short excerpts only in results/agf2_protocol_frozen.json)",
     "private_local_only"),
    ("zhu2025_src/main_text-clean.tex",
     None,
     "Extracted arXiv main text TeX (private working copy only)",
     "private_local_only"),
    ("zhu2025_src/Figure_1.png",
     None,
     "Original figure from arXiv bundle (private working copy only)",
     "private_local_only"),
    ("zhu2025_src/Figure_2.png",
     None,
     "Original figure from arXiv bundle (private working copy only)",
     "private_local_only"),
    ("zhu2025_src/Figure_3.png",
     None,
     "Original figure from arXiv bundle (private working copy only)",
     "private_local_only"),
    ("zhu2025_src/00README.json",
     None,
     "arXiv bundle metadata (private working copy only)",
     "private_local_only"),
]


def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    entries = []
    for rel, url, role, availability in DECLARED:
        p = ARCHIVE_DIR / rel
        entries.append(
            {
                "path": rel,
                "url": url,
                "role": role,
                "availability": availability,
                "local_archive": "sibling of the repository root (outside public staging; "
                                 "absolute local path NOT published)",
                "size_bytes": p.stat().st_size if p.is_file() else None,
                "sha256": sha256(p),
            }
        )
    manifest = {
        "schema_version": "stage1_source_manifest_v2",
        "description": "Declared SHA-256 registry of Zhu et al. PRL 135, 056801 (2025) sources "
                       "used by the Stage-1 review. Copyrighted artifacts are private local "
                       "archives (hashed for provenance, never redistributed in this repository). "
                       "Only URLs, hashes, short protocol excerpts and archive notes are public.",
        "download_date": "2026-10-20",
        "note": "git-tracked content contains no publisher PDFs, TeX sources, original figures "
                "or full-text SM dumps; see sources/README.md.",
        "files": entries,
    }
    out = ROOT / "results" / "SOURCES_MANIFEST.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
