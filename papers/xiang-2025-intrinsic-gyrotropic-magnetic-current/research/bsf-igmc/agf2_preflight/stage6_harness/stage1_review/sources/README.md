# Stage-1 source archive policy

The copyrighted source files (publisher/arXiv PDFs, arXiv TeX bundle, original
figures, full-text SM extraction) are a **private local archive** stored
**OUTSIDE the repository tree** (sibling of the repo root inside the session
workspace, directory `private-agf2-stage1-archives/`; the absolute path is not
published). They must never be committed to the public repository:

- `zhu2025_prl_published.pdf`, `zhu2025_prl_sm.pdf`, `zhu2025_arxiv_v2.pdf` — publisher/arXiv PDFs;
- `zhu2025_arxiv_v2.tar.gz` — arXiv TeX source bundle;
- `zhu2025_prl_sm.txt` — full-text extraction of the Supplemental Material;
- `zhu2025_src/` — extracted TeX and original figures from the arXiv bundle.

This `sources/` directory in the repository keeps only this policy note
(`README.md`) and a `.gitignore` that ignores everything else, so the public
staging area can never accidentally carry a copyrighted attachment (verified
by `scripts/verify_round9.py` of the IGMC validated-method gate).

Only the following are tracked in git:
- `SOURCES.md` / `results/SOURCES_MANIFEST.json` (URLs, SHA-256, short verbatim
  protocol excerpts, and this archive policy);
- short protocol excerpts inside `results/agf2_protocol_frozen.json`.

Rationale: arXiv/publisher downloads are not a redistribution license. See the
review on PR #8 (item 1) and `SOURCES.md` in this review directory.
