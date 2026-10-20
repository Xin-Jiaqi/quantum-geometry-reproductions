#!/usr/bin/env python3
"""Release-candidate integrity gate: offline, refuses a full-paper PASS claim.

Run from any cwd: python scripts/verify_round9.py.
Calls the frozen Round1–8 scientific verifier; this does not check GitHub Actions.
"""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_STATES = {
    "scientific_status": "PARTIAL",
    "publication_status": "LOCAL_DRAFT_NOT_PUBLISHED",
    "allow_empirical_multiplier": False,
    "publisher_pdf_or_figures_redistributed": False,
}
REQUIRED = [
    "results/RELEASE_MANIFEST.json",
    "results/frozen/round8_green_wavepacket_reference_v1.json",
    "data/toolkit_mpg_generators_v018.txt",
    "data/TOOLKIT_BSD3_LICENSE.txt",
    "docs/FINAL_ACCEPTANCE_MATRIX.md",
    "docs/ROUND8_SEMICLASSICAL_GREEN_CLOSURE.md",
    "results/raw/round7_metal_audit.json",
]


def preflight():
    cfg = json.loads((ROOT / "release/release_contract.json").read_text(encoding="utf-8"))
    for key, expected in EXPECTED_STATES.items():
        if cfg.get(key) != expected:
            raise AssertionError(f"Cannot promote this release gate: {key}: {cfg.get(key)!r}")
    if cfg.get("scope", {}).get("original_fig2d_absolute") != "BLOCKED_AUTHOR_DATA":
        raise AssertionError("Missing author Fig.2(d) data must stay BLOCKED")
    if cfg.get("scope", {}).get("eq14_antiunitary") != "PARTIAL_SOURCE_CONVENTION":
        raise AssertionError("Eq.14 source convention must remain PARTIAL")
    for name in REQUIRED:
        if not (ROOT/name).is_file():
            raise AssertionError(f"Missing required evidence: {name}")
    forbidden=[]
    for file in ROOT.rglob("*"):
        if not file.is_file():
            continue
        if file.suffix.lower() in {".pdf", ".pptx", ".docx", ".zip", ".tar", ".gz"}:
            forbidden.append(str(file.relative_to(ROOT)))
        if file.name in {'USER_PROTOCOL_SNAPSHOT.md'}:
            forbidden.append(str(file.relative_to(ROOT)))
    if forbidden:
        raise AssertionError(f"Source/proprietary attachment in public staging: {forbidden}")
    # Immutable release-candidate receipt is a separate file, not part of the
    # prior frozen model. It protects this round's status contract and Round8
    # reference against accidental edits that old acceptance scripts cannot see.
    receipt = json.loads((ROOT / "release/round9_receipt.json").read_text(encoding="utf-8"))
    if receipt.get("status") != "LOCAL_DRAFT_PARTIAL":
        raise AssertionError("Release receipt unexpectedly promoted")
    for name, expected in receipt["files_sha256"].items():
        found = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if found != expected:
            raise AssertionError(f"Round9 release receipt hash mismatch: {name}")
    found8 = hashlib.sha256((ROOT / cfg["frozen_baseline_round8"]).read_bytes()).hexdigest()
    if found8 != receipt["round8_frozen_reference_sha256"]:
        raise AssertionError("Round8 reference changed after Round9 receipt")
    workflow = ROOT.parents[1] / '.github/workflows/igmc.yml'
    if workflow.is_file() and hashlib.sha256(workflow.read_bytes()).hexdigest() != receipt["staged_workflow_sha256"]:
        raise AssertionError("Staged IGMC workflow changed after Round9 receipt")
    frozen = ROOT / cfg["frozen_baseline_round8"]
    if frozen.is_symlink():
        raise AssertionError("Frozen checkpoint must not be a symlink")
    return cfg


def main():
    c=preflight()
    run=subprocess.run([sys.executable,'scripts/verify_release.py'],cwd=ROOT,
        env=os.environ|{'PYTHONPATH':str(ROOT/'src')},text=True,capture_output=True)
    if run.returncode:
        sys.stderr.write(run.stdout[-10000:]+"\n"+run.stderr[-10000:])
        raise SystemExit(run.returncode)
    if 'PASS_MODEL_OFFLINE' not in run.stdout:
        raise AssertionError('Historical frozen model gate did not identify a passing subset')
    result={
        'release_gate':'PASS_LOCAL_DRAFT',
        'scientific_status':c['scientific_status'],
        'publication_status':c['publication_status'],
        'baseline':'Round1-8 frozen model checks and unittest suite',
        'independent_full_paper_acceptance':False,
        'github_ci':'NOT_RUN',
    }
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
