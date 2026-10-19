#!/usr/bin/env python3
"""Fail-closed reproducibility gate for the released Jia et al. NMEE model.

PASS here means model/software reproduction, not proof of the paper's
CuMnAs finite-T PT-degenerate driven-dc intrinsic coefficient.
No private datasets, publisher PDF, network or GPU required.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "fig2_km_451.json": ["--case", "fig2-printed-vs-static", "--mesh", "451", "--temperature", "10", "--mu", "0.025"],
    "fig1b_cumnas.json": ["--case", "fig1b", "--mesh", "101", "--temperature", "100"],
    "fig1c_static_diagnostic.json": ["--case", "fig1c-static-diagnostic", "--mesh", "81", "--temperature", "100", "--mu", "0.03", "--theta", "0"],
    "fig1d_cumnas_eq2.json": ["--case", "fig1d-printed", "--mesh", "151", "--temperature", "100", "--mu", "0.03", "--theta", "0.15", "--tau", "50", "--a-nm", "1"],
    "fig1_cumnas_BGK.json": ["--case", "fig1-bgk", "--mesh", "81", "--temperature", "100", "--mu", "0.03", "--theta", "0.15"],
}

def compare(expected, actual, path="root"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(expected) != set(actual):
            raise AssertionError(f"JSON keys differ at {path}")
        for key in expected:
            compare(expected[key], actual[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"JSON length differs at {path}")
        for i, (x,y) in enumerate(zip(expected,actual)):
            compare(x, y, f"{path}[{i}]")
    elif type(expected) in (float, int) and type(actual) in (float, int):
        if not math.isclose(expected, actual, abs_tol=1e-10, rel_tol=2e-7):
            raise AssertionError(f"Numerical drift at {path}: expected={expected}, actual={actual}")
    elif expected != actual or type(expected) is not type(actual):
        raise AssertionError(f"Value changed at {path}: {expected!r} -> {actual!r}")

def verify_manifest():
    manifest = json.loads((ROOT / "results/RELEASE_MANIFEST.json").read_text(encoding="utf-8"))
    entries=manifest["source_sha256"]
    if not entries:
        raise AssertionError("Empty source SHA256 inventory")
    for rel,digest in entries.items():
        path=ROOT / rel
        if not path.is_file():
            raise AssertionError(f"Missing frozen artifact: {rel}")
        if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            raise AssertionError(f"SHA256 mismatch: {rel}")
    if manifest.get("scientific_acceptance") != "PARTIAL_PENDING_DC_DEFINITION":
        raise AssertionError("The unresolved physical response definition must not be hidden")
    return len(entries)

def verify_cases():
    for filename,opts in CASES.items():
        frozen=json.loads((ROOT/'results'/filename).read_text(encoding='utf-8'))
        proc=subprocess.run([sys.executable,str(ROOT/'scripts/reproduce.py'),*opts],
                            cwd=ROOT,capture_output=True,text=True,check=True,timeout=180)
        compare(frozen,json.loads(proc.stdout),filename)
    return len(CASES)

def run_tests():
    proc=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-v'],
                        cwd=ROOT,capture_output=True,text=True,timeout=180)
    if proc.returncode:
        print(proc.stdout, file=sys.stderr)
        print(proc.stderr, file=sys.stderr)
        raise AssertionError("Unit-test suite failed")
    details=proc.stderr
    line=next((v for v in details.splitlines() if v.startswith('Ran ')),None)
    print(line or 'Unit tests passed')
    return line

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-tests',action='store_true',help='Validate only frozen data/manifest')
    options=parser.parse_args()
    count=verify_manifest()
    print(f'Integrity: PASS ({count} SHA256 paths)')
    cases=verify_cases()
    print(f'Frozen model cases: PASS ({cases} outputs)')
    if not options.skip_tests:run_tests()
    print('Software/model reproduction: PASS')
    print('Scientific acceptance: PARTIAL (CuMnAs finite-T PT-degenerate driven-DC intrinsic definition awaits clarification)')

if __name__=='__main__':main()
