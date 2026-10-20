#!/usr/bin/env python3
"""Scoped GitHub core gate. Does NOT certify the offline v30 research archive."""
from pathlib import Path
import subprocess, sys, json
root=Path(__file__).resolve().parents[1]
scientific=json.loads((root/"results/frozen/core_reference.json").read_text())
assert scientific["scientific_status"]=="PARTIAL"
assert scientific["scope"]=="published-core-subset"
run=subprocess.run([sys.executable,"-m","unittest","discover","-s","tests","-v"],
                   cwd=root,text=True,capture_output=True)
print(run.stdout);print(run.stderr,file=sys.stderr)
if run.returncode:raise SystemExit(run.returncode)
from qgmap.dirac import eq11_kernel
if abs(float(eq11_kernel(1/8))-scientific["dirac_beta_1_8_kernel"])>1e-13:
    raise RuntimeError("Frozen Dirac benchmark drift")
print("GitHub core software checks: PASS; overall scientific status: PARTIAL")
