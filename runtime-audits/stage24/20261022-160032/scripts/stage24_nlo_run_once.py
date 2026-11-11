#!/usr/bin/env python3
"""Stage 24 — ONE low-grid `nlo-response run --response shift_current` on an isolated package.

Per the PR #66 review: run exactly once, on an isolated temporary package, with
spin_degeneracy=1 EXPLICITLY frozen (the NLO planner would otherwise auto-assign
spin_degeneracy=2 for soc:false, which is wrong for a genuinely spinless model),
and collect raw_band_sum, spectrum.csv, run-manifest.json and the console log.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    work = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(".").resolve()
    nlo_root = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else None
    bridge = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else None
    if nlo_root is None or bridge is None:
        raise SystemExit("usage: stage24_nlo_run_once.py WORK NLO_ROOT BRIDGE_SCRIPT")
    models = work / "models"
    models.mkdir(parents=True, exist_ok=True)
    package = models / "stage24_haldane_synthetic"
    if package.exists():
        shutil.rmtree(package)
    mapping_json = work / "STAGE24_HALDANE_PACKAGE.json"
    build = subprocess.run(
        [sys.executable, str(bridge), "--package-dir", str(package),
         "--report", str(mapping_json)],
        cwd=work, text=True, capture_output=True,
    )
    build_log = work / "logs" / "stage24_package_build.log"
    build_log.parent.mkdir(parents=True, exist_ok=True)
    build_log.write_text(build.stdout + build.stderr, encoding="utf-8")
    if build.returncode != 0:
        print(json.dumps({"step": "package_build", "status": "FAIL",
                          "stderr": build.stderr[-2000:]}))
        return 1

    output_root = work / "nlo-run-output"
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True)
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONPATH"] = str(nlo_root / "src")
    code = "from nlo_response.cli import main; raise SystemExit(main())"
    command = [
        sys.executable, "-c", code, "run", str(package), "--structure", "L1",
        "--response", "shift_current", "--output-root", str(output_root),
        "--dimensionality", "2D", "--kmesh", "5x5x1",
        "--energy-min", "0", "--energy-max", "2", "--energy-points", "401",
        "--broadening", "0.05", "--spin-degeneracy", "1",
        "--sign-factor", "-1", "--no-mat", "--workers", "1",
    ]
    start = datetime.now(timezone.utc)
    run = subprocess.run(command, cwd=nlo_root, env=env, text=True,
                         capture_output=True, timeout=1800)
    elapsed_s = round((datetime.now(timezone.utc) - start).total_seconds(), 3)
    out_log = work / "logs" / "stage24_nlo_run_once.stdout.log"
    err_log = work / "logs" / "stage24_nlo_run_once.stderr.log"
    out_log.write_text(run.stdout, encoding="utf-8")
    err_log.write_text(run.stderr, encoding="utf-8")

    batch_dirs = list(output_root.glob("shift_current_*"))
    manifest = results = spectrum = None
    if batch_dirs:
        batch_dir = batch_dirs[0]
        for cand in batch_dir.iterdir():
            if cand.name == "run-manifest.json":
                manifest = json.loads(cand.read_text(encoding="utf-8"))
            elif cand.name == "results.json":
                results = json.loads(cand.read_text(encoding="utf-8"))
            elif cand.name == "spectrum.csv":
                spectrum = cand.read_text(encoding="utf-8")
    summary = {
        "schema": "stage24-nlo-run-once/v1",
        "created_at_utc": start.isoformat(),
        "elapsed_s": elapsed_s,
        "returncode": run.returncode,
        "status": "PASS" if run.returncode == 0 else "FAIL",
        "command": " ".join(command),
        "package_dir": str(package),
        "output_root": str(output_root),
        "spin_degeneracy_frozen": 1,
        "auto_spin_degeneracy_warning": (
            "NLO planning.py auto-assigns spin_degeneracy=2 when soc:false and no "
            "explicit value is given; this run passes --spin-degeneracy 1 explicitly "
            "because the synthetic Haldane model is genuinely spinless."
        ),
        "run_manifest": manifest,
        "results_summary": (
            {
                "raw_band_sum_shape": None if results is None else
                [len(results.get("raw_band_sum", [])),
                 len(results.get("raw_band_sum", [[]])[0]) if results.get("raw_band_sum") else 0],
                "tensor_codes": results.get("tensor_codes"),
                "photon_energy_eV": results.get("photon_energy_eV"),
                "response_cell_unit": (
                    results.get("response", {}).get("cell", {}).get("unit")
                    if results else None),
                "metadata_consistency": results.get("metadata", {}).get("consistency"),
                "metadata_gap": results.get("metadata", {}).get("gap"),
            }
            if results else None),
        "spectrum_csv_bytes": len(spectrum or ""),
        "spectrum_csv_tail": (spectrum or "").splitlines()[-4:],
        "stdout_tail": run.stdout[-1500:],
        "stderr_tail": run.stderr[-1500:],
    }
    summary_path = work / "STAGE24_NLO_RUN_ONCE.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
                            encoding="utf-8")
    print(json.dumps({"status": summary["status"], "returncode": run.returncode,
                      "batch_dir": str(batch_dirs[0]) if batch_dirs else None,
                      "elapsed_s": elapsed_s}, indent=2))
    return 0 if run.returncode == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
