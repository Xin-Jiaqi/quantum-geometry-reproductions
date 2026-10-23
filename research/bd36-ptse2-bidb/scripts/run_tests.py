#!/usr/bin/env python3
"""BD36 R1.2 test runner: reproducible PASS/FAIL suite.

Sections:
  U1 read_spn        : --selftest (format roundtrip, Wannier transform,
                       Hermiticity) + bad-fixture rejections
  U2 make_win_kpoints: VASP official explicit example (Cartesian) rejected;
                       Gamma 81-point full-BZ parse; index-regression fixture
  U3 verify_relaxed   : --selftest (z-only PASS; in-plane drift FAIL)
  S1 static checker  : POSITIVE tree (built INPUTS/) -> all PASS
  S2 static checker  : NEGATIVE trees (each with ONE targeted corruption) ->
                       checker must FAIL each (independence: not self-validation)

Usage: python3 run_tests.py [--log LOGS/r12_test_run.log]
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = ROOT / "INPUTS"
LOGS = ROOT / "LOGS"
PY = sys.executable


def sh(args, expect_rc=0):
    a = list(args)
    if a[0].endswith(".py") and "/" not in a[0]:
        a[0] = str(HERE / a[0])   # resolve bare script names to scripts/
    r = subprocess.run([PY] + a, capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    return r.returncode, out


def section(title, results):
    print(f"\n=== {title} ===")
    bad = 0
    for name, rc_ok, rc, msg in results:
        state = "PASS" if rc_ok else "FAIL"
        if not rc_ok:
            bad += 1
        print(f"  [{state}] {name}")
        if not rc_ok and msg:
            print(f"          {msg.strip().splitlines()[-1]}")
    return bad


def neg_tree(name, mutate):
    """Copy INPUTS to a fixture tree and apply one corruption via mutate().
    Returns the tree path."""
    dst = LOGS / "test_fixtures" / name
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(INPUTS, dst)
    mutate(dst)
    return dst


def _mk_ibzkpt(name, lines):
    f = LOGS / name
    f.write_text("\n".join(lines) + "\n")
    return f


def _mk_contcar(name, dz_top, dz_bot, drift=None):
    """Build a simulated relaxed CONTCAR (with Selective dynamics header)
    from the 2253 baseline POSCAR: z-only shifts + optional in-plane drift."""
    import numpy as np
    sys.path.insert(0, str(HERE))
    from read_poscar_utils import parse_poscar_file, basic_poscar
    pos = INPUTS / "2253_bilayer_baseline" / "POSCAR"
    lat, sp, cnt, coords, fl = parse_poscar_file(pos)
    c = coords.copy().astype(float)
    zs = sorted({round(float(c[i][2]), 4) for i in range(len(sp)) if sp[i] == "Pt"})
    tol = 0.5 * (zs[1] - zs[0])
    for i in range(len(c)):
        if sp[i] == "Pt":
            c[i, 2] += dz_top if abs(c[i, 2] - zs[1]) < tol else dz_bot
        elif sp[i] == "Se":
            c[i, 2] += 0.5 * dz_top if abs(c[i, 2] - zs[1]) < tol else 0.5 * dz_bot
        if drift is not None and sp[i] == "Pt" and abs(c[i, 2] - zs[1]) < tol:
            c[i, 0] += drift  # in-plane registry-breaking drift (frac x)
    f = LOGS / name
    basic_poscar(f, lat, sp, cnt, c, sd_flags=["F F T"] * len(c))
    return f


def main():
    results = []

    # --- generate all regression fixtures here (self-contained) ---
    (LOGS / "test_kpoints_official_ibzkpt").write_text(
        "Explicit k-point list\n4\nCartesian\n"
        "0.0  0.0  0.0   1\n0.0  0.0  0.5   1\n0.0  0.5  0.5   2\n"
        "0.5  0.5  0.5   4\n")   # VASP wiki official explicit example
    k81 = ["Automatic mesh", "81", "Gamma"]
    for i in range(9):
        for j in range(9):
            k81.append(f"{i/9.0:.6f} {j/9.0:.6f} 0.000000 1.000000")
    _mk_ibzkpt("test_ibzkpt_81", k81)
    _mk_ibzkpt("test_ibzkpt_3", ["Automatic mesh", "3", "Gamma",
                                 "0.000000 0.000000 0.000000 1",
                                 "0.333333 0.000000 0.000000 1",
                                 "0.333333 0.333333 0.000000 1"])
    _mk_contcar("test_contcar_zonly", dz_top=0.02, dz_bot=-0.02)
    _mk_contcar("test_contcar_drift", dz_top=0.02, dz_bot=-0.02, drift=0.05)

    # U1 read_spn
    rc, out = sh(["read_spn.py", "--selftest"])
    results.append(("read_spn selftest (roundtrip+Wannier+Hermiticity)", rc == 0, rc, out))
    # bad .spn: header garbage
    bad = LOGS / "test_fixtures" / "bad_spn" / "x.spn"
    bad.parent.mkdir(parents=True, exist_ok=True)
    bad.write_text("Written ...\nnope\n0.1 0.2\n")
    rc, out = sh(["read_spn.py", str(bad)])
    results.append(("read_spn rejects bad header", rc != 0, rc, out))

    # binary (unformatted) .spn parse: 4-band sigma_y via official packing
    sys.path.insert(0, str(HERE))
    import numpy as np
    from read_spn import (official_binary_spn, official_text_spn,
                          SX4, SY4, SZ4, parse_spn)
    nbnd4 = 4
    mats4 = np.zeros((1, nbnd4, nbnd4, 3), dtype=complex)
    mats4[0, :, :, 0] = SX4
    mats4[0, :, :, 1] = SY4
    mats4[0, :, :, 2] = SZ4
    (LOGS / "r13_spn_4band_binary.spn").write_bytes(official_binary_spn(mats4))
    (LOGS / "r13_spn_4band_official.spn").write_text(official_text_spn(mats4))
    rc, out = sh(["read_spn.py", str(LOGS / "r13_spn_4band_binary.spn"), "1"])
    ok = rc == 0 and "hermitian" in out
    results.append(("read_spn parses BINARY .spn (auto-detect)", ok, rc, out))
    rc, out = sh(["read_spn.py", str(LOGS / "r13_spn_4band_official.spn"), "1"])
    ok = rc == 0 and "hermitian" in out
    results.append(("read_spn parses FORMATTED 4-band official", ok, rc, out))

    # U2 make_win_kpoints
    rc, out = sh(["make_win_kpoints.py", str(LOGS / "test_kpoints_official_ibzkpt")])
    results.append(("win-kpts rejects VASP Cartesian example", rc != 0, rc, out))
    rc, out = sh(["make_win_kpoints.py", str(LOGS / "test_ibzkpt_81"), "9 9 1"])
    ok = rc == 0 and "count=81" in out and "weights sum=81" in out
    results.append(("win-kpts parses 81-pt Gamma full-BZ (count/weights)", ok, rc, out))
    rc, out = sh(["make_win_kpoints.py", str(LOGS / "test_ibzkpt_3")])
    ok = rc == 0 and out.splitlines()[1].startswith(" 0.000000")  # Gamma first!
    results.append(("win-kpts index regression (Gamma NOT skipped)", ok, rc, out))
    # truncated IBZKPT: 2 k-lines declared 3
    t = LOGS / "test_fixtures" / "bad_ibzkpt"
    t.mkdir(parents=True, exist_ok=True)
    (t / "IBZKPT").write_text("Automatic mesh\n3\nGamma\n0 0 0 1\n0 0.3 0 1\n")
    rc, out = sh(["make_win_kpoints.py", str(t / "IBZKPT")])
    results.append(("win-kpts rejects k-count mismatch", rc != 0, rc, out))

    # U3 verify_relaxed
    rc, out = sh(["verify_relaxed.py", "--selftest"])
    results.append(("verify_relaxed selftest (z-only PASS; drift FAIL)", rc == 0, rc, out))
    # SD-CONTCAR parse regression: fixture A has 'Selective dynamics' header
    rc, out = sh(["verify_relaxed.py", str(INPUTS / "2253_bilayer_baseline" / "POSCAR"),
                  str(LOGS / "test_contcar_zonly"), "0", "0"])
    results.append(("verify_relaxed parses real SD CONTCAR (z-only PASS)",
                    rc == 0, rc, out))
    rc, out = sh(["verify_relaxed.py", str(INPUTS / "2253_bilayer_baseline" / "POSCAR"),
                  str(LOGS / "test_contcar_drift"), "0", "0"])
    results.append(("verify_relaxed FAILs in-plane drift (0.187 A > 1e-3)",
                    rc != 0, rc, out))

    # S1 positive tree
    rc, out = sh(["static_check_inputs.py"])
    ok = rc == 0 and "== BUDGET Phase0 56 <= 56; Phase1 448 <= 448; total 504 <= 504 PASS" in out
    results.append(("static checker on built INPUTS (108 PASS, phased budget)", ok, rc, out))

    # S2 negative trees
    def mut_sd(dst):
        f = dst / "2253_bilayer_baseline" / "POSCAR_SD"
        lines = f.read_text().splitlines()
        lines = [l.replace("F F T", "T F F") if l.split()[-3:] == ["F", "F", "T"] else l
                 for l in lines]
        f.write_text("\n".join(lines) + "\n")
    rc, out = sh(["static_check_inputs.py", "--inputs", str(neg_tree("neg_sd_TFF", mut_sd))])
    results.append(("neg: SD flags T F F (R1.1 bug) detected", rc != 0, rc, out))

    def mut_incar(dst):
        f = dst / "2253_bilayer_baseline" / "INCAR_soc"
        s = f.read_text().replace("LASPH     = .TRUE.", "")
        f.write_text(s)
    rc, out = sh(["static_check_inputs.py", "--inputs", str(neg_tree("neg_no_lasph", mut_incar))])
    results.append(("neg: missing LASPH detected", rc != 0, rc, out))

    def mut_win(dst):
        f = dst / "2253_bilayer_baseline" / "wannier90.win"
        f.write_text(f.read_text() + "\nspin = true\n")
    rc, out = sh(["static_check_inputs.py", "--inputs", str(neg_tree("neg_win_spin", mut_win))])
    results.append(("neg: win spin=true detected", rc != 0, rc, out))

    def mut_budget(dst):
        caps = json.load(open(dst / "caps.json"))
        caps["total_core_h"] = 500   # old '<=500' claim, now invalid
        json.dump(caps, open(dst / "caps.json", "w"), indent=2)
    rc, out = sh(["static_check_inputs.py", "--inputs", str(neg_tree("neg_budget500", mut_budget))])
    results.append(("neg: total-cap 500 (unphased) detected", rc != 0, rc, out))

    def mut_isym(dst):
        f = dst / "2253_bilayer_baseline" / "INCAR_soc"
        f.write_text(f.read_text().replace("ISYM     = -1\n", ""))
    rc, out = sh(["static_check_inputs.py", "--inputs", str(neg_tree("neg_no_isym", mut_isym))])
    results.append(("neg: INCAR_soc missing ISYM=-1 detected", rc != 0, rc, out))

    def mut_lspn(dst):
        f = dst / "2253_bilayer_baseline" / "INCAR_soc_w90_vasp544"
        f.write_text(f.read_text().replace("LSPN_FMTED  = .TRUE.\n", ""))
    rc, out = sh(["static_check_inputs.py", "--inputs", str(neg_tree("neg_no_lspn", mut_lspn))])
    results.append(("neg: vasp544 missing LSPN_FMTED detected", rc != 0, rc, out))

    # phase0 script: bash syntax + hygiene (no internal absolute paths)
    r = subprocess.run(["bash", "-n",
                        str(INPUTS / "2253_bilayer_baseline" / "phase0_step0.sbatch")],
                       capture_output=True, text=True)
    results.append(("phase0_step0.sbatch bash syntax OK", r.returncode == 0,
                    r.returncode, r.stderr))
    txt = (INPUTS / "2253_bilayer_baseline" / "phase0_step0.sbatch").read_text()
    import re as _re
    rc = 1 if _re.search(r"/(home|share|Users|work|scratch)/", txt) else 0
    results.append(("phase0_step0.sbatch hygiene (no internal paths)", rc == 0, rc, ""))

    def mut_gga(dst):
        f = dst / "2253_bilayer_baseline" / "INCAR_relax"
        f.write_text(f.read_text().replace("GGA_COMPAT = .FALSE.", "GGA_COMPAT = .TRUE."))
    rc, out = sh(["static_check_inputs.py", "--inputs", str(neg_tree("neg_gga_compat", mut_gga))])
    results.append(("neg: GGA_COMPAT=.TRUE. detected", rc != 0, rc, out))

    bad = section("U1/U2/U3 unit + regression", results[:12])
    bad += section("S1/S2 checker positive + negative", results[12:])
    print(f"\nTOTAL: {len(results) - bad}/{len(results)} PASS")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
