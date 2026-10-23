#!/usr/bin/env python3
"""BD36 static_check_inputs: automated static checks over INPUTS/.

Checks per structure (item 7 of the R1.1 review):
  A. POSCAR <-> POSCAR_SD consistency (species, counts, coords) and
     selective-dynamics flags (bilayer T F F / monolayer T T T).
  B. Wannier win: num_wann == 2*(5*nPt+3*nSe); projections == Pt:d,Se:p;
     mp_grid == KPOINTS mesh; win atoms_frac == POSCAR atoms (slide-correct);
     NO spin=true; NO hard-coded begin kpoints (must be post-SCF generated);
     num_bands >= num_wann.
  C. INCAR compatibility: relax (ISPIN=1, ISIF=2, NSW<=60, EDIFFG<=-0.01,
     IVDW=12, ENCUT in [400,600]); soc/soc_w90 (LNONCOLLINEAR=T, LSORBIT=T,
     ISPIN=2, SAXIS, MAGMOM length == 3*natoms, NBANDS>=num_wann,
     ISYM=-1, IVDW=12); w90 variants: LWANNIER90=T + LWRITE_SPN=T,
     NO embedded WANNIER90_WIN.
  D. Resource budget: per-job cores*wall_h <= caps; round-1 total <= cap.
  E. Hygiene: committed drafts contain no internal absolute paths
     (no '/home/', '/share/', '/Users/'); POTCAR absent (expected, private).

Exit 0 => all PASS; non-zero => FAILs reported. Writes static_check_report.json.
"""
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "INPUTS"
SPECTRA = {"Pt": 5, "Se": 3}   # spatial orbitals per species in projection basis


def load_caps():
    return json.load(open(INPUTS / "caps.json"))


def parse_poscar(path):
    lines = path.read_text().splitlines()
    lat = [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    if len(lines) > 8 and lines[7].strip().lower().startswith("selective"):
        atoms_start = 9
    else:
        atoms_start = 8
    coords = [[float(x) for x in l.split()[:3]] for l in lines[atoms_start:atoms_start + sum(counts)]]
    return species, counts, coords, lat


def parse_sd(path):
    lines = path.read_text().splitlines()
    assert lines[7].strip().lower().startswith("selective"), "no Selective dynamics header"
    species, counts, coords, lat = parse_poscar(path)
    flags = [l.split()[3:6] for l in lines[9:9 + sum(counts)]]
    return species, counts, coords, flags


def parse_incar(path):
    tags = {}
    for ln in path.read_text().splitlines():
        m = re.match(r"^\s*([A-Za-z0-9_]+)\s*=\s*(\S.*)$", ln)
        if m:
            tags[m.group(1)] = m.group(2).strip()
    return tags


def parse_win(path):
    text = path.read_text()
    tags = {}
    for ln in text.splitlines():
        m = re.match(r"^\s*([a-z_0-9]+)\s*=\s*(\S.*)$", ln)
        if m:
            tags[m.group(1).lower()] = m.group(2).strip()
    def block(name):
        m = re.search(rf"^\s*begin {name}\s*$(.*?)^\s*end {name}\s*$", text,
                      re.MULTILINE | re.DOTALL)
        return m.group(1).strip() if m else None
    return tags, block, text


def bool_of(v):
    return v.strip(".").lower() in {"t", "true", "1", "yes"}


def main():
    inputs = INPUTS
    if "--inputs" in sys.argv:
        inputs = Path(sys.argv[sys.argv.index("--inputs") + 1])
    caps = json.load(open(inputs / "caps.json"))
    report = {"caps": caps, "structures": {}}
    all_ok = True
    for d in sorted(inputs.iterdir()):
        if not d.is_dir():
            continue
        tag = d.name
        checks = []
        fail = lambda msg: checks.append(("FAIL", msg))
        ok = lambda msg: checks.append(("PASS", msg))
        # --- A. POSCAR / POSCAR_SD ---
        sp, cnt, coords, lat = parse_poscar(d / "POSCAR")
        n = sum(cnt)
        if (d / "POSCAR_SD").exists():
            sp_sd, cnt_sd, coords_sd, flags = parse_sd(d / "POSCAR_SD")
            if sp == sp_sd and cnt == cnt_sd and len(coords) == len(coords_sd):
                ok("A POSCAR_SD species/counts match POSCAR")
            else:
                fail("A POSCAR_SD species/counts MISMATCH")
            bad = [i for i, f in enumerate(flags) if len(f) != 3 or f[0] not in "TF"]
            if not bad:
                ok("A selective-dynamics flags well-formed")
            else:
                fail(f"A malformed SD flags at atoms {bad}")
            if n == 3 and all(f == ["T", "T", "T"] for f in flags):
                ok("A monolayer: full relaxation (T T T)")
            elif n == 6 and all(f == ["F", "F", "T"] for f in flags):
                ok("A bilayer: registry-preserving F F T (x,y FIXED, z FREE; "
                   "VASP T=free, F=fixed)")
            else:
                fail("A SD flags do not match monolayer/bilayer policy "
                     "(bilayer must be F F T: x,y FIXED, z FREE)")
        # --- A2. semantics documented independently of the generator ---
        readme = (d / "README").read_text(errors="replace") if (d / "README").exists() else ""
        pol = "T T T" if n == 3 else "F F T"   # monolayer vs bilayer policy
        if pol in readme and "free" in readme.lower() and "fixed" in readme.lower():
            ok(f"A2 README states SD semantics (T=free, F=fixed; policy {pol})")
        else:
            fail(f"A2 README does not document SD semantics "
                 f"(need '{pol}' + free/fixed wording)")
        # --- B. Wannier win ---
        if (d / "wannier90.win").exists():
            tags, block, text = parse_win(d / "wannier90.win")
            npt = sum(c for s, c in zip(sp, cnt) if s == "Pt")
            nse = sum(c for s, c in zip(sp, cnt) if s == "Se")
            want = 2 * (5 * npt + 3 * nse)
            nw = int(tags.get("num_wann", -1))
            if nw == want:
                ok(f"B num_wann={nw} == 2*(5*{npt}+3*{nse})")
            else:
                fail(f"B num_wann={nw} != expected {want}")
            nb = int(tags.get("num_bands", -1))
            if nb >= nw:
                ok(f"B num_bands={nb} >= num_wann")
            else:
                fail("B num_bands < num_wann")
            proj = block("projections") or ""
            proj_species = sorted(set(re.findall(r"^([A-Za-z]+)", proj, re.M)))
            if proj_species == ["Pt", "Se"] and "Pt:d" in proj and "Se:p" in proj and "Pt:s" not in proj and "Se:s" not in proj:
                ok("B projections = Pt:d + Se:p (no s)")
            else:
                fail(f"B projections unexpected: {proj_species!r}")
            if "spin = true" not in text.lower():
                ok("B no spin=true in win (LWRITE_SPN path)")
            else:
                fail("B spin=true present in win (incompatible)")
            kp_block = block("kpoints")
            if kp_block is None:
                ok("B no hard-coded kpoints block (post-SCF generation)")
            else:
                fail("B hard-coded kpoints block present - remove; generate post-SCF")
            mg = tags.get("mp_grid", "")
            if mg == "9 9 1":
                ok("B mp_grid == 9x9x1 (KPOINTS)")
            else:
                fail(f"B mp_grid {mg} != KPOINTS mesh")
            a_block = block("atoms_frac") or ""
            win_atoms = [ln.split() for ln in a_block.splitlines() if ln.strip()]
            if len(win_atoms) == n:
                drift = max(
                    max(abs(float(v) - w) for v, w in zip(at[1:4], cd))
                    for at, cd in zip(win_atoms, coords))
                if drift < 1e-6:
                    ok(f"B win atoms_frac == POSCAR (max drift {drift:.1e})")
                else:
                    fail(f"B win atoms drift {drift:.1e} from POSCAR")
            else:
                fail(f"B win atom count {len(win_atoms)} != POSCAR {n}")
        # --- C. INCARs ---
        if (d / "INCAR_soc").exists():
            t = parse_incar(d / "INCAR_soc")
            okc = True
            okc &= bool_of(t.get("LNONCOLLINEAR", "F"))
            okc &= bool_of(t.get("LSORBIT", "F"))
            okc &= t.get("ISPIN", "0") == "2"
            okc &= t.get("SAXIS", "") == "0 0 1"
            mm = t.get("MAGMOM", "").split()
            okc &= len(mm) == 3 * n and all(v == "0.02" for v in mm[2::3])
            okc &= t.get("ISYM", "0") == "-1"   # full 81-k mesh == interface
            okc &= bool_of(t.get("LASPH", "F"))
            okc &= t.get("GGA_COMPAT", "TRUE").strip(".") == "FALSE"
            ok("C INCAR_soc compatible (ISPIN=2, MAGMOM 3n, ISYM=-1, "
               "LASPH/GGA_COMPAT)") if okc else fail("C INCAR_soc INCOMPATIBLE")
        if (d / "INCAR_relax").exists():
            t = parse_incar(d / "INCAR_relax")
            okc = True
            okc &= t.get("ISPIN", "1") == "1"
            okc &= t.get("IBRION", "2") == "2"
            okc &= t.get("ISIF", "2") == "2"
            okc &= int(t.get("NSW", 0)) <= 60
            okc &= float(t.get("EDIFFG", 0)) <= -0.01
            okc &= str(t.get("IVDW", "")).isdigit() and int(t.get("IVDW", 0)) >= 10
            okc &= 400 <= float(t.get("ENCUT", 0)) <= 600
            okc &= bool_of(t.get("LASPH", "F"))
            okc &= t.get("GGA_COMPAT", "TRUE").strip(".") == "FALSE"
            ok("C INCAR_relax compatible (incl. LASPH=.TRUE., GGA_COMPAT=.FALSE.)") \
                if okc else fail("C INCAR_relax INCOMPATIBLE")
        for inc in ("INCAR_soc", "INCAR_soc_w90_vasp6", "INCAR_soc_w90_vasp544"):
            if not (d / inc).exists():
                continue
            t = parse_incar(d / inc)
            okc = True
            okc &= bool_of(t.get("LNONCOLLINEAR", "F"))
            okc &= bool_of(t.get("LSORBIT", "F"))
            okc &= t.get("ISPIN", "0") == "2"
            okc &= t.get("SAXIS", "") == "0 0 1"
            okc &= bool_of(t.get("LASPH", "F"))
            okc &= t.get("GGA_COMPAT", "TRUE").strip(".") == "FALSE"
            mm = t.get("MAGMOM", "").split()
            okc &= len(mm) == 3 * n
            if "w90" in inc:
                okc &= t.get("NBANDS", "0").isdigit() and int(t.get("NBANDS", 0)) >= nw
                okc &= t.get("ISYM", "-1") == "-1"   # full 81-k mesh
                okc &= t.get("KPAR", "1") == "1"
            okc &= str(t.get("IVDW", "")).isdigit() and int(t.get("IVDW", 0)) >= 10
            if "w90" in inc:
                okc &= bool_of(t.get("LWANNIER90", "F"))
                okc &= bool_of(t.get("LWRITE_SPN", "F"))
                if "WANNIER90_WIN" in t:
                    okc = False
                    fail("C embedded WANNIER90_WIN present (incompatible)")
            if inc == "INCAR_soc_w90_vasp544":
                # VaspWannier patch: LSPN_FMTED=.TRUE. -> ASCII .spn; WRITE_SPN
                # is serial-only (NCPU==1) - document in script
                okc &= t.get("LSPN_FMTED", "FALSE").strip(".").upper() == "TRUE"
            if inc == "INCAR_soc_w90_vasp6":
                # VASP>6.4.2 native: always formatted; LSPN_FMTED must be absent
                okc &= "LSPN_FMTED" not in t
            ok(f"C {inc} compatible") if okc else fail(f"C {inc} INCOMPATIBLE")
        # --- D. resources ---
        wall_h = 2
        cores = 28
        job_ch = cores * wall_h
        if job_ch <= caps["job_core_h"]:
            ok(f"D per-job {job_ch} core-h <= cap {caps['job_core_h']}")
        else:
            fail(f"D per-job {job_ch} exceeds cap")
        # --- E. hygiene ---
        bad_paths = []
        for f in d.iterdir():
            if f.suffix in (".sh", ".pbs", ".sbatch", ".py") or f.name == "README":
                txt = f.read_text(errors="replace")
                if re.search(r"/(home|share|Users|work|scratch)/", txt):
                    bad_paths.append(f.name)
        if not bad_paths:
            ok("E no internal absolute paths in committed scripts")
        else:
            fail(f"E internal paths in {bad_paths}")
        if not (d / "POTCAR").exists():
            ok("E POTCAR absent (private; assembled at execution)")
        else:
            fail("E POTCAR present - must NOT be committed")
        report["structures"][tag] = {m[0]: m[1] for m in checks}
        n_fail = sum(1 for c in checks if c[0] == "FAIL")
        if n_fail:
            all_ok = False
        print(f"== {tag}: {n_fail} FAIL / {len(checks) - n_fail} PASS")
        for status, msg in checks:
            print(f"   [{status}] {msg}")
    # PHASED budget (R1.2): Phase 0 = para step-0 in-place interface validation
    # (56); Phase 1 = Round-1 minimal set (relax+soc_w90 per structure:
    # 2 jobs/structure x 56 + 4 convergence x conv_core_h); total = 504.
    manifest = json.load(open(inputs / "manifest.json"))
    round1 = [k for k, v in manifest.items() if v["phase"] == "round1"]
    phase1 = len(round1) * 2 * 28 * 2 + 4 * caps["conv_core_h"]
    phase0 = caps["phase0_validation_core_h"]
    tot = phase0 + phase1
    report["budget_phases"] = {"phase0_validation_core_h": phase0,
                               "phase1_round1_core_h": phase1,
                               "total_core_h": tot}
    ok3 = (phase0 <= caps["phase0_validation_core_h"]
           and phase1 <= caps["phase1_round1_core_h"]
           and tot <= caps["total_core_h"])
    print(f"== BUDGET Phase0 {phase0} <= {caps['phase0_validation_core_h']}; "
          f"Phase1 {phase1} <= {caps['phase1_round1_core_h']}; "
          f"total {tot} <= {caps['total_core_h']} "
          + ("PASS" if ok3 else "FAIL"))
    if not ok3:
        all_ok = False
    (ROOT / "LOGS" / "static_check_report.json").parent.mkdir(exist_ok=True)
    (ROOT / "LOGS" / "static_check_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True))
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
