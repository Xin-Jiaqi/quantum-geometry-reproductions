#!/usr/bin/env python3
"""read_poscar_utils: minimal POSCAR/CONTCAR read/write helpers shared by the
BD36 verification scripts (numpy-only, no external dependencies)."""

import numpy as np


def parse_poscar_file(path):
    """Returns (lat, species_expanded, counts, coords_frac, flags)."""
    lines = path.read_text().splitlines()
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)])
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    n = sum(counts)
    sd = lines[7].strip().lower().startswith("selective")
    start = 9 if sd else 8
    coords, flags = [], []
    for l in lines[start:start + n]:
        p = l.split()
        coords.append([float(x) for x in p[0:3]])
        if sd and len(p) >= 6:
            flags.append(p[3:6])
    expanded = [sp for sp, cnt in zip(species, counts) for _ in range(cnt)]
    return lat, expanded, counts, np.array(coords), (flags if sd and len(flags) == n else None)


def basic_poscar(path, lat, species_expanded, counts, coords_frac, title="fixture",
                 sd_flags=None):
    """Write a POSCAR/CONTCAR; if sd_flags given, emit the Selective dynamics
    header (line 7) + per-atom flags (a real VASP CONTCAR layout)."""
    n = sum(counts)
    species = []
    for sp in species_expanded:  # dedup consecutive runs (compact form)
        if not species or species[-1] != sp:
            species.append(sp)
    lines = [title, "1.0"]
    for row in lat:
        lines.append("  " + " ".join(f"{v: .12f}" for v in row))
    lines.append(" ".join(species))
    lines.append(" ".join(str(c) for c in counts))
    if sd_flags:
        lines.append("Selective dynamics")
    lines.append("Direct")
    for i in range(n):
        ln = " ".join(f"{v: .12f}" for v in coords_frac[i])
        if sd_flags:
            ln += "  " + " ".join(sd_flags[i])
        lines.append(ln)
    path.write_text("\n".join(lines) + "\n")
    return path
