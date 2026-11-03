#!/usr/bin/env python3
"""AgF2 B1 -> real-material IGMC: read-only feasibility analysis.

Reads ONLY existing local records (POSCAR.soc_input, RUN_RECEIPT.json,
BANDS_GAP.json, MOMENTS.csv) from the archived PR #38 SOC screen. No VASP,
no Wannier, no new material computation. Outputs a JSON evidence receipt.

Checks:
  1. {m_x | (1/2,1/2,0)} glide residual of the B1 screening geometry;
  2. what the prior SOC run produced and what it did NOT (operator chain);
  3. cost model from the prior run (rank-hours, iterations, memory);
  4. Wannier manifold estimate (spinor band counts, frozen-window guess).
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
S1 = (ROOT.parents[0] / "stage2_real_dft" / "b1_wsl_relax" / "idea_screen_1"
      / "soc_static_1")

POS = np.array([
    [1.924240, 2.943540, 15.496370], [3.816300, 5.883370, 15.496370],
    [0.012060, 1.991130, 11.440790], [1.901460, 4.930960, 11.440780],
    [0.431320, 1.830810, 14.530410], [2.347100, 1.126110, 16.488970],
    [3.393440, 4.065940, 16.488970], [1.482190, 4.770640, 14.530410],
    [2.364210, 0.862200, 12.412410], [0.439940, 0.183300, 10.446670],
    [1.473580, 3.123130, 10.446680], [3.376330, 3.802040, 12.412410],
])
CELL = np.array([
    [3.8270249700655650, 0, 0],
    [0, 5.8796657458083281, 0],
    [0, 0, 26.9385414400000000],
])


def glide_residual(pos, cell):
    frac = pos @ np.linalg.inv(cell)
    g = np.array([0.5 - frac[:, 0], 0.5 + frac[:, 1], frac[:, 2]]).T % 1.0
    worst, pairs = 0.0, []
    for i in range(len(pos)):
        d = g[i:i + 1] - frac
        d -= np.round(d)
        dist = np.sqrt(np.sum((d @ cell) ** 2, axis=1))
        j = int(np.argmin(dist))
        pairs.append([i + 1, j + 1, float(dist[j])])
        worst = max(worst, float(dist[j]))
    return worst, pairs


def main():
    rec = {
        "schema": "agf2-real-igmc-feasibility-v1",
        "date": "2026-10-24",
        "sources": {
            "poscar": "POSCAR.soc_input (sha256 "
                      "b4f1aa46c00faf455e4b0af4f323b3e1a502fea5b27692cabe583199ba2befd4)",
            "receipt": "RUN_RECEIPT.json (2026-10-22 SOC screen, PASS_SCREEN)",
        },
        "no_compute_executed": True,
    }
    # 1) glide symmetry of the screening geometry
    worst, pairs = glide_residual(POS, CELL)
    rec["glide_mx_half_half_max_residual_ang"] = worst
    rec["glide_pairs_atom1_to_atom2"] = [
        [p[0], p[1]] for p in pairs]
    # 2) inventory: what exists vs what is missing (operator chain)
    rec["existing"] = {
        "soc_energy_ev": -35.55001581,
        "soc_lowering_ev": -0.122071,
        "gap_ev": 1.3616, "homo_ev": -5.7746, "lumo_ev": -4.4130,
        "k_mesh": "12x8x1 Gamma", "smearing": "ISMEAR=0 SIGMA=0.05",
        "ag_spin_moments_mub": [-0.568, 0.568, 0.568, -0.568],
        "ag_orbital_moments_Lz_mub": 0.070,
        "occupations": "clean 0/1 (insulator at screening grid)",
        "geometry": "step-11 ionic screening geometry, max|F|=0.078 eV/A",
    }
    rec["missing_for_real_igmc"] = [
        "WAVECAR / wavefunctions (LWAVE=.FALSE.; no Bloch states saved)",
        "k-resolved spinor H(k) matrix (Bloch or Wannier)",
        "same-gauge Zeeman spin matrix S_j(k) (site moments != operator)",
        "position / Berry connection A_i(k) (Wannier r.dat, final gauge)",
        "VASP->Wannier90 amn/mmn/spn/eig/u/v files",
        "multiband Dxx/omega value and Wilson cross-check on the material",
    ]
    # 3) cost model from the prior run
    r = json.loads((S1 / "RUN_RECEIPT.json").read_text())
    sec = r["execution"]["elapsed_sec"]
    iters = r["execution"]["n_scf_iterations"]
    mem_kb = r["execution"]["max_memory_rank0_kb"]
    ranks = 8
    rec["prior_run_cost_model"] = {
        "wall_sec": sec, "wall_h": round(sec / 3600, 3),
        "rank_hours": round(sec * ranks / 3600, 2),
        "scf_iterations": iters,
        "sec_per_iteration_8ranks": round(sec / iters, 1),
        "max_memory_rank0_MB": round(mem_kb / 1024, 1),
        "n_bands_hint": "not recorded in receipt; NBANDS is VASP default "
                        "for the noncollinear run (50 occupied spinor bands "
                        "for 100 e-)",
    }
    # 4) Wannier manifold estimate
    nelec = 100                      # Ag 4d10 5s1 (11e) x4 + F 2p5 (7e) x8
    n_occ_spinor = nelec // 2        # 50 occupied spinor bands (noncollinear)
    wf_options = {
        "d+p (Ag-4d 20 + F-2p 24)": 44,
        "d+p+s (add Ag-5s 4)": 48,
        "d+p+s+Ag-4p (deep window)": 60,
    }
    rec["wannier_manifold_estimate"] = {
        "nelec": nelec, "n_occupied_spinor_bands": n_occ_spinor,
        "num_wann_options": wf_options,
        "frozen_window_guess_ev": [-7.0, -3.0],
        "note": "window/disentanglement finalized only in the authorized "
                "run against the actual band structure; acceptance: Wannier "
                "interpolation within 10 meV of the direct eigenvalues in "
                "the frozen window around HOMO/LUMO",
    }
    out = ROOT / "AGF2_REAL_IGMC_FEASIBILITY_EVIDENCE.json"
    out.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
