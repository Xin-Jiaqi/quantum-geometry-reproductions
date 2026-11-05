"""Static longitudinal layer and total spin response at fixed filling and stacking."""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(STUDY), str(STUDY.parents[1] / "papers/jia-2026-nonlin-magnetoelectric-edelstein")]
from bsf_nmee.bilayer_audit import AuditBilayer, AuditSpec, integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift
from bsf_nmee.sliding_tight_binding import P_TOP, P_BOTTOM
from nmee_lab.clean_core.core import MU_B


def fixed_stack_matrices(model, k):
    # The existing hopping is linear in t. Cache the k-dependent matrices once.
    h = model.matrices(k)
    v = model.matrices(k, axis=0)
    return h, v


def response(model, h_unit, v_unit, weights, t_mev, bias_mev, step):
    """Global Zeeman field, layer-resolved readout, zero frequency, half filling."""
    response_parts = np.zeros(3)
    positive = negative = 0.
    # Track a common insulating chemical potential at 0 and both derivative fields.
    vmax, cmin = np.full(3, -np.inf), np.full(3, np.inf)
    for start in range(0, len(weights), 512):
        chunk = slice(start, start + 512)
        h, v = h_unit[chunk].copy(), v_unit[chunk].copy()
        for matrix in (h, v):
            matrix[:, :4, 4:] *= t_mev / 1000
            matrix[:, 4:, :4] *= t_mev / 1000
        h += (bias_mev / 2000) * (P_TOP - P_BOTTOM)
        bottom_k = np.zeros(len(h))
        for index, direction in enumerate((0, 1, -1)):
            if direction:
                field_h = h - direction * step * model.mono.g_factor * MU_B * model.spins["total_z"]
                e, u = np.linalg.eigh(field_h)
            else:
                e = np.linalg.eigvalsh(h)
            vmax[index] = max(vmax[index], e[:, 3].max())
            cmin[index] = min(cmin[index], e[:, 4].min())
            if direction:
                occ, empty = u[:, :, :4], u[:, :, 4:]
                velocity = (empty.conj().transpose(0, 2, 1) @ v @ occ).transpose(0, 2, 1)
                kernel = velocity / (e[:, :4, None] - e[:, None, 4:]) ** 2
                for i, obs in enumerate(("bottom_x", "top_x", "total_x")):
                    spin = occ.conj().transpose(0, 2, 1) @ model.spins[obs] @ empty
                    chi = -2 * (spin * kernel).imag.sum(axis=(1, 2))
                    response_parts[i] += direction * (weights[chunk] @ chi) / (2 * step)
                    if i == 0:
                        bottom_k += direction * chi / (2 * step)
        positive += weights[chunk] @ np.maximum(bottom_k, 0)
        negative += weights[chunk] @ np.minimum(bottom_k, 0)
    gap = cmin[0] - vmax[0]
    common_gap = min(cmin) - max(vmax)
    valid = common_gap > 0
    bottom, top, total = response_parts if valid else (None, None, None)
    return dict(t_meV=t_mev, U_top_minus_bottom_meV=bias_mev, gap_meV=1000*gap,
                common_gap_meV=1000*common_gap, bottom=bottom, top=top, total=total,
                layer_sum_error=bottom+top-total if valid else None,
                bottom_positive_k=positive if valid else None,
                bottom_negative_k=negative if valid else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mass", type=float, default=10., help="sublattice onsite term in meV")
    parser.add_argument("--fraction", type=float, default=.25)
    parser.add_argument("--grid", nargs=3, type=int, default=(41, 45, 72))
    parser.add_argument("--step", type=float, default=.0625, help="Zeeman step, tesla")
    parser.add_argument("--points", required=True, help="t:U pairs in meV, comma separated")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not np.isfinite(args.step) or args.step <= 0:
        parser.error("--step must be positive and finite")
    k, w, _ = integration_grid(AuditBilayer(), *args.grid)
    # The only costly k dependence of tunneling is computed in bounded batches.
    model = AuditBilayer(AuditSpec(shift=ab_ba_shift(args.fraction),
                                 staggered_eV=args.mass / 1000, tunneling_eV=1.0))
    pieces = [fixed_stack_matrices(model, k[i:i+512])
              for i in range(0, len(k), 512)]
    h = np.concatenate([part[0] for part in pieces])
    v = np.concatenate([part[1] for part in pieces])
    del pieces
    with args.output.open("w") as stream:
        writer = None
        for point in args.points.split(","):
            t, bias = map(float, point.split(":"))
            row = dict(fraction=args.fraction, mass_meV=args.mass, grid="x".join(map(str, args.grid)), step_T=args.step)
            row.update(response(model, h, v, w, t, bias, args.step))
            if writer is None:
                writer = csv.DictWriter(stream, fieldnames=row.keys(), lineterminator="\n")
                writer.writeheader()
            writer.writerow(row)
            stream.flush()
            print(t, bias, row["gap_meV"], row["bottom"], row["total"], flush=True)


if __name__ == "__main__":
    main()
