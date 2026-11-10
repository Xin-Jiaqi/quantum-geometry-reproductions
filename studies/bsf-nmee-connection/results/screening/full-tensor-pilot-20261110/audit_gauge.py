"""Audit WAVECAR/MMN/cache gauge agreement across a derivative stencil.

The audit is deliberately standalone: it does not import ``export_tensor``.
For one direction from every exact +/- pair it checks one interior link and
one link crossing a reciprocal-cell boundary.  The latter explicitly
exercises the plane-wave index shift and atomic PAW phase used at a periodic
boundary.  A legacy cache therefore has four cases; a six-link hexagonal
cache has six, including the added diagonal pair.

Run only inside a PBS compute allocation, for example::

    python audit_gauge.py /path/to/seed /path/to/tensor_cache.npz \
        --mesh 6 6 --output gauge_audit.json
"""

from __future__ import annotations

import argparse
import json
import os
import socket
import sys
from pathlib import Path

import numpy as np


def adjoint(array: np.ndarray) -> np.ndarray:
    return np.asarray(array).conj().swapaxes(-1, -2)


def _serial(value: object) -> object:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def _cache_path(path: Path) -> Path:
    return path / "tensor_cache.npz" if path.is_dir() else path


def _choose_cases(
    kpoints: np.ndarray,
    neighbors: np.ndarray,
    mesh: np.ndarray,
    steps: np.ndarray,
) -> list[dict[str, object]]:
    """Choose interior/boundary links for one member of every +/- pair."""

    cases: list[dict[str, object]] = []
    lookup = {tuple(map(int, step)): direction for direction, step in enumerate(steps)}
    if len(lookup) != len(steps) or (0, 0) in lookup:
        raise ValueError("cache stencil steps must be unique and nonzero")
    representatives = []
    seen: set[tuple[int, int]] = set()
    for direction, step in enumerate(steps):
        key = tuple(map(int, step))
        reverse = tuple(-value for value in key)
        if reverse not in lookup:
            raise ValueError(f"cache stencil step {key} has no negative partner")
        if key in seen:
            continue
        representatives.append(direction)
        seen.update((key, reverse))

    for direction in representatives:
        step = np.asarray(steps[direction], dtype=int)
        target = np.zeros(3)
        target[:2] = step / mesh
        candidates: dict[str, dict[str, object]] = {}
        for source in range(len(kpoints)):
            destination = int(neighbors[source, direction])
            shift_float = kpoints[source] + target - kpoints[destination]
            shift = np.rint(shift_float).astype(int)
            if np.max(np.abs(shift_float - shift)) > 1.0e-7:
                raise ValueError(
                    f"cache step {tuple(step)} neighbor has inconsistent displacement"
                )
            kind = "interior" if not np.any(shift) else "boundary"
            candidate = {
                "step": step,
                "direction": direction,
                "kind": kind,
                "source": source,
                "destination": destination,
                "shift": shift,
            }
            # For a diagonal direction, prefer the corner link that wraps both
            # reduced axes and therefore tests the combined atomic PAW phase.
            if kind not in candidates or np.count_nonzero(shift) > np.count_nonzero(
                candidates[kind]["shift"]
            ):
                candidates[kind] = candidate
        for kind in ("interior", "boundary"):
            if kind not in candidates:
                raise ValueError(
                    f"mesh has no step {tuple(step)} {kind} link for gauge audit"
                )
            cases.append(candidates[kind])
    return cases


def _read_target_mmn(
    path: Path,
    nbands: int,
    nkpoints: int,
    cases: list[dict[str, object]],
) -> dict[int, np.ndarray]:
    """Read matrix data only for the selected MMN records."""

    wanted: dict[tuple[int, int, int, int, int], int] = {}
    for index, case in enumerate(cases):
        shift = np.asarray(case["shift"], dtype=int)
        key = (
            int(case["source"]) + 1,
            int(case["destination"]) + 1,
            int(shift[0]),
            int(shift[1]),
            int(shift[2]),
        )
        if key in wanted:
            raise ValueError(f"duplicate gauge-audit target {key}")
        wanted[key] = index

    found: dict[int, np.ndarray] = {}
    with path.open() as handle:
        next(handle)
        bands, points, neighbor_count = map(int, next(handle).split())
        if (bands, points) != (nbands, nkpoints):
            raise ValueError(
                "MMN dimensions disagree with WAVECAR/cache: "
                f"{(bands, points)} != {(nbands, nkpoints)}"
            )
        for _ in range(points * neighbor_count):
            header = tuple(map(int, next(handle).split()))
            if len(header) != 5:
                raise ValueError(f"invalid MMN link header: {header}")
            target_index = wanted.get(header)
            if target_index is None:
                for _ in range(bands * bands):
                    next(handle)
                continue
            if target_index in found:
                raise ValueError(f"duplicate target MMN record: {header}")
            raw = np.fromstring(
                " ".join(next(handle) for _ in range(bands * bands)), sep=" "
            )
            if raw.size != 2 * bands * bands:
                raise ValueError(f"incomplete MMN matrix for {header}")
            pairs = raw.reshape(-1, 2)
            found[target_index] = (
                pairs[:, 0] + 1j * pairs[:, 1]
            ).reshape(bands, bands).T

    missing = sorted(set(range(len(cases))).difference(found))
    if missing:
        labels = [
            f"{tuple(cases[index]['step'])}/{cases[index]['kind']}"
            for index in missing
        ]
        raise ValueError(f"target MMN records not found: {', '.join(labels)}")
    return found


def _wave_state(wave: object, index: int, nbands: int) -> tuple[np.ndarray, np.ndarray]:
    gvectors = wave.gvectors(index + 1)
    coefficients = np.stack(
        [
            wave.readBandCoeff(ikpt=index + 1, iband=band + 1).reshape(2, -1)
            for band in range(nbands)
        ],
        axis=-1,
    )
    return gvectors, coefficients


def _reconstruct_overlap(
    source: int,
    destination: int,
    shift: np.ndarray,
    *,
    atoms: object,
    wave: object,
    pots: list[object],
    nonlq: object,
    state_cache: dict[int, tuple[np.ndarray, np.ndarray, object]],
) -> np.ndarray:
    """Reconstruct one full-KS PAW overlap directly from WAVECAR data."""

    def state(index: int) -> tuple[np.ndarray, np.ndarray, object]:
        if index not in state_cache:
            gvectors, coefficients = _wave_state(wave, index, wave._nbands)
            projectors = nonlq(
                atoms, wave._encut, potcar=pots, k=wave._kvecs[index]
            )
            if not np.array_equal(gvectors, projectors.Gvec):
                raise ValueError(f"PAW plane-wave ordering mismatch at k={index}")
            state_cache[index] = (gvectors, coefficients, projectors)
        return state_cache[index]

    source_g, source_w, source_p = state(source)
    target_g, target_w, target_p = state(destination)
    target_lookup = {tuple(vector): index for index, vector in enumerate(target_g)}
    source_indices = [
        index
        for index, vector in enumerate(source_g)
        if tuple(vector + shift) in target_lookup
    ]
    target_indices = [
        target_lookup[tuple(source_g[index] + shift)] for index in source_indices
    ]
    if not source_indices:
        raise ValueError("no common plane waves for reconstructed overlap")

    overlap = sum(
        adjoint(source_w[spin, source_indices])
        @ target_w[spin, target_indices]
        for spin in range(2)
    )
    positions = atoms.get_scaled_positions()
    for atom, element_type in enumerate(source_p.element_idx):
        projections = [
            projectors.qproj[element_type]
            * projectors.cqfak[element_type][:, None]
            * projectors.crexp[None, :, atom]
            for projectors in (source_p, target_p)
        ]
        phase = np.exp(-2j * np.pi * shift @ positions[atom])
        qij = pots[element_type].get_Qij()
        for spin in range(2):
            overlap += (
                phase
                * adjoint(projections[0] @ source_w[spin])
                @ qij
                @ (projections[1] @ target_w[spin])
            )
    return overlap


def audit(
    seed: Path,
    cache_path: Path,
    mesh: np.ndarray,
    output: Path,
    tolerance: float,
) -> dict[str, object]:
    from ase.io import read
    from paw import nonlq, pawpotcar
    from vaspwfc import vaspwfc

    if mesh.shape != (2,) or np.any(mesh < 2):
        raise ValueError("--mesh requires two integers of at least 2")
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("--tolerance must be positive and finite")

    cache_path = _cache_path(cache_path)
    cache = np.load(cache_path)
    wave = vaspwfc(str(seed / "WAVECAR"), lsorbit=True)
    atoms = read(seed / "POSCAR")
    pots = [
        pawpotcar(dataset)
        for dataset in (seed / "POTCAR").read_text().split("End of Dataset")[:-1]
    ]
    kpoints = np.asarray(wave._kvecs, dtype=float)
    nkpoints, nbands = len(kpoints), int(wave._nbands)
    if nkpoints != int(np.prod(mesh)):
        raise ValueError("WAVECAR k count does not match --mesh")
    if not np.array_equal(np.asarray(cache["mesh"], dtype=int), mesh):
        raise ValueError("cache mesh does not match --mesh")
    if not np.allclose(cache["k"], kpoints, atol=1.0e-10, rtol=0.0):
        raise ValueError("cache and WAVECAR k-point order differ")
    links = np.asarray(cache["links"])
    neighbors = np.asarray(cache["neighbors"], dtype=int)
    if links.ndim != 4 or links.shape[0] != nkpoints or links.shape[2:] != (
        nbands,
        nbands,
    ):
        raise ValueError("cache links do not contain the full WAVECAR band space")
    directions = links.shape[1]
    if neighbors.shape != (nkpoints, directions):
        raise ValueError("cache neighbors has the wrong shape")
    if "stencil_steps" in cache:
        steps = np.asarray(cache["stencil_steps"])
        if steps.shape != (directions, 2):
            raise ValueError("cache stencil_steps has the wrong shape")
        rounded_steps = np.rint(steps)
        if np.max(np.abs(steps - rounded_steps)) > 1.0e-7:
            raise ValueError("cache stencil_steps must contain integers")
        steps = rounded_steps.astype(int)
    else:
        if directions != 4:
            raise ValueError("non-four-link cache must provide stencil_steps")
        steps = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]], dtype=int)

    nocc = int(np.asarray(cache["nocc"]).reshape(()))
    window = slice(max(0, nocc - 8), min(nbands, nocc + 8))
    cases = _choose_cases(kpoints, neighbors, mesh, steps)
    raw_links = _read_target_mmn(
        seed / "wannier90.mmn", nbands, nkpoints, cases
    )
    state_cache: dict[int, tuple[np.ndarray, np.ndarray, object]] = {}
    rows = []
    for index, case in enumerate(cases):
        source = int(case["source"])
        destination = int(case["destination"])
        direction = int(case["direction"])
        shift = np.asarray(case["shift"], dtype=int)
        raw = raw_links[index]
        reconstructed = _reconstruct_overlap(
            source,
            destination,
            shift,
            atoms=atoms,
            wave=wave,
            pots=pots,
            nonlq=nonlq,
            state_cache=state_cache,
        )
        raw_window = raw[window, window]
        denominator = float(np.linalg.norm(raw_window))
        if denominator <= np.finfo(float).tiny:
            raise ValueError("selected MMN window has zero norm")
        wave_error = float(
            np.linalg.norm((reconstructed - raw)[window, window]) / denominator
        )
        cache_error = float(
            np.linalg.norm((links[source, direction] - raw)[window, window])
            / denominator
        )
        rows.append(
            {
                "label": (
                    f"step_{int(case['step'][0]):+d}_{int(case['step'][1]):+d}_"
                    f"{case['kind']}"
                ),
                "step_mesh_indices": case["step"],
                "source_zero_based": source,
                "destination_zero_based": destination,
                "reciprocal_shift": shift,
                "wave_mmn_relative_frobenius": wave_error,
                "cache_mmn_relative_frobenius": cache_error,
                "cache_mmn_max_abs": float(
                    np.max(np.abs(links[source, direction] - raw))
                ),
                "passed": bool(wave_error <= tolerance and cache_error <= 1.0e-12),
            }
        )

    result = {
        "passed": all(row["passed"] for row in rows),
        "scope": "WAVECAR+PAW versus raw MMN and exported cache; interior and reciprocal-boundary links for every +/- stencil pair",
        "source_seed": str(seed),
        "cache": str(cache_path),
        "mesh": mesh,
        "bands": nbands,
        "nocc": nocc,
        "band_window_zero_based_stop_exclusive": [window.start, window.stop],
        "wave_mmn_relative_tolerance": tolerance,
        "cache_mmn_relative_tolerance": 1.0e-12,
        "wave_mmn_max_relative_frobenius": max(
            row["wave_mmn_relative_frobenius"] for row in rows
        ),
        "cache_mmn_max_relative_frobenius": max(
            row["cache_mmn_relative_frobenius"] for row in rows
        ),
        "cases": rows,
        "node": socket.gethostname(),
        "job": os.environ["PBS_JOBID"],
    }
    output.write_text(json.dumps(result, default=_serial, indent=2) + "\n")
    return result


def main() -> None:
    if (
        not os.getenv("PBS_JOBID")
        or not os.getenv("PBS_NODEFILE")
        or socket.gethostname().startswith(("super", "cms-ssh"))
    ):
        raise SystemExit("PBS compute allocation required")

    parser = argparse.ArgumentParser()
    parser.add_argument("seed", type=Path)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--mesh", nargs=2, type=int, required=True)
    parser.add_argument("--output", type=Path, default=Path("gauge_audit.json"))
    parser.add_argument("--tolerance", type=float, default=0.08)
    args = parser.parse_args()

    parent = Path("/work/sunq/xjq/mosse-7269-mechanism-20261104")
    sys.path[:0] = [str(parent), str(parent / "pawlib")]
    result = audit(
        args.seed,
        args.cache,
        np.asarray(args.mesh, dtype=int),
        args.output,
        args.tolerance,
    )
    print("GAUGE_AUDIT_PASSED", result["passed"], flush=True)
    if not result["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
