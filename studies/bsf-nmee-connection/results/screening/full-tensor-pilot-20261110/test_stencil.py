"""Geometry-only tests for the directed MMN finite-difference stencil."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_stencil import complete_stencil  # noqa: E402


def _write_case(
    root: Path,
    cell: np.ndarray,
    steps: tuple[tuple[int, int], ...],
    mesh: tuple[int, int] = (3, 3),
) -> tuple[Path, Path]:
    """Write a one-band cache and exact identity-link MMN fixture."""

    n1, n2 = mesh
    points = [(i, j) for i in range(n1) for j in range(n2)]
    lookup = {point: index for index, point in enumerate(points)}
    kpoints = np.array([[i / n1, j / n2, 0.0] for i, j in points])
    cache = root / "cache.npz"
    np.savez_compressed(
        cache,
        energies=np.zeros((len(points), 1)),
        k=kpoints,
        mesh=np.asarray(mesh),
        cell=cell,
    )

    seed = root / "seed"
    seed.mkdir()
    lines = ["synthetic MMN\n", f"1 {len(points)} {len(steps)}\n"]
    for source, (i, j) in enumerate(points):
        for step_i, step_j in steps:
            raw_i, raw_j = i + step_i, j + step_j
            destination = lookup[(raw_i % n1, raw_j % n2)]
            shift_i, shift_j = raw_i // n1, raw_j // n2
            lines.append(
                f"{source + 1} {destination + 1} {shift_i} {shift_j} 0\n"
            )
            lines.append("1.0 0.0\n")
    (seed / "wannier90.mmn").write_text("".join(lines))
    return seed, cache


class CompleteStencilTests(unittest.TestCase):
    def test_six_directed_hex_links_have_exact_cartesian_first_moment(self) -> None:
        root_context = tempfile.TemporaryDirectory()
        self.addCleanup(root_context.cleanup)
        root = Path(root_context.name)
        lattice = 4.1
        cell = np.array(
            [
                [lattice, 0.0, 0.0],
                [0.5 * lattice, np.sqrt(3.0) * lattice / 2.0, 0.0],
                [0.0, 0.0, 22.0],
            ]
        )
        directed_hex = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1))
        seed, cache = _write_case(root, cell, directed_hex)
        output = root / "full-shell.npz"

        report = complete_stencil(seed, cache, output)
        with np.load(output) as data:
            steps = np.asarray(data["stencil_steps"])
            coefficients = np.asarray(data["derivative_weights"])
            links = np.asarray(data["links"])
            neighbors = np.asarray(data["neighbors"])
            mesh = np.asarray(data["mesh"])

        b_reduced = np.column_stack(
            [steps / mesh, np.zeros(len(steps), dtype=float)]
        )
        b_cart = (b_reduced @ (2.0 * np.pi * np.linalg.inv(cell).T))[:, :2]
        first_moment = np.einsum("di,dj->ij", coefficients, b_cart)

        np.testing.assert_allclose(first_moment, np.eye(2), atol=2e-14)
        np.testing.assert_allclose(coefficients.sum(axis=0), 0.0, atol=2e-14)
        self.assertLess(report["first_moment_error"], 2e-14)
        self.assertEqual(len(report["shell_weights_A2"]), 1)
        self.assertTrue(np.all(np.asarray(report["shell_weights_A2"]) > 0.0))
        np.testing.assert_allclose(links, 1.0, atol=0.0)

        step_index = {tuple(step): index for index, step in enumerate(steps)}
        for direction, step in enumerate(steps):
            reverse = step_index[tuple(-step)]
            np.testing.assert_array_equal(
                neighbors[neighbors[:, direction], reverse],
                np.arange(len(neighbors)),
            )
            np.testing.assert_allclose(
                coefficients[reverse], -coefficients[direction], atol=2e-14
            )

        saved_report = json.loads(output.with_suffix(".json").read_text())
        self.assertEqual(saved_report["stencil_steps"], steps.tolist())

    def test_unsupported_oblique_four_link_shell_fails_closed(self) -> None:
        root_context = tempfile.TemporaryDirectory()
        self.addCleanup(root_context.cleanup)
        root = Path(root_context.name)
        cell = np.array(
            [[4.0, 0.0, 0.0], [1.2, 5.0, 0.0], [0.0, 0.0, 22.0]]
        )
        primitive_only = ((1, 0), (-1, 0), (0, 1), (0, -1))
        seed, cache = _write_case(root, cell, primitive_only)

        with self.assertRaisesRegex(
            ValueError, "No positive complete first-derivative shell"
        ):
            complete_stencil(seed, cache, root / "unsupported.npz")


if __name__ == "__main__":
    unittest.main(verbosity=2)
