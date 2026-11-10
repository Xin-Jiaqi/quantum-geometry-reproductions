"""Self-contained numerical tests for :mod:`tensor_response`.

The fixtures are finite-dimensional periodic Bloch models built in an orbital
basis and exported to an eigenstate cache.  They require no DFT files and are
small enough to run with ``python3 -m unittest``.
"""

from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parent))
from tensor_response import (  # noqa: E402
    C_UNITS,
    E_LABELS,
    LINK_DIRECTION_LABELS,
    SCOPE_LABELS,
    SPIN_LABELS,
    TensorResponse,
)


def adjoint(array: np.ndarray) -> np.ndarray:
    return np.asarray(array).conj().swapaxes(-1, -2)


def random_hermitian(
    rng: np.random.Generator, size: int, scale: float = 1.0
) -> np.ndarray:
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    return scale * (raw + raw.conj().T) / (2.0 * np.sqrt(size))


def random_unitary(rng: np.random.Generator, size: int) -> np.ndarray:
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    q, r = np.linalg.qr(raw)
    phases = np.diag(r)
    phases = np.where(np.abs(phases) > 0.0, phases / np.abs(phases), 1.0)
    return q * phases.conj()[None, :]


def make_cache(
    mesh: tuple[int, int] = (3, 3),
    shift: tuple[float, float] = (0.071, -0.113),
    *,
    seed: int = 813,
    cell: np.ndarray | None = None,
    stencil_steps: np.ndarray | None = None,
    derivative_weights: np.ndarray | None = None,
) -> dict[str, np.ndarray | int]:
    """Build a periodic five-band cache with a robust two-band gap."""

    rng = np.random.default_rng(seed)
    n1, n2 = mesh
    nb, nocc = 5, 2
    if cell is None:
        cell = np.array(
            [[4.2, 0.35, 0.0], [0.8, 5.1, 0.0], [0.0, 0.0, 22.0]]
        )
    custom_stencil = stencil_steps is not None
    if stencil_steps is None:
        steps = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]], dtype=int)
    else:
        steps = np.asarray(stencil_steps, dtype=int)
        if derivative_weights is None:
            raise ValueError("custom stencil requires derivative_weights")

    points = [(i, j) for i in range(n1) for j in range(n2)]
    lookup = {point: index for index, point in enumerate(points)}
    k = np.array(
        [[shift[0] + i / n1, shift[1] + j / n2, 0.0] for i, j in points]
    )
    neighbors = np.empty((n1 * n2, len(steps)), dtype=int)
    for index, (i, j) in enumerate(points):
        neighbors[index] = [
            lookup[((i + step[0]) % n1, (j + step[1]) % n2)]
            for step in steps
        ]

    diagonal = np.diag([-2.4, -1.25, 0.9, 1.75, 2.7])
    fourier = [random_hermitian(rng, nb, 0.12) for _ in range(6)]
    orbital_h = []
    for kx, ky, _ in k:
        orbital_h.append(
            diagonal
            + np.cos(2.0 * np.pi * kx) * fourier[0]
            + np.sin(2.0 * np.pi * kx) * fourier[1]
            + np.cos(2.0 * np.pi * ky) * fourier[2]
            + np.sin(2.0 * np.pi * ky) * fourier[3]
            + np.cos(2.0 * np.pi * (kx + ky)) * fourier[4]
            + np.sin(2.0 * np.pi * (kx - ky)) * fourier[5]
        )
    energies, eigenvectors = np.linalg.eigh(np.asarray(orbital_h))

    links = np.empty((len(k), len(steps), nb, nb), dtype=complex)
    for index in range(len(k)):
        for direction in range(len(steps)):
            links[index, direction] = (
                adjoint(eigenvectors[index])
                @ eigenvectors[neighbors[index, direction]]
            )

    z_orbital = random_hermitian(rng, nb, 0.7)
    lower_orbital = np.stack(
        [random_hermitian(rng, nb, 0.55) for _ in range(3)]
    )
    upper_orbital = np.stack(
        [random_hermitian(rng, nb, 0.45) for _ in range(3)]
    )
    spin_orbital = np.stack(
        [lower_orbital + upper_orbital, lower_orbital, upper_orbital]
    )

    z = np.empty((len(k), nb, nb), dtype=complex)
    spin = np.empty((len(k), 3, 3, nb, nb), dtype=complex)
    for index, basis in enumerate(eigenvectors):
        z[index] = adjoint(basis) @ z_orbital @ basis
        for scope in range(3):
            for component in range(3):
                spin[index, scope, component] = (
                    adjoint(basis) @ spin_orbital[scope, component] @ basis
                )

    cache = {
        "energies": energies,
        "links": links,
        "neighbors": neighbors,
        "mesh": np.asarray(mesh, dtype=int),
        "cell": np.asarray(cell, dtype=float),
        "k": k,
        "nocc": nocc,
        "spin_tensor": spin,
        "Z": z,
    }
    if custom_stencil:
        cache["stencil_steps"] = steps
        cache["derivative_weights"] = np.asarray(derivative_weights, dtype=float)
    return cache


def stencil_geometry(
    cell: np.ndarray, mesh: tuple[int, int], steps: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return Cartesian reciprocal displacements and isotropic weights."""

    reciprocal_rows = 2.0 * np.pi * np.linalg.inv(cell).T
    b_cartesian = (
        np.asarray(steps, dtype=float) / np.asarray(mesh, dtype=float)
    ) @ reciprocal_rows[:2]
    b_cartesian = b_cartesian[:, :2]
    moment = b_cartesian.T @ b_cartesian
    weights = b_cartesian @ np.linalg.inv(moment)
    return b_cartesian, weights


def make_c3_scalar_cache(
    steps: np.ndarray, *, explicit_stencil: bool
) -> tuple[dict[str, object], float]:
    """Make a C3-covariant model where a z-spin/in-plane-E vector is forbidden."""

    rng = np.random.default_rng(404)
    size = 7
    mesh = (size, size)
    cell = np.array(
        [[4.0, 0.0, 0.0], [2.0, 2.0 * np.sqrt(3.0), 0.0], [0.0, 0.0, 20.0]]
    )
    points = [(i, j) for i in range(size) for j in range(size)]
    lookup = {point: index for index, point in enumerate(points)}
    k = np.array([[i / size, j / size, 0.0] for i, j in points])

    # In reduced reciprocal coordinates this is a physical +120 degree turn.
    c3_k = np.array([[0, -1], [1, -1]], dtype=int)
    omega = np.exp(2j * np.pi / 3.0)
    c3_orbital = np.diag([1.0, omega, omega**2, 1.0])
    diagonal = np.diag([-2.4, -1.2, 0.9, 2.2])
    harmonics = [random_hermitian(rng, 4, 0.22) for _ in range(4)]

    def unsymmetrized(q: np.ndarray) -> np.ndarray:
        x, y = q
        return (
            diagonal
            + np.cos(2.0 * np.pi * x) * harmonics[0]
            + np.sin(2.0 * np.pi * x) * harmonics[1]
            + np.cos(2.0 * np.pi * y) * harmonics[2]
            + np.sin(2.0 * np.pi * y) * harmonics[3]
        )

    def hamiltonian(q: np.ndarray) -> np.ndarray:
        result = np.zeros((4, 4), dtype=complex)
        rotated_q = np.asarray(q)
        orbital_power = np.eye(4, dtype=complex)
        for _ in range(3):
            result += (
                adjoint(orbital_power)
                @ unsymmetrized(rotated_q)
                @ orbital_power
            )
            rotated_q = c3_k @ rotated_q
            orbital_power = orbital_power @ c3_orbital
        return result / 3.0

    orbital_h = np.array([hamiltonian(point[:2]) for point in k])
    covariance_error = max(
        float(
            np.max(
                np.abs(
                    hamiltonian(c3_k @ point[:2])
                    - c3_orbital
                    @ hamiltonian(point[:2])
                    @ adjoint(c3_orbital)
                )
            )
        )
        for point in k
    )
    energies, eigenvectors = np.linalg.eigh(orbital_h)

    neighbors = np.array(
        [
            [
                lookup[((i + step[0]) % size, (j + step[1]) % size)]
                for step in steps
            ]
            for i, j in points
        ]
    )
    links = np.array(
        [
            [adjoint(eigenvectors[index]) @ eigenvectors[neighbor] for neighbor in row]
            for index, row in enumerate(neighbors)
        ]
    )

    # Both vertices commute with C3, so the resulting (C_zx, C_zy) must be a
    # C3-invariant in-plane vector and hence vanish for a C3-complete stencil.
    z_orbital = np.diag([0.7, -0.3, 0.2, -0.6])
    spin_z_orbital = np.diag([0.4, -0.7, 0.9, 0.1])
    z = np.array([adjoint(basis) @ z_orbital @ basis for basis in eigenvectors])
    spin = np.zeros((len(k), 3, 3, 4, 4), dtype=complex)
    for index, basis in enumerate(eigenvectors):
        spin_z = adjoint(basis) @ spin_z_orbital @ basis
        spin[index, 0, 2] = spin_z
        spin[index, 1, 2] = 0.45 * spin_z
        spin[index, 2, 2] = 0.55 * spin_z

    cache: dict[str, object] = {
        "energies": energies,
        "links": links,
        "neighbors": neighbors,
        "mesh": np.asarray(mesh),
        "cell": cell,
        "k": k,
        "nocc": 2,
        "spin_tensor": spin,
        "Z": z,
    }
    if explicit_stencil:
        _, weights = stencil_geometry(cell, mesh, steps)
        cache["stencil_steps"] = np.asarray(steps)
        cache["derivative_weights"] = weights
    return cache, covariance_error


def legacy_x_response(cache: dict[str, object], b: float) -> np.ndarray:
    """Literal NumPy transcription of the legacy square-mesh x kernel."""

    energies = np.asarray(cache["energies"])
    nk, nb = energies.shape
    nocc = int(cache["nocc"])
    scalar_mesh = round(nk**0.5)
    if scalar_mesh * scalar_mesh != nk:
        raise ValueError("legacy reference only supports a square mesh")
    h = energies[:, :, None] * np.eye(nb)
    z = np.asarray(cache["Z"])
    links = np.asarray(cache["links"])
    neighbors = np.asarray(cache["neighbors"])
    cell = np.asarray(cache["cell"])
    spin = np.asarray(cache["spin_tensor"])
    area = np.linalg.norm(np.cross(cell[0], cell[1]))
    factor = 1.0e10 * (2.00231930436256 * 5.7883818060e-5 / 2.0) / area

    def chi(matrix: np.ndarray) -> np.ndarray:
        eig, vectors = np.linalg.eigh(matrix)
        occupied = vectors[:, :, :nocc]
        empty = vectors[:, :, nocc:]
        derivative = np.zeros_like(occupied, dtype=complex)
        for direction in range(4):
            transported = links[:, direction] @ occupied[neighbors[:, direction]]
            overlap = adjoint(occupied) @ transported
            derivative += (
                (-1) ** direction
                * scalar_mesh
                * cell[direction // 2, 0]
                / (4.0 * np.pi)
                * transported
                @ np.linalg.inv(overlap)
            )
        dipole = adjoint(empty) @ (1j * derivative)
        denominator = eig[:, :nocc, None] - eig[:, None, nocc:]
        values = []
        for operator in (spin[:, 0, 0], spin[:, 1, 0]):
            cross = adjoint(occupied) @ operator @ empty
            values.append(
                float(
                    np.mean(
                        np.sum(
                            2.0 * (cross * dipole.swapaxes(1, 2)).real
                            / denominator,
                            axis=(1, 2),
                        )
                    )
                    * factor
                )
            )
        return np.array([values[0], values[1], values[0] - values[1]])

    return (chi(h + b * z) - chi(h - b * z)) / (2.0 * b)


def direct_tensor_response(cache: dict[str, object], b: float) -> np.ndarray:
    """Independent loop reference for rectangular meshes and all directions."""

    energies = np.asarray(cache["energies"])
    nk, nb = energies.shape
    nocc = int(cache["nocc"])
    h = energies[:, :, None] * np.eye(nb)
    z = np.asarray(cache["Z"])
    links = np.asarray(cache["links"])
    neighbors = np.asarray(cache["neighbors"])
    mesh = np.asarray(cache["mesh"])
    cell = np.asarray(cache["cell"])
    operators = np.asarray(cache["spin_tensor"])
    factor = (
        1.0e10
        * (2.00231930436256 * 5.7883818060e-5 / 2.0)
        / np.linalg.norm(np.cross(cell[0], cell[1]))
    )

    def chi(matrix: np.ndarray) -> np.ndarray:
        eig, vectors = np.linalg.eigh(matrix)
        occupied = vectors[:, :, :nocc]
        empty = vectors[:, :, nocc:]
        derivatives = np.zeros((2, nk, nb, nocc), dtype=complex)
        for direction in range(4):
            axis = direction // 2
            sign = 1 if direction % 2 == 0 else -1
            transported = links[:, direction] @ occupied[neighbors[:, direction]]
            overlap = adjoint(occupied) @ transported
            parallel = transported @ np.linalg.inv(overlap)
            for electric in range(2):
                derivatives[electric] += (
                    sign
                    * mesh[axis]
                    * cell[axis, electric]
                    / (4.0 * np.pi)
                    * parallel
                )
        denominator = eig[:, :nocc, None] - eig[:, None, nocc:]
        result = np.empty((3, 3, 2))
        for scope in range(3):
            for component in range(3):
                cross = adjoint(occupied) @ operators[:, scope, component] @ empty
                for electric in range(2):
                    dipole = adjoint(empty) @ (1j * derivatives[electric])
                    result[scope, component, electric] = (
                        np.mean(
                            np.sum(
                                2.0 * (cross * dipole.swapaxes(1, 2)).real
                                / denominator,
                                axis=(1, 2),
                            )
                        )
                        * factor
                    )
        return result

    return (chi(h + b * z) - chi(h - b * z)) / (2.0 * b)


def gauge_transform(
    cache: dict[str, object], gauges: np.ndarray
) -> tuple[dict[str, object], np.ndarray]:
    """Change the full cached basis independently at each k point."""

    transformed = copy.deepcopy(cache)
    energies = np.asarray(cache["energies"])
    nk, nb = energies.shape
    h = np.zeros((nk, nb, nb), dtype=complex)
    z = np.empty_like(np.asarray(cache["Z"]), dtype=complex)
    spin = np.empty_like(np.asarray(cache["spin_tensor"]), dtype=complex)
    links = np.empty_like(np.asarray(cache["links"]), dtype=complex)
    neighbors = np.asarray(cache["neighbors"])
    for index in range(nk):
        h0 = np.diag(energies[index])
        h[index] = adjoint(gauges[index]) @ h0 @ gauges[index]
        z[index] = (
            adjoint(gauges[index]) @ np.asarray(cache["Z"])[index] @ gauges[index]
        )
        for scope in range(3):
            for component in range(3):
                spin[index, scope, component] = (
                    adjoint(gauges[index])
                    @ np.asarray(cache["spin_tensor"])[index, scope, component]
                    @ gauges[index]
                )
        for direction in range(4):
            neighbor = neighbors[index, direction]
            links[index, direction] = (
                adjoint(gauges[index])
                @ np.asarray(cache["links"])[index, direction]
                @ gauges[neighbor]
            )
    transformed["Z"] = z
    transformed["spin_tensor"] = spin
    transformed["links"] = links
    return transformed, h


class TensorResponseTests(unittest.TestCase):
    def test_legacy_x_regression_is_nonzero(self) -> None:
        cache = make_cache(mesh=(3, 3), seed=10)
        b = 0.0017
        expected = legacy_x_response(cache, b)
        actual = TensorResponse(cache).response(b=b)["C"][:, 0, 0]
        self.assertGreater(float(np.max(np.abs(expected))), 1.0e-3)
        np.testing.assert_allclose(actual, expected, rtol=2.0e-11, atol=2.0e-9)

    def test_shifted_rectangular_mesh_matches_direct_reference(self) -> None:
        cache = make_cache(mesh=(3, 5), shift=(0.037, 0.081), seed=23)
        kernel = TensorResponse(cache)
        result = kernel.response(b=0.0023)
        expected = direct_tensor_response(cache, 0.0023)

        self.assertEqual(result["C"].shape, (3, 3, 2))
        self.assertEqual(tuple(kernel.mesh), (3, 5))
        self.assertGreater(float(np.linalg.norm(result["C"][:, :, 0])), 1.0e-3)
        self.assertGreater(float(np.linalg.norm(result["C"][:, :, 1])), 1.0e-3)
        np.testing.assert_allclose(result["C"], expected, rtol=3.0e-11, atol=3.0e-9)

    def test_full_gauge_covariance_with_degenerate_occupied_space(self) -> None:
        cache = make_cache(mesh=(3, 4), seed=51)
        # Exact occupied degeneracy exercises subspace, rather than band, gauge
        # covariance.  Empty levels remain separated from the occupied pair.
        cache["energies"] = np.tile(
            np.array([-1.8, -1.8, 0.7, 1.4, 2.2]),
            (np.asarray(cache["k"]).shape[0], 1),
        )
        reference = TensorResponse(cache).response(b=0.0011)["C"]

        rng = np.random.default_rng(992)
        gauges = np.stack(
            [random_unitary(rng, 5) for _ in range(np.asarray(cache["k"]).shape[0])]
        )
        transformed, transformed_h = gauge_transform(cache, gauges)
        actual = TensorResponse(transformed).response(h=transformed_h, b=0.0011)["C"]

        self.assertGreater(float(np.linalg.norm(reference)), 1.0e-3)
        np.testing.assert_allclose(actual, reference, rtol=3.0e-9, atol=3.0e-7)

    def test_cartesian_rotation_covariance(self) -> None:
        cache = make_cache(mesh=(4, 3), seed=74)
        reference = TensorResponse(cache).response(b=0.0015)["C"]
        angle = 0.413
        cosine, sine = np.cos(angle), np.sin(angle)
        rotation = np.array(
            [[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]]
        )

        rotated = copy.deepcopy(cache)
        rotated["cell"] = np.asarray(cache["cell"]) @ rotation.T
        rotated["spin_tensor"] = np.einsum(
            "ip,kspab->ksiab",
            rotation,
            np.asarray(cache["spin_tensor"]),
            optimize=True,
        )
        actual = TensorResponse(rotated).response(b=0.0015)["C"]
        expected = np.einsum(
            "ip,spq,jq->sij",
            rotation,
            reference,
            rotation[:2, :2],
            optimize=True,
        )
        np.testing.assert_allclose(actual, expected, rtol=4.0e-10, atol=4.0e-8)

    def test_six_link_stencil_is_polynomial_exact_and_strictly_validated(self) -> None:
        mesh = (5, 5)
        cell = np.array(
            [[4.0, 0.0, 0.0], [2.0, 2.0 * np.sqrt(3.0), 0.0], [0.0, 0.0, 21.0]]
        )
        steps = np.array(
            [[1, 0], [-1, 0], [0, 1], [0, -1], [-1, -1], [1, 1]]
        )
        displacements, weights = stencil_geometry(cell, mesh, steps)
        cache = make_cache(
            mesh=mesh,
            shift=(0.019, -0.043),
            seed=78,
            cell=cell,
            stencil_steps=steps,
            derivative_weights=weights,
        )
        kernel = TensorResponse(cache)
        response = kernel.response(b=0.0014)

        moment = weights.T @ displacements
        np.testing.assert_allclose(moment, np.eye(2), rtol=0.0, atol=3.0e-16)
        gradient = np.array([1.7, -0.63])
        linear_samples = 2.4 + displacements @ gradient
        np.testing.assert_allclose(
            weights.T @ linear_samples, gradient, rtol=0.0, atol=5.0e-16
        )
        self.assertEqual(response["stencil_steps"].shape, (6, 2))
        self.assertEqual(
            response["neighbor_singular_by_direction_minus_zero_plus"].shape,
            (3, 6),
        )
        self.assertLess(response["stencil_first_moment_max_abs"], 3.0e-16)
        self.assertLess(response["stencil_odd_max_abs"], 3.0e-16)

        bad_moment = copy.deepcopy(cache)
        bad_moment["derivative_weights"] = 1.01 * weights
        with self.assertRaisesRegex(ValueError, "first-moment normalization"):
            TensorResponse(bad_moment)

        bad_odd = copy.deepcopy(cache)
        bad_odd["derivative_weights"] = weights.copy()
        bad_odd["derivative_weights"][[0, 1]] += np.array([0.03, -0.02])
        with self.assertRaisesRegex(ValueError, "not odd under step reversal"):
            TensorResponse(bad_odd)

    def test_complete_c3_star_removes_forbidden_z_spin_vector(self) -> None:
        four_steps = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])
        six_steps = np.array(
            [[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [-1, -1]]
        )
        four_cache, four_covariance_error = make_c3_scalar_cache(
            four_steps, explicit_stencil=False
        )
        six_cache, six_covariance_error = make_c3_scalar_cache(
            six_steps, explicit_stencil=True
        )

        four_forbidden = TensorResponse(four_cache).response(b=0.001)["C"][0, 2]
        six_forbidden = TensorResponse(six_cache).response(b=0.001)["C"][0, 2]

        self.assertLess(four_covariance_error, 1.0e-14)
        self.assertLess(six_covariance_error, 1.0e-14)
        self.assertGreater(float(np.linalg.norm(four_forbidden)), 0.1)
        self.assertLess(float(np.linalg.norm(six_forbidden)), 1.0e-7)

    def test_rotated_explicit_stencil_is_cartesian_covariant(self) -> None:
        mesh = (4, 4)
        cell = np.array(
            [[4.3, 0.0, 0.0], [2.15, 2.15 * np.sqrt(3.0), 0.0], [0.0, 0.0, 20.0]]
        )
        steps = np.array(
            [[1, 0], [-1, 0], [0, 1], [0, -1], [-1, -1], [1, 1]]
        )
        _, weights = stencil_geometry(cell, mesh, steps)
        cache = make_cache(
            mesh=mesh,
            shift=(0.0, 0.0),
            seed=83,
            cell=cell,
            stencil_steps=steps,
            derivative_weights=weights,
        )
        reference = TensorResponse(cache).response(b=0.0016)["C"]

        angle = -0.371
        cosine, sine = np.cos(angle), np.sin(angle)
        rotation = np.array(
            [[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]]
        )
        rotated = copy.deepcopy(cache)
        rotated["cell"] = cell @ rotation.T
        rotated["derivative_weights"] = weights @ rotation[:2, :2].T
        rotated["spin_tensor"] = np.einsum(
            "ip,kspab->ksiab",
            rotation,
            np.asarray(cache["spin_tensor"]),
            optimize=True,
        )
        actual = TensorResponse(rotated).response(b=0.0016)["C"]
        expected = np.einsum(
            "ip,spq,jq->sij",
            rotation,
            reference,
            rotation[:2, :2],
            optimize=True,
        )
        np.testing.assert_allclose(actual, expected, rtol=6.0e-10, atol=6.0e-8)

    def test_diagnostics_layer_sum_and_b_parity(self) -> None:
        cache = make_cache(mesh=(3, 4), seed=91)
        result = TensorResponse(cache).response(b=0.0018)
        reverse_cache = copy.deepcopy(cache)
        reverse_cache["Z"] = -np.asarray(cache["Z"])
        reversed_result = TensorResponse(reverse_cache).response(b=0.0018)

        np.testing.assert_allclose(
            result["even_chi"],
            (result["chi_plus"] + result["chi_minus"]) / 2.0,
            atol=0.0,
            rtol=0.0,
        )
        np.testing.assert_allclose(
            result["b_even_residual"],
            result["even_chi"] - result["zero_chi"],
            atol=0.0,
            rtol=0.0,
        )
        np.testing.assert_allclose(reversed_result["C"], -result["C"], rtol=1e-11)
        np.testing.assert_allclose(
            reversed_result["even_chi"], result["even_chi"], rtol=1e-12
        )
        self.assertLess(result["layer_sum_max_abs"], 1.0e-8)
        self.assertLess(result["operator_layer_sum_max_abs"], 1.0e-12)
        self.assertGreater(result["gap_eV"], 1.0)
        self.assertGreater(result["min_neighbor_singular"], 0.9)
        self.assertEqual(result["units"], C_UNITS)
        self.assertEqual(result["scope_labels"], SCOPE_LABELS)
        self.assertEqual(result["spin_labels"], SPIN_LABELS)
        self.assertEqual(result["electric_field_labels"], E_LABELS)
        self.assertEqual(result["link_direction_labels"], LINK_DIRECTION_LABELS)

    def test_rejects_bad_shapes_neighbors_hermiticity_and_layer_sum(self) -> None:
        cache = make_cache(mesh=(3, 3), seed=112)

        bad_shape = copy.deepcopy(cache)
        bad_shape["spin_tensor"] = np.asarray(cache["spin_tensor"])[:, :, :2]
        with self.assertRaisesRegex(ValueError, "spin_tensor must have shape"):
            TensorResponse(bad_shape)

        bad_neighbors = copy.deepcopy(cache)
        bad_neighbors["neighbors"][0, 0] = 0
        with self.assertRaisesRegex(ValueError, "neighbors direction"):
            TensorResponse(bad_neighbors)

        bad_z = copy.deepcopy(cache)
        bad_z["Z"][0, 0, 1] += 0.1j
        with self.assertRaisesRegex(ValueError, "Z is not Hermitian"):
            TensorResponse(bad_z)

        bad_layer = copy.deepcopy(cache)
        bad_layer["spin_tensor"][0, 2, 1, 0, 0] += 0.01
        with self.assertRaisesRegex(ValueError, r"total = lower \+ upper"):
            TensorResponse(bad_layer)

    def test_reports_singular_transport_and_closed_indirect_gap(self) -> None:
        cache = make_cache(mesh=(3, 3), seed=131)
        singular = copy.deepcopy(cache)
        swap = np.eye(5, dtype=complex)
        swap[[0, 2]] = swap[[2, 0]]
        index, direction = 0, 0
        neighbor = singular["neighbors"][index, direction]
        singular["links"][index, direction] = swap
        singular["links"][neighbor, direction ^ 1] = adjoint(swap)
        with self.assertRaisesRegex(ValueError, "occupied neighboring spaces are singular"):
            TensorResponse(singular).chi()

        kernel = TensorResponse(cache)
        h = np.zeros_like(kernel.h)
        diagonal = np.arange(5)
        values = np.tile(np.array([-2.0, -1.0, 1.0, 2.0, 3.0]), (kernel.nk, 1))
        values[0] = [-2.0, 0.9, 1.0, 2.0, 3.0]
        values[1] = [-2.0, -1.0, -0.8, 2.0, 3.0]
        h[:, diagonal, diagonal] = values
        with self.assertRaisesRegex(ValueError, "global insulating gap closed"):
            kernel.chi(h=h)

    def test_rejects_nonhermitian_h_override(self) -> None:
        kernel = TensorResponse(make_cache(mesh=(3, 3), seed=151))
        h = kernel.h.copy()
        h[0, 0, 1] = 0.1
        with self.assertRaisesRegex(ValueError, "h is not Hermitian"):
            kernel.response(h=h)


if __name__ == "__main__":
    unittest.main(verbosity=2)
