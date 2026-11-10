"""Gauge-covariant frozen-potential spin response on a 2D k mesh.

The cache contract is deliberately small and NumPy-only::

    energies[nk, nb]
    links[nk, nd, nb, nb]
    neighbors[nk, nd]
    mesh[2]                    # (N1, N2), with N1 * N2 == nk
    cell[3, 3]                 # direct-lattice vectors as rows, in Angstrom
    k[nk, 3]                  # reduced reciprocal coordinates
    nocc
    spin_tensor[nk, 3, 3, nb, nb]  # total/lower/upper; sx/sy/sz
    Z[nk, nb, nb]              # dimensionless Zeeman probe vertex

An old cache with ``nd == 4`` needs no further fields and uses the exact legacy
``(+a1, -a1, +a2, -a2)`` stencil.  A general stencil additionally supplies::

    stencil_steps[nd, 2]       # integer increments in mesh-index coordinates
    derivative_weights[nd, 2]  # Cartesian x/y derivative weights, in Angstrom

Every step must have its negative partner.  The weights must be odd under that
pairing and obey ``sum_d weights[d,j] * b_cart[d,k] == delta[j,k]``.

``TensorResponse.response`` returns ``C[scope, spin, E]`` with shape
``(3, 3, 2)``.  It preserves the legacy convention

    C = (chi(H + b Z) - chi(H - b Z)) / (2 b)

and the legacy units/factor.  A generic Hermitian ``h[nk, nb, nb]`` may be
passed to ``chi`` or ``response``; the default Hamiltonian is diagonal in the
cached eigenbasis with eigenvalues ``energies``.
"""

from __future__ import annotations

from typing import Mapping, Optional

import numpy as np


SCOPE_LABELS = ("total", "lower", "upper")
SPIN_LABELS = ("x", "y", "z")
E_LABELS = ("x", "y")
LINK_DIRECTION_LABELS = ("+a1", "-a1", "+a2", "-a2")
C_UNITS = "1e-9 (nm V T)^-1, Pauli spin density"

_ELECTRON_G = 2.00231930436256
_MU_B_EV_PER_T = 5.7883818060e-5


def adjoint(array: np.ndarray) -> np.ndarray:
    """Conjugate transpose on the final two axes."""

    return np.asarray(array).conj().swapaxes(-1, -2)


def _max_abs(array: np.ndarray) -> float:
    return float(np.max(np.abs(array), initial=0.0))


def _finite(name: str, array: np.ndarray) -> None:
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values")


def _real_array(name: str, value: object, tolerance: float) -> np.ndarray:
    array = np.asarray(value)
    _finite(name, array)
    if np.iscomplexobj(array) and _max_abs(array.imag) > tolerance:
        raise ValueError(f"{name} must be real")
    return np.asarray(array.real, dtype=float)


class TensorResponse:
    """Full spin/E-direction frozen-potential response kernel.

    Scope and Cartesian-axis order is fixed by the module-level labels.  The
    upper-layer vertex is consumed directly from ``spin_tensor``; it is never
    reconstructed as total minus lower.  The latter relation is instead an
    independently checked cache invariant and an output diagnostic.
    """

    def __init__(
        self,
        cache: Mapping[str, object],
        nb: Optional[int] = None,
        *,
        hermitian_tolerance: float = 1.0e-8,
        link_tolerance: float = 1.0e-6,
        mesh_tolerance: float = 1.0e-7,
        stencil_tolerance: float = 1.0e-7,
        layer_sum_tolerance: float = 1.0e-7,
        singular_tolerance: float = 1.0e-5,
        gap_tolerance: float = 0.0,
    ) -> None:
        required = {
            "energies",
            "links",
            "neighbors",
            "mesh",
            "cell",
            "k",
            "nocc",
            "spin_tensor",
            "Z",
        }
        missing = sorted(required.difference(cache.keys()))
        if missing:
            raise ValueError(f"cache missing required fields: {', '.join(missing)}")

        tolerances = (
            hermitian_tolerance,
            link_tolerance,
            mesh_tolerance,
            stencil_tolerance,
            layer_sum_tolerance,
            singular_tolerance,
        )
        if any(not np.isfinite(x) or x < 0.0 for x in tolerances):
            raise ValueError("validation tolerances must be finite and non-negative")
        if not np.isfinite(gap_tolerance):
            raise ValueError("gap_tolerance must be finite")

        energies_all = _real_array("energies", cache["energies"], hermitian_tolerance)
        if energies_all.ndim != 2:
            raise ValueError("energies must have shape (nk, nb)")
        self.nk, available_nb = energies_all.shape
        if self.nk == 0 or available_nb < 2:
            raise ValueError("energies must contain at least one k point and two bands")

        if nb is None:
            self.nb = available_nb
        else:
            if isinstance(nb, bool) or int(nb) != nb:
                raise ValueError("nb must be an integer")
            self.nb = int(nb)
            if not 1 <= self.nb <= available_nb:
                raise ValueError(f"nb must lie in [1, {available_nb}]")
        self.energies = energies_all[:, : self.nb].copy()

        nocc_array = np.asarray(cache["nocc"])
        if nocc_array.size != 1 or not np.isfinite(nocc_array).all():
            raise ValueError("nocc must be one finite integer")
        nocc_float = float(nocc_array.reshape(()))
        if not nocc_float.is_integer():
            raise ValueError("nocc must be an integer")
        self.nocc = int(nocc_float)
        if not 0 < self.nocc < self.nb:
            raise ValueError("nocc must leave at least one occupied and one empty band")

        mesh_raw = _real_array("mesh", cache["mesh"], mesh_tolerance)
        if mesh_raw.shape != (2,) or np.any(mesh_raw <= 0.0):
            raise ValueError("mesh must have shape (2,) with positive entries")
        mesh_round = np.rint(mesh_raw)
        if _max_abs(mesh_raw - mesh_round) > mesh_tolerance:
            raise ValueError("mesh entries must be integers")
        self.mesh = mesh_round.astype(int)
        if int(np.prod(self.mesh)) != self.nk:
            raise ValueError("mesh product must equal nk")

        self.cell = _real_array("cell", cache["cell"], hermitian_tolerance)
        if self.cell.shape != (3, 3):
            raise ValueError("cell must have shape (3, 3)")
        self.k = _real_array("k", cache["k"], mesh_tolerance)
        if self.k.shape != (self.nk, 3):
            raise ValueError("k must have shape (nk, 3)")

        neighbors_raw = np.asarray(cache["neighbors"])
        _finite("neighbors", neighbors_raw)
        links_all = np.asarray(cache["links"], dtype=complex)
        _finite("links", links_all)
        if (
            links_all.ndim != 4
            or links_all.shape[0] != self.nk
            or links_all.shape[1] == 0
            or links_all.shape[2] < self.nb
            or links_all.shape[3] < self.nb
        ):
            raise ValueError("links must have shape (nk, nd, nb, nb)")
        self.nd = links_all.shape[1]
        self.links = links_all[:, :, : self.nb, : self.nb].copy()

        if neighbors_raw.shape != (self.nk, self.nd):
            raise ValueError("neighbors must have shape (nk, nd)")
        neighbors_round = np.rint(neighbors_raw)
        if _max_abs(neighbors_raw - neighbors_round) > mesh_tolerance:
            raise ValueError("neighbors must contain integer indices")
        self.neighbors = neighbors_round.astype(int)
        if np.any(self.neighbors < 0) or np.any(self.neighbors >= self.nk):
            raise ValueError("neighbors contains an out-of-range index")

        z_all = np.asarray(cache["Z"], dtype=complex)
        _finite("Z", z_all)
        if (
            z_all.ndim != 3
            or z_all.shape[0] != self.nk
            or z_all.shape[1] < self.nb
            or z_all.shape[2] < self.nb
        ):
            raise ValueError("Z must have shape (nk, nb, nb)")
        self.z = z_all[:, : self.nb, : self.nb].copy()

        spin_all = np.asarray(cache["spin_tensor"], dtype=complex)
        _finite("spin_tensor", spin_all)
        if (
            spin_all.ndim != 5
            or spin_all.shape[:3] != (self.nk, 3, 3)
            or spin_all.shape[3] < self.nb
            or spin_all.shape[4] < self.nb
        ):
            raise ValueError(
                "spin_tensor must have shape (nk, 3 scopes, 3 spins, nb, nb)"
            )
        self.spin_tensor = spin_all[:, :, :, : self.nb, : self.nb].copy()

        self.hermitian_tolerance = float(hermitian_tolerance)
        self.link_tolerance = float(link_tolerance)
        self.mesh_tolerance = float(mesh_tolerance)
        self.stencil_tolerance = float(stencil_tolerance)
        self.layer_sum_tolerance = float(layer_sum_tolerance)
        self.singular_tolerance = float(singular_tolerance)
        self.gap_tolerance = float(gap_tolerance)

        self.area_A2 = float(np.linalg.norm(np.cross(self.cell[0], self.cell[1])))
        if not np.isfinite(self.area_A2) or self.area_A2 <= 0.0:
            raise ValueError("cell must have nonzero in-plane area")
        self.factor = (
            1.0e10 * (_ELECTRON_G * _MU_B_EV_PER_T / 2.0) / self.area_A2
        )

        has_steps = "stencil_steps" in cache
        has_weights = "derivative_weights" in cache
        if has_steps != has_weights:
            raise ValueError(
                "stencil_steps and derivative_weights must be supplied together"
            )
        if has_steps:
            steps_raw = _real_array(
                "stencil_steps", cache["stencil_steps"], mesh_tolerance
            )
            if steps_raw.shape != (self.nd, 2):
                raise ValueError("stencil_steps must have shape (nd, 2)")
            steps_round = np.rint(steps_raw)
            if _max_abs(steps_raw - steps_round) > self.mesh_tolerance:
                raise ValueError("stencil_steps must contain integers")
            self.stencil_steps = steps_round.astype(int)
            self.derivative_weights = _real_array(
                "derivative_weights", cache["derivative_weights"], hermitian_tolerance
            )
            if self.derivative_weights.shape != (self.nd, 2):
                raise ValueError("derivative_weights must have shape (nd, 2)")
            self.uses_legacy_four_stencil = False
        else:
            if self.nd != 4:
                raise ValueError(
                    "non-four-link caches require stencil_steps and derivative_weights"
                )
            self.stencil_steps = np.array(
                [[1, 0], [-1, 0], [0, 1], [0, -1]], dtype=int
            )
            self.derivative_weights = np.array(
                [
                    self.mesh[0] * self.cell[0, :2] / (4.0 * np.pi),
                    -self.mesh[0] * self.cell[0, :2] / (4.0 * np.pi),
                    self.mesh[1] * self.cell[1, :2] / (4.0 * np.pi),
                    -self.mesh[1] * self.cell[1, :2] / (4.0 * np.pi),
                ]
            )
            self.uses_legacy_four_stencil = True

        step_lookup = {}
        for direction, step in enumerate(self.stencil_steps):
            key = tuple(int(value) for value in step)
            if key == (0, 0):
                raise ValueError("stencil_steps cannot contain the zero step")
            if key in step_lookup:
                raise ValueError(f"duplicate stencil step {key}")
            step_lookup[key] = direction
        reverse_directions = []
        for step in self.stencil_steps:
            reverse_key = tuple(int(-value) for value in step)
            if reverse_key not in step_lookup:
                raise ValueError(
                    f"stencil step {tuple(step)} has no exact negative partner"
                )
            reverse_directions.append(step_lookup[reverse_key])
        self.reverse_directions = np.asarray(reverse_directions, dtype=int)
        if self.uses_legacy_four_stencil:
            self.link_direction_labels = LINK_DIRECTION_LABELS
        else:
            self.link_direction_labels = tuple(
                f"({step[0]:+d},{step[1]:+d})" for step in self.stencil_steps
            )

        reciprocal_rows = 2.0 * np.pi * np.linalg.inv(self.cell).T
        reduced_steps = self.stencil_steps / self.mesh[None, :]
        self.b_cartesian = (reduced_steps @ reciprocal_rows[:2])[:, :2]
        stencil_first_moment = np.einsum(
            "dj,dk->jk", self.derivative_weights, self.b_cartesian
        )
        stencil_first_moment_error = _max_abs(stencil_first_moment - np.eye(2))
        reverse_weights = self.derivative_weights[self.reverse_directions]
        stencil_odd_error = _max_abs(self.derivative_weights + reverse_weights)
        if stencil_first_moment_error > self.stencil_tolerance:
            raise ValueError(
                "derivative_weights fail Cartesian first-moment normalization "
                f"(max error {stencil_first_moment_error:.3e})"
            )
        if stencil_odd_error > self.stencil_tolerance:
            raise ValueError(
                "derivative_weights are not odd under step reversal "
                f"(max error {stencil_odd_error:.3e})"
            )

        self._mesh_indices, mesh_shift = self._validate_mesh_and_neighbors()
        z_hermitian_error = _max_abs(self.z - adjoint(self.z))
        spin_hermitian_error = _max_abs(self.spin_tensor - adjoint(self.spin_tensor))
        if z_hermitian_error > self.hermitian_tolerance:
            raise ValueError(
                f"Z is not Hermitian (max error {z_hermitian_error:.3e})"
            )
        if spin_hermitian_error > self.hermitian_tolerance:
            raise ValueError(
                "spin_tensor is not Hermitian "
                f"(max error {spin_hermitian_error:.3e})"
            )

        link_reverse_error = 0.0
        for direction in range(self.nd):
            reverse_direction = self.reverse_directions[direction]
            reverse = self.links[
                self.neighbors[:, direction], reverse_direction
            ]
            link_reverse_error = max(
                link_reverse_error,
                _max_abs(self.links[:, direction] - adjoint(reverse)),
            )
        if link_reverse_error > self.link_tolerance:
            raise ValueError(
                "links are not reciprocal adjoints "
                f"(max error {link_reverse_error:.3e})"
            )

        self.operator_layer_sum = (
            self.spin_tensor[:, 0]
            - self.spin_tensor[:, 1]
            - self.spin_tensor[:, 2]
        )
        operator_layer_sum_error = _max_abs(self.operator_layer_sum)
        if operator_layer_sum_error > self.layer_sum_tolerance:
            raise ValueError(
                "independent spin scopes do not satisfy total = lower + upper "
                f"(max error {operator_layer_sum_error:.3e})"
            )

        self.h = np.zeros((self.nk, self.nb, self.nb), dtype=complex)
        diagonal = np.arange(self.nb)
        self.h[:, diagonal, diagonal] = self.energies

        self.validation = {
            "mesh_shift_reduced": mesh_shift,
            "z_hermitian_max_abs": z_hermitian_error,
            "spin_hermitian_max_abs": spin_hermitian_error,
            "link_reverse_max_abs": link_reverse_error,
            "stencil_first_moment_max_abs": stencil_first_moment_error,
            "stencil_odd_max_abs": stencil_odd_error,
            "operator_layer_sum_max_abs": operator_layer_sum_error,
            "area_A2": self.area_A2,
        }

    def _validate_mesh_and_neighbors(self) -> tuple[np.ndarray, np.ndarray]:
        """Return integer grid coordinates after validating a shifted mesh."""

        base = self.k[0]
        grid_axes = []
        shifts = []
        for axis, count in enumerate(self.mesh):
            scaled = count * (self.k[:, axis] - base[axis])
            nearest = np.rint(scaled)
            if _max_abs(scaled - nearest) > self.mesh_tolerance:
                raise ValueError(
                    "k is not a uniform shifted mesh in reduced coordinates"
                )
            indices = np.mod(nearest.astype(int), count)
            if np.unique(indices).size != count:
                raise ValueError(f"k mesh does not span all {count} points on axis {axis}")
            grid_axes.append(indices)
            shifts.append(float(np.mod(base[axis], 1.0 / count)))

        z_delta = self.k[:, 2] - base[2]
        if _max_abs(z_delta - np.rint(z_delta)) > self.mesh_tolerance:
            raise ValueError("k must be a two-dimensional mesh with fixed reduced kz")

        grid = np.column_stack(grid_axes)
        if np.unique(grid, axis=0).shape[0] != self.nk:
            raise ValueError("k does not contain each 2D mesh point exactly once")

        for direction, step in enumerate(self.stencil_steps):
            expected = grid.copy()
            expected = np.mod(expected + step[None, :], self.mesh[None, :])
            actual = grid[self.neighbors[:, direction]]
            if not np.array_equal(actual, expected):
                raise ValueError(
                    f"neighbors direction {self.link_direction_labels[direction]} "
                    "does not follow mesh geometry"
                )

            reverse_direction = self.reverse_directions[direction]
            reverse = self.neighbors[
                self.neighbors[:, direction], reverse_direction
            ]
            if not np.array_equal(reverse, np.arange(self.nk)):
                raise ValueError(
                    f"neighbors direction {self.link_direction_labels[direction]} "
                    "does not invert with its reverse"
                )

            target = np.zeros(3)
            target[:2] = step / self.mesh
            residual = self.k[self.neighbors[:, direction]] - self.k - target
            residual -= np.rint(residual)
            if _max_abs(residual) > self.mesh_tolerance:
                raise ValueError(
                    f"neighbors direction {self.link_direction_labels[direction]} "
                    "has an inconsistent reduced-k displacement"
                )
        return grid, np.asarray(shifts)

    def _validated_h(self, h: Optional[np.ndarray]) -> np.ndarray:
        if h is None:
            return self.h
        array = np.asarray(h, dtype=complex)
        _finite("h", array)
        if array.shape != (self.nk, self.nb, self.nb):
            raise ValueError(f"h must have shape ({self.nk}, {self.nb}, {self.nb})")
        error = _max_abs(array - adjoint(array))
        if error > self.hermitian_tolerance:
            raise ValueError(f"h is not Hermitian (max error {error:.3e})")
        return array

    def _chi_core(self, h: np.ndarray) -> dict[str, object]:
        energies, eigenvectors = np.linalg.eigh(h)
        occupied = eigenvectors[:, :, : self.nocc]
        empty = eigenvectors[:, :, self.nocc :]

        derivative = np.zeros((self.nk, 2, self.nb, self.nocc), dtype=complex)
        singular_by_direction = np.empty(self.nd)
        for direction in range(self.nd):
            transported = (
                self.links[:, direction]
                @ occupied[self.neighbors[:, direction]]
            )
            occupied_overlap = adjoint(occupied) @ transported
            singular_values = np.linalg.svd(occupied_overlap, compute_uv=False)
            singular_by_direction[direction] = float(np.min(singular_values))
            if singular_by_direction[direction] < self.singular_tolerance:
                raise ValueError(
                    "occupied neighboring spaces are singular in direction "
                    f"{self.link_direction_labels[direction]} "
                    f"(minimum singular value {singular_by_direction[direction]:.3e})"
                )
            parallel = transported @ np.linalg.inv(occupied_overlap)
            coefficients = self.derivative_weights[direction]
            derivative += parallel[:, None, :, :] * coefficients[None, :, None, None]

        indirect_gap = float(
            np.min(energies[:, self.nocc])
            - np.max(energies[:, self.nocc - 1])
        )
        direct_gaps = energies[:, self.nocc] - energies[:, self.nocc - 1]
        min_direct_gap = float(np.min(direct_gaps))
        if indirect_gap <= self.gap_tolerance:
            raise ValueError(
                "global insulating gap closed; frozen-potential insulating response "
                f"is invalid (gap {indirect_gap:.6g} eV)"
            )

        dipole = adjoint(empty)[:, None] @ (1j * derivative)
        dipole = dipole.swapaxes(-1, -2)  # (k, E, occupied, empty)
        spin_cross = np.einsum(
            "kbo,ksibc,kce->ksioe",
            occupied.conj(),
            self.spin_tensor,
            empty,
            optimize=True,
        )
        denominators = (
            energies[:, : self.nocc, None]
            - energies[:, None, self.nocc :]
        )
        products = (
            spin_cross[:, :, :, None, :, :]
            * dipole[:, None, None, :, :, :]
        )
        chi_k = np.sum(
            2.0
            * products.real
            / denominators[:, None, None, None, :, :],
            axis=(-2, -1),
        )
        chi = np.mean(chi_k, axis=0) * self.factor
        layer_sum = chi[0] - chi[1] - chi[2]
        return {
            "chi": chi,
            "gap_eV": indirect_gap,
            "min_direct_gap_eV": min_direct_gap,
            "min_neighbor_singular": float(np.min(singular_by_direction)),
            "min_neighbor_singular_by_direction": singular_by_direction,
            "layer_sum": layer_sum,
            "layer_sum_max_abs": _max_abs(layer_sum),
        }

    def chi(self, h: Optional[np.ndarray] = None) -> dict[str, object]:
        """Return electric spin susceptibility and numerical diagnostics."""

        return self._chi_core(self._validated_h(h))

    def response(
        self, h: Optional[np.ndarray] = None, b: float = 0.002
    ) -> dict[str, object]:
        """Return the central Zeeman derivative ``C[scope, spin, E]``.

        ``b`` is the energy coefficient (eV) multiplying the dimensionless
        cached vertex ``Z``.  Positive and negative probes are evaluated with
        unchanged links and vertices, matching the legacy frozen-potential
        convention.
        """

        if not np.isfinite(b) or b <= 0.0:
            raise ValueError("b must be a positive finite energy in eV")
        base_h = self._validated_h(h)
        zero = self._chi_core(base_h)
        plus = self._chi_core(base_h + b * self.z)
        minus = self._chi_core(base_h - b * self.z)

        chi_plus = plus["chi"]
        chi_minus = minus["chi"]
        c_tensor = (chi_plus - chi_minus) / (2.0 * b)
        even_chi = (chi_plus + chi_minus) / 2.0
        b_even_residual = even_chi - zero["chi"]
        layer_sum = c_tensor[0] - c_tensor[1] - c_tensor[2]
        even_layer_sum = even_chi[0] - even_chi[1] - even_chi[2]

        gaps = np.asarray([minus["gap_eV"], zero["gap_eV"], plus["gap_eV"]])
        direct_gaps = np.asarray(
            [
                minus["min_direct_gap_eV"],
                zero["min_direct_gap_eV"],
                plus["min_direct_gap_eV"],
            ]
        )
        singulars = np.asarray(
            [
                minus["min_neighbor_singular"],
                zero["min_neighbor_singular"],
                plus["min_neighbor_singular"],
            ]
        )
        singulars_by_direction = np.stack(
            [
                minus["min_neighbor_singular_by_direction"],
                zero["min_neighbor_singular_by_direction"],
                plus["min_neighbor_singular_by_direction"],
            ]
        )

        return {
            "C": c_tensor,
            "zero_chi": zero["chi"],
            "chi_plus": chi_plus,
            "chi_minus": chi_minus,
            "even_chi": even_chi,
            "b_even_residual": b_even_residual,
            "b_even_residual_max_abs": _max_abs(b_even_residual),
            "gap_eV": float(np.min(gaps)),
            "gaps_eV_minus_zero_plus": gaps,
            "min_direct_gap_eV": float(np.min(direct_gaps)),
            "direct_gaps_eV_minus_zero_plus": direct_gaps,
            "min_neighbor_singular": float(np.min(singulars)),
            "neighbor_singular_minus_zero_plus": singulars,
            "neighbor_singular_by_direction_minus_zero_plus": singulars_by_direction,
            "layer_sum": layer_sum,
            "layer_sum_max_abs": _max_abs(layer_sum),
            "zero_layer_sum_max_abs": zero["layer_sum_max_abs"],
            "even_layer_sum_max_abs": _max_abs(even_layer_sum),
            "operator_layer_sum_max_abs": self.validation[
                "operator_layer_sum_max_abs"
            ],
            "b_eV": float(b),
            "units": C_UNITS,
            "scope_labels": SCOPE_LABELS,
            "spin_labels": SPIN_LABELS,
            "electric_field_labels": E_LABELS,
            "link_direction_labels": self.link_direction_labels,
            "stencil_steps": self.stencil_steps.copy(),
            "derivative_weights_A": self.derivative_weights.copy(),
            "stencil_first_moment_max_abs": self.validation[
                "stencil_first_moment_max_abs"
            ],
            "stencil_odd_max_abs": self.validation["stencil_odd_max_abs"],
            "probe_convention": "[chi(H+bZ)-chi(H-bZ)]/(2b)",
        }


# A small compatibility convenience for analysis scripts descended from the
# legacy scalar module.  The API and returned tensor are the new ones.
Response = TensorResponse


__all__ = [
    "C_UNITS",
    "E_LABELS",
    "LINK_DIRECTION_LABELS",
    "Response",
    "SCOPE_LABELS",
    "SPIN_LABELS",
    "TensorResponse",
    "adjoint",
]
