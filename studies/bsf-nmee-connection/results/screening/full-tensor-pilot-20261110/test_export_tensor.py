"""Small algebraic tests for :mod:`export_tensor`; no VASP files required."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parent))
from export_tensor import adj, masked_wave, spin_blocks, vertices  # noqa: E402


class _FakePawPotential:
    def __init__(self, qij: np.ndarray) -> None:
        self._qij = np.asarray(qij)

    def get_Qij(self) -> np.ndarray:
        return self._qij


class _FakeProjectors:
    def __init__(self) -> None:
        self.element_idx = np.array([0, 0])
        self.qproj = [
            np.array(
                [
                    [1.0, 0.15j, -0.08, 0.04j],
                    [0.06, -0.12j, 0.9, 0.11],
                ],
                dtype=complex,
            )
        ]
        self.cqfak = [np.array([1.1, 0.75])]
        self.crexp = np.array(
            [
                [1.0, np.exp(0.23j)],
                [np.exp(-0.17j), np.exp(0.41j)],
                [np.exp(0.31j), np.exp(-0.29j)],
                [np.exp(-0.37j), np.exp(0.13j)],
            ]
        )
        self.pawpp = [
            _FakePawPotential(
                np.array([[0.24, 0.03 + 0.02j], [0.03 - 0.02j, 0.17]])
            )
        ]


class ExportTensorTests(unittest.TestCase):
    def test_spin_blocks_uses_standard_pauli_y_sign(self) -> None:
        rng = np.random.default_rng(1902)
        wave = rng.normal(size=(2, 4, 3)) + 1j * rng.normal(size=(2, 4, 3))
        projected = rng.normal(size=(2, 4, 3)) + 1j * rng.normal(
            size=(2, 4, 3)
        )
        pauli = np.array(
            [
                [[0, 1], [1, 0]],
                [[0, -1j], [1j, 0]],
                [[1, 0], [0, -1]],
            ],
            dtype=complex,
        )
        expected = np.array(
            [
                sum(
                    pauli[component, left, right]
                    * (adj(wave[left]) @ projected[right])
                    for left in range(2)
                    for right in range(2)
                )
                for component in range(3)
            ]
        )
        np.testing.assert_allclose(spin_blocks(wave, projected), expected)

        plus_y = np.array([1.0, 1.0j], dtype=complex)[:, None, None] / np.sqrt(2)
        np.testing.assert_allclose(
            spin_blocks(plus_y, plus_y)[:, 0, 0], [0.0, 1.0, 0.0], atol=1e-15
        )

    def test_masked_wave_is_periodic_fourier_projection(self) -> None:
        rng = np.random.default_rng(9021)
        wave = rng.normal(size=(2, 5, 3)) + 1j * rng.normal(size=(2, 5, 3))
        g = np.array(
            [
                [0, 0, -2],
                [0, 0, 0],
                [0, 0, 3],
                [1, 0, -1],
                [1, 0, 2],
            ]
        )
        # The two periodic boundaries lie exactly on FFT planes and therefore
        # carry half weight, matching the legacy z-partition convention.
        mask = np.array([0.5, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0, 0.0])

        expected = np.empty_like(wave)
        real_grid = np.arange(len(mask)) / len(mask)
        for xy in np.unique(g[:, :2], axis=0):
            ids = np.flatnonzero(np.all(g[:, :2] == xy, axis=1))
            block = np.empty((len(ids), len(ids)), dtype=complex)
            for row, left in enumerate(ids):
                for column, right in enumerate(ids):
                    block[row, column] = np.mean(
                        mask
                        * np.exp(
                            2j
                            * np.pi
                            * (g[right, 2] - g[left, 2])
                            * real_grid
                        )
                    )
            for spin in range(2):
                expected[spin, ids] = block @ wave[spin, ids]

        actual = masked_wave(wave, g, mask)
        np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=2e-15)
        np.testing.assert_allclose(
            actual + masked_wave(wave, g, 1.0 - mask),
            wave,
            rtol=2e-15,
            atol=2e-15,
        )

    def test_vertices_assigns_each_paw_augmentation_to_its_own_layer(self) -> None:
        rng = np.random.default_rng(711)
        wave = rng.normal(size=(2, 4, 3)) + 1j * rng.normal(size=(2, 4, 3))
        g = np.array([[0, 0, -1], [0, 0, 0], [0, 0, 2], [0, 0, 3]])
        lower_mask = np.array([0.5, 1.0, 1.0, 0.5, 0.0, 0.0])
        masks = [lower_mask, 1.0 - lower_mask]
        projectors = _FakeProjectors()
        atom_lower = np.array([True, False])

        smooth_total = spin_blocks(wave, wave)
        smooth_layers = [
            spin_blocks(wave, masked_wave(wave, g, mask)) for mask in masks
        ]
        augmentations = []
        expected_metric = sum(adj(component) @ component for component in wave)
        for atom in range(2):
            p = (
                projectors.qproj[0]
                * projectors.cqfak[0][:, None]
                * projectors.crexp[None, :, atom]
            )
            beta = np.array([p @ component for component in wave])
            qb = np.array(
                [projectors.pawpp[0].get_Qij() @ component for component in beta]
            )
            augmentations.append(spin_blocks(beta, qb))
            expected_metric += sum(
                adj(component) @ corrected
                for component, corrected in zip(beta, qb)
            )

        expected = np.array(
            [
                smooth_total + augmentations[0] + augmentations[1],
                smooth_layers[0] + augmentations[0],
                smooth_layers[1] + augmentations[1],
            ]
        )
        expected = (expected + adj(expected)) / 2

        actual, metric, hermiticity_error = vertices(
            wave, g, projectors, masks, atom_lower
        )
        np.testing.assert_allclose(actual, expected, rtol=3e-15, atol=3e-15)
        np.testing.assert_allclose(metric, expected_metric, rtol=3e-15, atol=3e-15)
        np.testing.assert_allclose(actual[0], actual[1] + actual[2], atol=5e-15)
        np.testing.assert_allclose(
            actual[1] - smooth_layers[0], augmentations[0], atol=5e-15
        )
        np.testing.assert_allclose(
            actual[2] - smooth_layers[1], augmentations[1], atol=5e-15
        )
        self.assertLess(hermiticity_error, 5e-14)


if __name__ == "__main__":
    unittest.main(verbosity=2)
