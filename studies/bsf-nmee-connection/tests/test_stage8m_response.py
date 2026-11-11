#!/usr/bin/env python3
"""Stage8m response-pipeline unit tests: .win parsing, formatted .spn reader,
frozen Hamiltonian builder. Offline (no VASP files)."""
import unittest
import tempfile, os, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "hpc-stage8i"))

from stage8m_response import (read_win_kpoints, read_win_cell,
                              read_spn_formatted, build_hamiltonian, MU_B)


class Stage8mResponseTests(unittest.TestCase):

    def test_win_kpoints_and_cell_parse(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "wannier90.win"
            p.write_text(
                "begin unit_cell_cart\n"
                "  3.3190066 0.0 0.0\n"
                " -1.6595033 2.8743441 0.0\n"
                "  0.0 0.0 27.7171056\n"
                "end unit_cell_cart\n"
                "mp_grid = 21 21 1\n"
                "begin kpoints\n"
                "  0.0 0.0 0.0\n"
                "  0.047619047619 0.0 0.0\n"
                "end kpoints\n")
            k = read_win_kpoints(p)
            c = read_win_cell(p)
            self.assertEqual(k.shape, (2, 2))
            self.assertAlmostEqual(k[1, 0], 0.047619047619)
            self.assertEqual(c.shape, (3, 3))
            self.assertAlmostEqual(c[1, 1], 2.8743441)

    def test_spn_roundtrip(self):
        nb, nk = 4, 2
        rng = np.random.default_rng(7)
        S = np.zeros((nk, nb, nb, 3), complex)
        for ik in range(nk):
            for ip in range(3):
                A = rng.normal(size=(nb, nb)) + 1j * rng.normal(size=(nb, nb))
                H = (A + A.T.conj()) / 2
                np.fill_diagonal(H, np.abs(np.diag(H)))   # real diagonal
                S[ik, :, :, ip] = H
        Ntri = nb * (nb + 1) // 2
        indn, indm = np.tril_indices(nb)
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "wannier90.spn"
            with open(p, "w") as f:
                f.write("Written by TEST\n")
                f.write(f"{nb} {nk}\n")
                for ik in range(nk):
                    A3 = np.zeros((3, nb, nb), complex)
                    A3[:, indn, indm] = S[ik].transpose(2, 0, 1)[:, indn, indm]
                    for z in A3[:, indn, indm].T.ravel():
                        f.write(f"  {z.real:.12e}  {z.imag:.12e}\n")
            got = read_spn_formatted(str(p), nb, nk)
            self.assertLess(np.max(np.abs(got - S)), 1e-10)

    def test_spn_short_row_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "wannier90.spn"
            p.write_text("Written by TEST\n2 1\n  1.0 0.0\n")  # too few rows
            with self.assertRaises(ValueError):
                read_spn_formatted(str(p), 2, 1)

    def test_build_hamiltonian_frozen(self):
        rng = np.random.default_rng(11)
        E = np.sort(rng.normal(size=6))
        Sx = rng.normal(size=(6, 6))
        Sx = (Sx + Sx.T) / 2
        Sz = np.diag(np.linspace(-1, 1, 6))
        H = build_hamiltonian(E, Sx, Sz, 0.0, 0.0, 6)
        self.assertTrue(np.allclose(H, np.diag(E)))  # lam=B=0 limit check
        H2 = build_hamiltonian(E, Sx, Sz, 0.001, 0.25, 6)
        ref = np.diag(E) + 0.001 * 0.5 * Sx - MU_B * 0.25 * Sz
        self.assertLess(np.max(np.abs(H2 - ref)), 1e-14)


if __name__ == "__main__":
    unittest.main()
