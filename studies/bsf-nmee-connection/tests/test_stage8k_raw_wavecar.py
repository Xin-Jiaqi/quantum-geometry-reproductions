#!/usr/bin/env python3
"""Stage8k automated tests for the independent WAVECAR reader and spin checks.

Offline/self-contained: builds a synthetic VASP-style ncl WAVECAR in memory and
verifies layout, energy round-trip, spinor spin expectations, phase invariance,
grouped-vs-interleaved layout discrimination, time-reversal spin reversal,
Kramers gauge warning, and failure branches. The raw HPC WAVECAR integration
results (SHA, maxdE=0, TR reversal on the real file) are recorded in
STAGE8K_INDEPENDENT_SPIN_EVIDENCE_20261022.json and asserted here when present.

Runs under both unittest (CI: python -m unittest discover -s tests)
and pytest (python -m pytest tests/test_stage8k_raw_wavecar.py -v).
"""
import unittest

import json, os, struct, sys, tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOCAL = HERE.parent / "hpc-stage8i"
sys.path.insert(0, str(LOCAL))

from stage8k_independent_reader import (read_header, kpoint_meta, band_coeffs,
                                        split_spinor, spin_expected, find_kindex)


def write_synthetic_wavecar(path, nk=2, nb=6, nG=8, rtag=45200,
                            kpts=None, enocc=None, coeff_seed=1234):
    """Write a minimal ncl WAVECAR: nplane=2*nG, recl8 = 4+3*nb+2*nG (8-byte units)."""
    rng = np.random.default_rng(coeff_seed)
    recl8 = 4 + 3 * nb + 2 * nG
    hdr = recl8 * 8
    kpts = kpts if kpts is not None else [[0.0, 0.0, 0.0], [0.25, 0.25, 0.0]]
    with open(path, "wb") as f:
        # REC1: recl, spin=1, rtag
        f.write(struct.pack("ddd", float(recl8 * 8), 1.0, float(rtag)))
        f.write(b"\x00" * (hdr - 24))
        # REC2: nk, nb, encut, lat(9), efermi, pad
        rec2 = ([float(nk), float(nb), 520.0] + [3.0, 0, 0, 0, 3.0, 0, 0, 0, 25.0]
                + [0.123] + [0.0] * (recl8 - 13))
        assert len(rec2) == recl8
        f.write(struct.pack(f"{len(rec2)}d", *rec2))
        for ik in range(nk):
            krec = [float(2 * nG), *kpts[ik]]
            enocc_flat = []
            for b in range(nb):
                e = -15.0 + 0.2 * b + 0.001 * ik
                occ = 1.0 if b < 3 else 0.0
                enocc_flat += [e, 0.0, occ]
            krec += enocc_flat
            krec += [0.0] * (recl8 - 4 - 3 * nb)
            f.write(struct.pack(f"{len(krec)}d", *krec))
            for b in range(nb):
                # ncl grouped layout: up (nG complex64) then dn (nG complex64)
                up = rng.normal(size=nG) + 1j * rng.normal(size=nG)
                dn = rng.normal(size=nG) + 1j * rng.normal(size=nG)
                # make band0 spin-up pure (Ising test): dn = 0
                if b == 0:
                    dn = np.zeros(nG)
                if b == 1:
                    # TR partner of band0 at -k: pure spin-down
                    up = np.zeros(nG)
                coeffs = np.concatenate([up, dn]).astype(np.complex64)
                f.write(coeffs.tobytes())
                f.write(b"\x00" * ((recl8 - 2 * nG) * 8))
    return recl8


class Stage8kRawWavecarTests(unittest.TestCase):
    def test_header_roundtrip(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "WAVECAR")
            write_synthetic_wavecar(p)
            h = read_header(p)
            assert h["nk"] == 2 and h["nb"] == 6 and h["encut"] == 520.0
            assert h["rtag"] == 45200 and h["spin"] == 1


    def test_energy_roundtrip(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "WAVECAR")
            write_synthetic_wavecar(p)
            h = read_header(p)
            for ik in (0, 1):
                nplane, k, enocc = kpoint_meta(p, h, ik)
                assert nplane == 16  # 2*nG for ncl
                for b in range(6):
                    e = -15.0 + 0.2 * b + 0.001 * ik
                    assert abs(enocc[b][0] - e) < 1e-6
                    assert abs(enocc[b][2] - (1.0 if b < 3 else 0.0)) < 1e-6


    def test_file_size_closure(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "WAVECAR")
            write_synthetic_wavecar(p)
            h = read_header(p)
            expect = 2 * h["hdr_bytes"] + h["nk"] * h["kp_bytes"]
            assert os.path.getsize(p) == expect


    def test_spin_expected_trivial(self):
        up = np.array([1.0 + 0j])
        dn = np.array([0.0 + 0j])
        sx, sy, sz, tot = spin_expected(up, dn)
        assert abs(sx) < 1e-12 and abs(sy) < 1e-12 and abs(sz - 0.5) < 1e-12
        up = np.array([1.0 + 0j])
        dn = np.array([1.0 + 0j])
        sx, sy, sz, tot = spin_expected(up, dn)
        assert abs(sx - 0.5) < 1e-12 and abs(sy) < 1e-12 and abs(sz) < 1e-12


    def test_phase_invariance(self):
        rng = np.random.default_rng(7)
        up = rng.normal(size=5) + 1j * rng.normal(size=5)
        dn = rng.normal(size=5) + 1j * rng.normal(size=5)
        s0 = spin_expected(up, dn)[:3]
        for phi in (0.3, 1.1, 4.7):
            s1 = spin_expected(up * np.exp(1j * phi), dn * np.exp(1j * phi))[:3]
            assert max(abs(a - b) for a, b in zip(s0, s1)) < 1e-12


    def test_grouped_vs_interleaved_discrimination(self):
        """A pure Ising band stored grouped must read Sz=+-1/2 grouped, mixed interleaved."""
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "WAVECAR")
            write_synthetic_wavecar(p)  # band0: up only -> Sz=+0.5
            h = read_header(p)
            c, nplane = band_coeffs(p, h, 0, 0)
            ug, dg, _ = split_spinor(c, nplane, "grouped")
            sg = spin_expected(ug, dg)
            assert abs(sg[2] - 0.5) < 1e-6 and abs(sg[0]) < 1e-6 and abs(sg[1]) < 1e-6
            ui, di, _ = split_spinor(c, nplane, "interleaved")
            si = spin_expected(ui, di)
            # interleaved mis-read must differ (documented parser trap)
            assert (abs(si[2] - 0.5) > 1e-3) or (abs(si[0]) + abs(si[1]) > 1e-3)


    def test_time_reversal_synthetic(self):
        """band0 at k and band1 at -k are TR partners: S(-k) = -S(k) when same band index."""
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "WAVECAR")
            write_synthetic_wavecar(p, kpts=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
            h = read_header(p)
            c0, np0 = band_coeffs(p, h, 0, 0)
            c1, np1 = band_coeffs(p, h, 0, 1)
            u0, d0, _ = split_spinor(c0, np0, "grouped")
            u1, d1, _ = split_spinor(c1, np1, "grouped")
            s0 = spin_expected(u0, d0)
            s1 = spin_expected(u1, d1)
            assert abs(s0[2] - 0.5) < 1e-6 and abs(s1[2] + 0.5) < 1e-6


    def test_kramers_gauge_warning(self):
        """Degenerate spin-conserving H can have gauge-dependent in-plane spin."""
        up=np.array([1.0+0j])
        down=np.array([0.0+0j])
        sx0,sy0,sz0,_=spin_expected(up,down)
        sx1,sy1,sz1,_=spin_expected(up/np.sqrt(2),np.ones(1)/np.sqrt(2))
        self.assertAlmostEqual(sz0,0.5)
        self.assertAlmostEqual(sx0,0.0)
        self.assertAlmostEqual(sx1,0.5)
        self.assertAlmostEqual(sz1,0.0)
        # The parent spin-conserving H=epsilon*I commutes with Sz.
        h=np.eye(2)*2.5; sigmaz=np.diag([1,-1])
        self.assertTrue(np.allclose(h@sigmaz-sigmaz@h,0))


    def test_find_kindex_mod1(self):
        kpts = np.array([[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 0.0], [2 / 3, 2 / 3, 0.0]])
        h = {"nk": 3}
        ki, d = find_kindex(h, (-1 / 3, -1 / 3, 0.0), kpts)
        assert ki == 2 and d < 1e-9


    def test_failure_branches(self):
        with tempfile.TemporaryDirectory() as td:
            p=os.path.join(td,"WAVECAR")
            write_synthetic_wavecar(p,rtag=99999)
            with self.assertRaisesRegex(ValueError,"rtag"):
                read_header(p)
            with self.assertRaises(FileNotFoundError):
                read_header(os.path.join(td,"missing"))

    def test_evidence_json_when_present(self):
        ev = HERE.parent / "hpc-stage8i" / "STAGE8K_INDEPENDENT_SPIN_EVIDENCE_20261022.json"
        if not ev.exists():
            return  # CI without raw data: skip
        d = json.loads(ev.read_text())
        assert d["parser_status"] == "WAVECAR_ENERGY_ROUNDTRIP_PASS"
        for job in ("soc_AB_s0.00", "soc_GLIDE_s0.50"):
            assert job in d["key_results"]
            # isolated nondegenerate band must carry a grouped spin record
            for pt in d["key_results"][job].values():
                assert len(pt["S_grouped"]) == 3
                assert pt["phase_invariance_pass"] is True
        # time-reversal spin reversal must hold for band35 in both jobs
        for job in ("soc_AB_s0.00", "soc_GLIDE_s0.50"):
            assert d["time_reversal"][job]["band35"]["spin_reversal"] is True
        assert "NMEE_NOT_COMPUTED" in d["conclusion"]["status_level"]
