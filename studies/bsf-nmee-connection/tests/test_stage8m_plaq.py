#!/usr/bin/env python3
"""Stage8m regressions: DeterminantPlaquette engine must reproduce the
independently validated Stage7i Wilson/Kubo static-response benchmark, and
the wannier90 file readers must fail closed on malformed inputs.

Runs offline (synthetic tight-binding model, no VASP files).
"""
import unittest
import tempfile, os, sys
from pathlib import Path
from dataclasses import replace
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "hpc-stage8i"))

from bsf_nmee.sliding_tight_binding import SlidingBilayer, SlidingSpec
from bsf_nmee.mixed_wilson_geometry import wilson_static_kernel
from nmee_lab.clean_core.core import MU_B

from stage8m_plaq import (DeterminantPlaquette, read_eig, read_mmn_spn,
                          link_matrix, polar_unitary)


def spec_u(u, t):
    return replace(SlidingSpec(), shift=(.23 * u, .41 * u), tunneling_eV=t)


class ModelProvider:
    """Exact occupied frames for the 8-band model (4 occupied bands)."""

    def __init__(self, model, spin_axis="x", nb_occ=4):
        self.model = model
        self.spin = model.spins["total_" + spin_axis]
        self.nb_occ = nb_occ

    def frames(self, kk, lam, B):
        k = np.atleast_2d(np.asarray(kk, float))
        H = self.model.matrices(k)          # (..., 8, 8) with batch first
        hh = H + lam * self.spin - (
            self.model.mono.g_factor * MU_B * B * self.model.spins["total_z"])
        e, v = np.linalg.eigh(hh)
        mu = self.model.spec.mu_eV
        if v.ndim == 3:
            g = min(float(np.min(e[:, 4] - mu)), float(np.min(mu - e[:, 3])))
            occ = v[:, :, :self.nb_occ]
            if occ.shape[0] == 1:
                occ = occ[0]
        else:
            g = min(float(e[4] - mu), float(mu - e[3]))
            occ = v[:, :self.nb_occ]
        return np.ascontiguousarray(occ), g


class Stage8mPlaquetteTests(unittest.TestCase):

    def test_engine_matches_stage7i_wilson_kernel_local(self):
        model = SlidingBilayer(spec_u(.4, .005))
        kc = np.array([[1 / 3, 2 / 3]]) @ model.model().reciprocal + np.array([[.05, .02]])
        prov = ModelProvider(model, "x")
        ref, gap_ref = wilson_static_kernel(model, kc, dk=.005, dl=.000125)
        eng = DeterminantPlaquette(
            lambda k, lam, B: prov.frames(k, lam, B)[0],
            mode="frame", dk=.005, dl=.000125, axis=0)
        gam = eng.Gamma_k(kc[0], 0.0, 0.0, dB=.25)
        self.assertGreater(gap_ref, .003)
        self.assertLess(abs(gam - ref[0]) / abs(ref[0]), 5e-8)

    def test_mmn_chain_matches_frame_for_smooth_overlap(self):
        """mmn mode builds each k edge as the two-step chain
        k1 -> k_c -> k2 from single-step bonds; for a smooth Bloch overlap
        this approximates the direct centred bond to O(dk^2), so it must
        track the exact frame-mode value at small dk."""
        model = SlidingBilayer(spec_u(.4, .005))
        kc = np.array([[1 / 3, 2 / 3]]) @ model.model().reciprocal + np.array([[.05, .02]])
        prov = ModelProvider(model, "x")
        ref, _ = wilson_static_kernel(model, kc, dk=.005, dl=.000125)
        eye = np.eye(8)
        # single-step bonds: k1=(0), k_c=(1), k2=(2); chain needs (k1,kc),(kc,k2)
        mats = {(0, 1): eye, (1, 0): eye, (1, 2): eye, (2, 1): eye,
                (0, 0): eye, (1, 1): eye, (2, 2): eye}

        def ik_of(k_c, sgn):
            return k_c + (1 if sgn > 0 else -1)

        def frames_idx(ik, lam, B):
            kk = kc[0] + np.array([.0025, 0.0]) * (ik - 1)
            return prov.frames(kk, lam, B)[0]

        eng_m = DeterminantPlaquette(
            frames_idx, mode="mmn", links=mats, ik_of=ik_of,
            dk=.005, dl=.000125, axis=0)
        g1 = eng_m.Gamma_k(1, 0.0, 0.0, .25)
        eng_f = DeterminantPlaquette(
            lambda k, lam, B: prov.frames(k, lam, B)[0],
            mode="frame", dk=.005, dl=.000125, axis=0)
        g2 = eng_f.Gamma_k(kc[0], 0.0, 0.0, .25)
        # O(dk^2) chain error; dk = 5e-3 so expect well under 1e-4 relative
        self.assertLess(abs(g1 - g2) / abs(g2), 1e-4)
        self.assertLess(abs(g2 - ref[0]) / abs(ref[0]), 5e-8)

    def test_full_mesh_engine_agrees_with_stage7i_kubo_sign(self):
        vals = {}
        for u in (.4, -.4):
            model = SlidingBilayer(spec_u(u, .005))
            rec = model.model().reciprocal
            n = 6
            base = (np.arange(n) - n / 2) / n
            gx, gy = np.meshgrid(base, base)
            kpts = np.stack([gx.ravel(), gy.ravel()], axis=-1) @ rec
            prov = ModelProvider(model, "x")
            eng = DeterminantPlaquette(
                lambda k, lam, B: prov.frames(k, lam, B)[0],
                mode="frame", dk=.005, dl=.000125, axis=0)
            total = sum(eng.Gamma_k(k, 0.0, 0.0, .25) for k in kpts)
            vals[u] = total / len(kpts)
        self.assertGreater(vals[.4], 0)
        self.assertLess(vals[-.4], 0)
        self.assertLess(abs(vals[.4] + vals[-.4]) / abs(vals[.4]), 5e-6)

    def test_source_step_gap_guard_fails_closed(self):
        model = SlidingBilayer(spec_u(.4, .005))
        kc = np.array([[1 / 3, 2 / 3]]) @ model.model().reciprocal + np.array([[.05, .02]])
        prov = ModelProvider(model, "x")

        def guarded(k, lam, B):
            v, g = prov.frames(k, lam, B)
            if g <= 1e-6:
                raise ValueError("gap closed by source")
            return v

        eng = DeterminantPlaquette(guarded, mode="frame", dk=.005,
                                   dl=.000125, axis=0)
        g = eng.plaquette_phase(kc[0], 1e-4, 0.0)
        self.assertTrue(np.isfinite(g))

        def closed(k, lam, B):
            raise ValueError("gap closed by source")

        bad = DeterminantPlaquette(closed, mode="frame", dk=.005,
                                   dl=.000125, axis=0)
        with self.assertRaises(ValueError):
            bad.plaquette_phase(kc[0], 1e-4, 0.0)

    def test_polar_unitary_singular_value_tracking(self):
        rng = np.random.default_rng(3)
        m = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        u, smin = polar_unitary(m)
        self.assertLess(abs(np.linalg.det(u)) - 1, 1e-10)
        self.assertGreater(smin, 0)

    def test_eig_reader_roundtrip_and_fail_closed(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "x.eig")
            with open(p, "w") as f:
                for k in range(2):
                    for b in range(3):
                        f.write(f"{b+1} {k+1} {-1.0 + 0.1*b + 0.01*k}\n")
            E = read_eig(p)
            self.assertEqual(E.shape, (2, 3))
            self.assertAlmostEqual(E[1, 2], -1.0 + 0.2 + 0.01)
            with open(p, "w") as f:
                f.write("1 1 not_a_number\n")
            with self.assertRaises(ValueError):
                read_eig(p)
            with open(p, "w") as f:
                f.write("")
            with self.assertRaises(ValueError):
                read_eig(p)

    def test_mmn_reader_incomplete_matrix_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "x.mmn")
            with open(p, "w") as f:
                f.write("2 1 1\n")
                f.write("1 1 1 1 1.0 0.0\n")
            with self.assertRaises(ValueError):
                read_mmn_spn(p)

    def test_link_matrix_reverse_conjugate(self):
        rng = np.random.default_rng(5)
        m = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        mats = {(1, 0): m}
        got = link_matrix(mats, 0, 1, 3)
        self.assertLess(np.max(np.abs(got - m.T.conj())), 1e-12)
        with self.assertRaises(ValueError):
            link_matrix(mats, 2, 1, 3)


def sgn_lookup(k, base):
    return np.array([.0025, 0.0]) if np.allclose(k, base) else np.zeros(2)


if __name__ == "__main__":
    unittest.main()
