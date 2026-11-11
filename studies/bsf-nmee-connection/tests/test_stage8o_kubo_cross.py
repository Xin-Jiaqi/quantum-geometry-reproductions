"""Stage8o regression tests for the corrected Kubo/projector cross-check.

The synthetic 2H-AB-type bilayer (SlidingBilayer) has TR+inversion symmetry,
so its static response iTr P[dP/dl, dP/dk] is EXACTLY zero for every spin
source axis at every k. This physical zero is used as the strongest
available local consistency test of the two independent Kubo routes
(resolvent vs projector finite differences): both must return ~0 and must
agree with each other. A nonzero benchmark was attempted (parity-broken
layer potential) but the toy layer orbital ordering made the perturbation
non-layer-diagonal; this limitation is recorded, not papered over.
"""
import unittest
import numpy as np
from dataclasses import replace

from bsf_nmee.sliding_tight_binding import SlidingBilayer, SlidingSpec
from stage8o_kubo_cross import (resolvent_gamma_k, projector_gamma_k,
                                dM_covariant)


class Stage8oKuboCrossTests(unittest.TestCase):
    def _model(self):
        return SlidingBilayer(replace(SlidingSpec(),
                                      shift=(.4 * .23, .4 * .41),
                                      tunneling_eV=.005))

    def _gather(self, model, k, dk, B):
        rec = model.model().reciprocal
        H = model.matrices(np.atleast_2d(k))
        Hm = H[0] if H.ndim == 3 else H
        e, _ = np.linalg.eigh(Hm)
        Sz = model.spins["total_z"]
        Sz = Sz[0] if Sz.ndim == 3 else Sz
        Sx = model.spins["total_x"]
        Sx = Sx[0] if Sx.ndim == 3 else Sx
        kp = np.atleast_2d(k) + np.array([dk, 0.])
        km = np.atleast_2d(k) - np.array([dk, 0.])
        Mp = model.matrices(kp)
        Mm = model.matrices(km)
        Mp = Mp[0] if Mp.ndim == 3 else Mp
        Mm = Mm[0] if Mm.ndim == 3 else Mm
        g = model.mono.g_factor
        Hp = np.diag(e).astype(complex) - g * 5.7884e-5 * B * Sz
        eB, U = np.linalg.eigh(Hp)
        S2_zb = U.T.conj() @ (Sx / 2.0) @ U
        dM = dM_covariant(U, Mp, Mm, dk)
        return eB, S2_zb, dM, Sz, Sx, e, Mp, Mm

    def test_zero_response_all_spin_axes_and_resolvent_projector_agree(self):
        model = self._model()
        rec = model.model().reciprocal
        k = (np.array([[.19, .07]]) @ rec)[0]
        dk = .005
        for axis_name, src in (("x", "total_x"), ("y", "total_y"),
                               ("z", "total_z")):
            e, S2, dM, Sz, Sx, e0, Mp, Mm = self._gather(model, k, dk, .0)
            g_res, gmin = resolvent_gamma_k(e, S2, dM, 4)
            self.assertLess(abs(g_res), 1e-10,
                            f"resolvent not zero for source {axis_name}")
            g_proj, _ = projector_gamma_k(
                e0, Sz, model.spins[src][0] if model.spins[src].ndim == 3
                else model.spins[src], Mp, Mm, 4, dk, 1e-4, 0.0,
                mu_b=5.7884e-5, g=model.mono.g_factor)
            self.assertLess(abs(g_proj), 1e-8,
                            f"projector not zero for source {axis_name}")
            self.assertLess(abs(g_res - g_proj),
                            1e-6 + 1e-4 * max(abs(g_res), abs(g_proj)),
                            "resolvent and projector disagree")


if __name__ == "__main__":
    unittest.main()
