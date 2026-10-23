"""Stage8p: nonzero real-link-compatible projector-vs-resolvent validation."""
import importlib.util
import sys
from pathlib import Path
import unittest
import numpy as np

BASE=Path(__file__).resolve().parents[1]/"hpc-stage8i"
sys.path.insert(0,str(BASE))
import stage8p_nonzero_benchmark as model
import stage8p_projector_kubo_core as core


def make_band_input(k):
    """Synthetic 2-band band-basis H0 and spin-Pauli operator matrices."""
    E,U=np.linalg.eigh(model.ham(k,0,0))
    return E,U.conj().T@model.SIGX@U,U.conj().T@model.SIGZ@U,U


def make_state(k,dk=.002):
    ec,sxc,szc,Uc=make_band_input(k)
    ep,sxp,szp,Up=make_band_input(k+dk)
    em,sxm,szm,Um=make_band_input(k-dk)
    return ((ec,sxc,szc),(ep,sxp,szp),(em,sxm,szm),
            (Uc.conj().T@Up,Uc.conj().T@Um))


class ProjectorKuboInterbandTests(unittest.TestCase):
    def test_h_source_conventions_are_spin_pauli_not_extra_g(self):
        E,sx,sz,_=make_band_input(.4)
        H=core.h_source(E,sx,sz,.002,.25,magnetic_prefactor=.4)
        expected=np.diag(E)+.001*sx-.1*sz
        self.assertTrue(np.allclose(H,expected))

    def test_nonzero_analytic_two_routes_same_sign_and_magnitude(self):
        k=.4;dk=.002;dl=2e-4;dB=.01
        center,plus,minus,links=make_state(k,dk)
        F0=core.compare_local_three_way(center,plus,minus,links,
                                        B=0,occ=1,dk=dk,dl=dl,magnetic_prefactor=.4)
        self.assertAlmostEqual(F0["F_resolvent"],model.analytic_F(k,0),places=6)
        self.assertAlmostEqual(F0["F_projector_fd"],model.analytic_F(k,0),places=6)
        p=core.compare_local_three_way(center,plus,minus,links,
                                       B=dB,occ=1,dk=dk,dl=dl,magnetic_prefactor=.4)
        m=core.compare_local_three_way(center,plus,minus,links,
                                       B=-dB,occ=1,dk=dk,dl=dl,magnetic_prefactor=.4)
        for name in ("F_resolvent","F_projector_fd"):
            gamma=(p[name]-m[name])/(2*dB)
            self.assertGreater(gamma,.01)
            self.assertAlmostEqual(gamma,model.analytic_Gamma(k),delta=3e-5)
        self.assertLess(F0["MMN_subspace_leakage"]["max_link_nonunitarity_fro"],1e-12)
        self.assertLess(F0["MMN_subspace_leakage"]["max_neighbour_projector_idempotency_error"],1e-12)

    def test_neighbor_basis_rephasing_has_no_effect(self):
        k=.43;dk=.003;center,plus,minus,links=make_state(k,dk)
        old=core.compare_local_three_way(center,plus,minus,links,
                                         B=.01,occ=1,dk=dk,dl=2e-4,magnetic_prefactor=.4)
        dp=np.diag(np.exp(1j*np.array([.4,1.7])))
        dm=np.diag(np.exp(1j*np.array([-.5,2.0])))
        def rephase(t,D):
            e,sx,sz=t
            return (e,D.conj().T@sx@D,D.conj().T@sz@D)
        new=core.compare_local_three_way(center,rephase(plus,dp),rephase(minus,dm),
                                         (links[0]@dp,links[1]@dm),
                                         B=.01,occ=1,dk=dk,dl=2e-4,magnetic_prefactor=.4)
        self.assertAlmostEqual(old["F_resolvent"],new["F_resolvent"],places=11)
        self.assertAlmostEqual(old["F_projector_fd"],new["F_projector_fd"],places=11)

    def test_nonunitary_truncated_links_are_quantified_not_hidden(self):
        center,plus,minus,links=make_state(.4)
        bad=(links[0]*.8,links[1])
        x=core.compare_local_three_way(center,plus,minus,bad,
                                       B=0,occ=1,dk=.002,dl=.0002,magnetic_prefactor=.4)
        self.assertGreater(x["MMN_subspace_leakage"]["max_link_nonunitarity_fro"],.1)
        self.assertGreater(x["MMN_subspace_leakage"]["max_neighbour_projector_idempotency_error"],.01)

    def test_invalid_gap_and_invalid_dk_fail(self):
        center,plus,minus,links=make_state(.4)
        with self.assertRaisesRegex(ValueError,"step"):
            core.compare_local_three_way(center,plus,minus,links,
                                         B=0,occ=1,dk=0,dl=.0002,magnetic_prefactor=.4)
        with self.assertRaisesRegex(ValueError,"occupied and unoccupied"):
            core.frame_and_projector(np.diag([0.,1.]),2)

if __name__=="__main__":unittest.main()
