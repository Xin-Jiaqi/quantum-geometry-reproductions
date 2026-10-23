"""Stage8f: source-backed sliding-glide response gate and model existence check."""
import unittest
import numpy as np
from bsf_nmee.sliding_tight_binding import (
    SlidingBilayer,SlidingSpec,M_GLOBAL)
from bsf_nmee.valley_kp import static_kernel
from nmee_lab.clean_core.core import MU_B
import importlib.util
from pathlib import Path

FILE=Path(__file__).resolve().parents[1]/"material-screening"/"stage8f_glide_selection.py"
spec=importlib.util.spec_from_file_location("source_glide",FILE)
glide=importlib.util.module_from_spec(spec);spec.loader.exec_module(glide)


class Stage8FGlideAndLayerSpinTests(unittest.TestCase):
    def test_three_paths_equal_and_exact_tensor_ranks(self):
        result=glide.gate()
        np.testing.assert_allclose(result["three_C3_related_arm_lengths_A"],
                                   [3.319006624468/np.sqrt(3)]*3,rtol=1e-9)
        for arm,rows in result["paths"].items():
            self.assertEqual([r["order"] for r in rows],[6,2,4,2,6])
            self.assertEqual([r["total_spin_tensor_dim"] for r in rows],[1,2,0,2,1])
            self.assertEqual([r["layer_odd_tensor_dim"] for r in rows],[1,2,2,2,1])
            self.assertEqual([r["is_layer_exchange_glide"] for r in rows],
                             [False,False,True,False,False])
            self.assertTrue(all(r["mirror_maps_to_one_minus_s"] for r in rows))

    def test_all_path_mirror_pairs_not_just_sample_nodes(self):
        cell,src=glide.read_source()
        for arm in glide.ARMS:
            for s in (.07,.12,.33,.41,.50,.67,.88):
                self.assertTrue(glide.mirror_maps(
                    glide.slide(src,arm,s),
                    glide.slide(src,arm,1-s),cell))

    def test_toy_glide_even_forbidden_odd_allowed_and_finite(self):
        # Independently demonstrated existence at the toy m_z-symmetric
        # reference u=0. NOT at WSe2 structure, not actual material NMEE.
        spec=SlidingSpec(shift=(0.,0.),tunneling_eV=.005)
        obj=SlidingBilayer(spec)
        spin_even=obj.spins["total_x"]
        spin_odd=obj.spins["top_x"]-obj.spins["bottom_x"]
        self.assertLess(np.max(abs(M_GLOBAL@spin_even@M_GLOBAL+spin_even)),1e-12)
        self.assertLess(np.max(abs(M_GLOBAL@spin_odd@M_GLOBAL-spin_odd)),1e-12)
        gauss,w=np.polynomial.legendre.leggauss(12)
        radius=.75
        r=(gauss+1)*radius/2
        rw=w*radius/2
        theta=2*np.pi*(np.arange(36)+.5)/36
        q=(r[:,None,None]*np.stack([np.cos(theta),np.sin(theta)],
                                   axis=-1)[None,:,:]).reshape(-1,2)
        weights=np.repeat(r*rw,36)/(2*np.pi*36)
        centers=np.array([[1/3,2/3],[2/3,1/3]])@obj.model().reciprocal
        value=np.zeros(2)
        for center in centers:
            k=q+center
            h=obj.matrices(k)
            v=np.stack([obj.matrices(k,axis=j) for j in (0,1)])
            self.assertLess(np.max(abs(M_GLOBAL@h@M_GLOBAL-h)),2e-11)
            self.assertLess(np.max(abs(M_GLOBAL@v[0]@M_GLOBAL-v[0])),2e-11)
            ker=static_kernel(h,v,[spin_even,spin_odd],
                     obj.spins["total_z"],spec.mu_eV,
                     g_factor=obj.mono.g_factor,step_T=.25)
            value+=np.array([weights@ker[:,j,0] for j in (0,1)])
        self.assertLess(abs(value[0]),1e-9)
        self.assertLess(value[1],-.015)


if __name__=="__main__":
    unittest.main()
