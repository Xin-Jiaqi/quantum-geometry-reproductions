"""Stage8l: synthetic pseudo occupied-empty spin operator tests, NO actual VASP."""
import importlib.util
import sys
from pathlib import Path
from unittest.mock import patch
import unittest
import numpy as np

HERE=Path(__file__).resolve().parents[1]/"hpc-stage8i"
sys.path.insert(0,str(HERE))
spec=importlib.util.spec_from_file_location("stage8l_pseudo_interband_gate",
                                             HERE/"stage8l_pseudo_interband_gate.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class Stage8lPseudoInterbandTests(unittest.TestCase):
    def test_spin_flip_pair_nonzero_sx_and_sy(self):
        up=np.array([1+0j]);zero=np.array([0+0j])
        v=m.source_vertices((up,zero),(zero,up))
        self.assertAlmostEqual(abs(v[0]),.5)
        self.assertAlmostEqual(abs(v[1]),.5)
        self.assertAlmostEqual(abs(v[2]),0)
    def test_spin_conserving_same_up_has_zero_transverse_vertex(self):
        up=np.array([1+0j]);zero=np.array([0+0j])
        v=m.source_vertices((up,zero),(up,zero))
        self.assertAlmostEqual(abs(v[0]),0)
        self.assertAlmostEqual(abs(v[1]),0)
        self.assertAlmostEqual(v[2].real,.5)
    def test_interband_amplitude_independent_of_each_band_phase(self):
        a=(np.array([1+2j,3j]),np.array([2-1j,1j]))
        b=(np.array([2j,-1+1j]),np.array([1+3j,.5j]))
        x=m.source_vertices(a,b)
        y=m.source_vertices(tuple(z*np.exp(1j*.7) for z in a),
                            tuple(z*np.exp(1j*1.6) for z in b))
        self.assertTrue(np.allclose(np.abs(x),np.abs(y),atol=1e-12))
    def test_zero_job_screen_uses_exact_grid_and_one_based_occupied_empty(self):
        coords=list(m.TARGETS.values())
        kpts=np.array(coords,float)
        energies=[[[(i-37)*.25, 1 if i<36 else 0] for i in range(50)]
                  for _ in range(4)]
        def km(w,h,index):
            arr=np.zeros((50,3))
            arr[:,0]=[v[0] for v in energies[index]]
            return 2,kpts[index],arr
        def coeff(w,h,index,b0):
            data=np.array([1+0j,0j],dtype=np.complex64) if b0<36 else np.array([0j,1+0j],dtype=np.complex64)
            return data,2
        with patch.object(m,"read_header",return_value={"nk":4,"nb":50}), \
             patch.object(m,"parse_eigenval",return_value=(kpts.tolist(),energies)), \
             patch.object(m,"kpoint_meta",side_effect=km), \
             patch.object(m,"band_coeffs",side_effect=coeff):
            ans=m.screen("WAVECAR","EIGENVAL")
        self.assertEqual(ans["new_PBS_jobs"],0)
        self.assertIsNone(ans["material_static_NMEE"])
        self.assertEqual(len(ans["k_points"]),4)
        for p in ans["k_points"].values():
            self.assertEqual(p["candidate_pair_count"],32)
            self.assertTrue(all(x["occupied_band_1based"]<=36
                                and x["empty_band_1based"]>=37 for x in p["pairs"]))
            self.assertGreater(p["pseudo_crossgap_transverse_spin_block_frobenius_norm"],0)
    def test_bad_pseudo_norm_is_fail_closed(self):
        zero=np.array([0j])
        with self.assertRaisesRegex(ValueError,"norm"):
            m.source_vertices((zero,zero),(zero,zero))

if __name__=="__main__":unittest.main()
