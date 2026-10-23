"""Stage8o mock-data tests of real crystal fractional symmetry SPN validator.

No spglib, WAVECAR, source files or VASP are needed for synthetic tests.
"""
import importlib.util
import sys
import numpy as np
from pathlib import Path
import unittest

path=Path(__file__).resolve().parents[1]/"hpc-stage8i"/"stage8o_actual_spn_symmetry_screen.py"
sys.path.insert(0,str(path.parent))
spec=importlib.util.spec_from_file_location("stage8o_actual_spn_symmetry_screen",path)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fake_hex_lattice(a=3.319):
    return np.array([[a,0,0],[-a/2,np.sqrt(3)*a/2,0],[0,0,25.]])

def full_grid(n=3):
    return np.array([(i/n,j/n) for i in range(n) for j in range(n)])

def spin4(nk):
    S=np.zeros((nk,4,4,3),dtype=complex)
    for k in range(nk):
        S[k,0,1,0]=S[k,1,0,0]=1
        S[k,2,3,0]=S[k,3,2,0]=1
        S[k,0,1,1]=-1j;S[k,1,0,1]=1j
        S[k,2,3,1]=-1j;S[k,3,2,1]=1j
        S[k,:,:,2]=np.diag([1,-1,1,-1])
    return S


class Stage8oRealSymmetryTests(unittest.TestCase):
    def test_actual_hex_mirror_maps_all_nine_k_without_cartesian_rounding(self):
        R=np.eye(3,dtype=int)
        R[:2,:2]=np.array([[0,1],[1,0]])
        mapping=m.pairs_from_fractional_symmetry(full_grid(),R,n1=3,n2=3)
        self.assertEqual(len(mapping),9)
        self.assertEqual(len(set(j for _,j in mapping)),9)
        self.assertEqual(mapping[1],(1,3))

    def test_glide_full_mesh_is_identity_in_k_and_requires_sewing(self):
        R=np.diag([1,1,-1])
        pairs=m.pairs_from_fractional_symmetry(full_grid(),R,n1=3,n2=3)
        self.assertEqual(pairs,[(i,i) for i in range(9)])
        E=np.tile(np.array([0.,1.,2.,3.]),(9,1))
        data=m.symmetry_summary(full_grid(),E,spin4(9),fake_hex_lattice(),
                                R,[.5,.5,0.1],n1=3,n2=3)
        self.assertTrue(data["is_horizontal_glide"])
        self.assertEqual(data["mesh_pairs_exact"],9)
        self.assertLess(data["max_spin_gram_covariance_residual"],1e-12)
        self.assertLess(data["max_bandwise_spectrum_symmetry_residual_eV"],1e-12)
        self.assertEqual(data["glide_nondegenerate_transverse_spin_check"]["number_of_above_threshold"],0)

    def test_actual_vertical_mirror_rotates_spin_gram_and_enforces_spectrum(self):
        R=np.eye(3,dtype=int)
        R[:2,:2]=np.array([[0,1],[1,0]])
        E=np.tile(np.array([0.,1.,2.,3.]),(9,1))
        data=m.symmetry_summary(full_grid(),E,spin4(9),fake_hex_lattice(),
                                R,[0,0,0],n1=3,n2=3)
        self.assertEqual(data["mesh_pairs_exact"],9)
        self.assertLess(data["max_spin_gram_covariance_residual"],1e-12)
        self.assertLess(data["max_bandwise_spectrum_symmetry_residual_eV"],1e-12)

    def test_injected_spin_gram_symmetry_break_detected(self):
        R=np.diag([1,1,-1])
        E=np.tile(np.array([0.,1.,2.,3.]),(9,1))
        s=spin4(9)
        s[2,:,:,0]+=0.5*np.diag([1.,-1.,1.,-1.])
        data=m.symmetry_summary(full_grid(),E,s,fake_hex_lattice(),R,
                                [.5,.5,0.1],n1=3,n2=3)
        self.assertGreater(data["max_spin_gram_covariance_residual"],.01)

    def test_missing_grid_duplicates_fail_closed(self):
        coords=full_grid()
        coords[1]=coords[0]
        with self.assertRaisesRegex(ValueError,"Duplicate folded"):
            m.reciprocal_grid(coords,3,3)

    def test_non_crystalline_reciprocal_matrix_rejected(self):
        R=np.eye(3,dtype=int)
        R[0,2]=1
        with self.assertRaisesRegex(ValueError,"mixes"):
            m.pairs_from_fractional_symmetry(full_grid(),R,3,3)


if __name__=="__main__":
    unittest.main()
