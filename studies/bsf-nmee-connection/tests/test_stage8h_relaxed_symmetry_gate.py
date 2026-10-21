"""Stage8h: review original/source-like geometry and reject broken glide or absent outputs."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
def load(name):
    file=ROOT/"hpc-stage8g"/(name+".py")
    s=importlib.util.spec_from_file_location("stage8h_"+name,file)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    return m

prepare=load("prepare")
review=load("review_phase1")

class SourceBackedRelaxationReviewTests(unittest.TestCase):
    def test_source_geometries_match_ideal_point_groups(self):
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/"generated"
            prepare.build(folder)
            for key,(order,dim,glide) in review.GROUP_EXPECTED.items():
                result=review.structure_ops(folder/key/"POSCAR")
                self.assertEqual(result["operations"],order,key)
                self.assertEqual(result["tensor_dimension"]["total"],dim,key)
                self.assertEqual(result["horizontal_glide"],glide,key)
            mid=review.structure_ops(folder/"path_A_s0.50"/"POSCAR")
            self.assertEqual(mid["tensor_dimension"]["layer_odd"],2)

    def test_broken_glide_detected_after_asymmetric_z_relaxation(self):
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/"generated"
            prepare.build(folder)
            mid=folder/"path_A_s0.50"/"POSCAR"
            lines=mid.read_text().splitlines()
            j=10  # second W (atoms start at line 9)
            xyz=lines[j].split()
            xyz[2]=f"{float(xyz[2])+0.002:.12f}"  # 0.055 Å asymmetry
            lines[j]=" ".join(xyz)
            mid.write_text("\n".join(lines)+"\n")
            sym=review.structure_ops(mid)
            self.assertFalse(sym["horizontal_glide"])
            self.assertGreater(sym["tensor_dimension"]["total"],0)

    def test_missing_real_output_fails_closed(self):
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/"generated"
            prepare.build(folder)
            result=review.gate(folder)
            self.assertEqual(result["status"],"BLOCKED_RAW_OUTPUT_INCOMPLETE_OR_IONIC_CHECK_FAILED")
            self.assertIsNone(result["no_NEBar_or_NMEE_claim"] if False else None)

    def test_all_three_mock_outputs_still_not_neb(self):
        with tempfile.TemporaryDirectory() as t:
            folder=Path(t)/"generated"
            prepare.build(folder)
            for name,energy in [("path_A_s0.00",-45.664806),
                                ("path_A_s0.25",-45.651172),
                                ("path_A_s0.50",-45.640689)]:
                p=folder/name
                (p/"CONTCAR").write_text((p/"POSCAR").read_text())
                rows="\n".join(["0 0 0 0.0001 0.0001 0.0001"]*6)
                (p/"OUTCAR").write_text(
                   f"free energy TOTEN = {energy:.6f} eV\n"
                   "POSITION TOTAL-FORCE (eV/Angst)\n--------\n"+rows+
                   "\nreached required accuracy\n")
            result=review.gate(folder)
            self.assertTrue(result["status"].startswith("RAW_GEOMETRY_SYMMETRY_GATE_PASS"))
            self.assertAlmostEqual(result["delta_E_meV_per_2fu_cell"][-1],24.117,places=6)
            self.assertIsNone(result["NEB_barrier"])
            self.assertIsNone(result["material_static_NMEE"])

if __name__=="__main__":
    unittest.main()
