"""Stage8g: local-only research pipeline tests. NEVER submit HPC jobs."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]/"hpc-stage8g"

def read_script(name):
    spec=importlib.util.spec_from_file_location("stage8g_"+name,ROOT/(name+".py"))
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

prepare=read_script("prepare")
onsite=read_script("onsite_pbs")

class Stage8GAutonomousHpcTests(unittest.TestCase):
    def test_deterministic_all_paths_and_three_job_pilot(self):
        with tempfile.TemporaryDirectory() as t:
            dest=Path(t)/"generated"
            m=prepare.build(dest)
            self.assertEqual(len(m["jobs"]),15)
            self.assertEqual([x["name"] for x in m["jobs"]
                              if x["selected_for_pilot"]],list(onsite.PILOT))
            self.assertAlmostEqual(m["shortest_slide_A"],3.319006624468/3**.5,places=8)
            self.assertEqual(len(list(dest.rglob("POSCAR"))),15)
            self.assertEqual(len(list(dest.rglob("INCAR"))),15)
            self.assertFalse(list(dest.rglob("POTCAR")))
            root,pilot=onsite.load_inputs(dest)
            self.assertEqual(len(pilot),3)
            self.assertTrue(all((root/x["name"]/"POSCAR").is_file() for x in pilot))
            self.assertTrue(all(y["status"].startswith("NOT_RUN") for y
                                in onsite.audit(dest)["rows"]))

    def test_rejects_duplicate_runs_and_input_tampering(self):
        with tempfile.TemporaryDirectory() as t:
            dest=Path(t)/"work"
            prepare.build(dest)
            with self.assertRaisesRegex(RuntimeError,"overwrite"):
                prepare.build(dest)
            f=dest/"path_A_s0.25"/"INCAR"
            f.write_text(f.read_text()+"\nENCUT=200\n")
            with self.assertRaisesRegex(RuntimeError,"Scientific input changed"):
                onsite.load_inputs(dest)

    def test_audit_requires_authentic_forces_and_frozen_lateral_coordinates(self):
        with tempfile.TemporaryDirectory() as t:
            dest=Path(t)/"pilot"
            prepare.build(dest)
            folder=dest/"path_A_s0.00"
            (folder/"CONTCAR").write_text((folder/"POSCAR").read_text())
            force="\n".join(["0.0 0.0 0.0 0.0001 0.0001 0.0001"]*6)
            (folder/"OUTCAR").write_text(
                "free energy TOTEN = -123.456789 eV\n"
                "POSITION TOTAL-FORCE (eV/Angst)\n----\n"+force+
                "\nreached required accuracy\n")
            row=onsite.audit(dest)["rows"][0]
            self.assertEqual(row["status"],"CONVERGED_REQUIRES_SCIENTIFIC_REVIEW")
            self.assertAlmostEqual(row["final_TOTEN_eV"],-123.456789)
            self.assertLess(row["max_lateral_xy_drift_A"],1e-8)
            (folder/"OUTCAR").write_text("free energy TOTEN = -123.456789 eV\n")
            self.assertEqual(onsite.audit(dest)["rows"][0]["status"],
                             "INVALID_GEOMETRY_OR_FORCES")

    def test_mock_pbs_submission_is_bounded_and_idempotent(self):
        from unittest.mock import patch
        import subprocess
        with tempfile.TemporaryDirectory() as t:
            dest=Path(t)/"pilot"
            prepare.build(dest)
            pot="TITEL = PAW_PBE W 06Sep2000\nENMAX = 340.0;\n"
            pot+="TITEL = PAW_PBE Se 06Sep2000\nENMAX = 280.0;\n"
            for name in onsite.PILOT:
                (dest/name/"POTCAR").write_text(pot)
            site=Path(t)/"local_site.json"
            site.write_text("{}")
            cfg={"scheduler":"pbs","queue":"CP_001","ncpus":4,
                 "mem_gb":8,"walltime":"02:00:00",
                 "vasp_binary":"/usr/bin/true","mpi_command":"mpirun",
                 "license_confirmed":True,"validated":True}
            seq=iter(("1001.test","1002.test","1003.test"))
            def fake_submit(*args,**kwargs):
                return subprocess.CompletedProcess(
                    args[0],0,stdout=next(seq)+"\n",stderr="")
            with patch.object(onsite,"load_site",return_value=cfg), \
                 patch.object(onsite.subprocess,"run",side_effect=fake_submit) as mock:
                a=onsite.dispatch(dest,site,execute=False)
                self.assertEqual(a["status"],"DRY_RUN_VALIDATED")
                b=onsite.dispatch(dest,site,execute=True)
                self.assertEqual(len(b["pbs_ids"]),3)
                again=onsite.dispatch(dest,site,execute=True)
                self.assertEqual(b["pbs_ids"],again["pbs_ids"])
                self.assertEqual(mock.call_count,3)
            self.assertEqual(len(json.loads((dest/"pbs_ids.json").read_text())),3)

    def test_source_sha_and_separate_material_provenance(self):
        lat,pos=prepare.source_atoms()
        self.assertEqual(len(pos),6)
        self.assertEqual(prepare.sha(prepare.SOURCE),prepare.EXPECTED_SHA)
        self.assertEqual(prepare.SOURCE_UID,"WSe2-1cfbe6183886-2-1_0_0_1--0.67_-0.33")

    def test_fail_closed_scheduler_and_pseudopotential(self):
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/"fake_potcar"
            f.write_text("TITEL = PAW_PBE Se\nENMAX= 400.0;\n"
                         "TITEL = PAW_PBE W\nENMAX= 400.0;\n")
            with self.assertRaises(RuntimeError):
                onsite.potcar_gate(f)
            f.write_text("TITEL = PAW_PBE W\nENMAX= 460.0;\n"
                         "TITEL = PAW_PBE Se\nENMAX= 350.0;\n")
            with self.assertRaisesRegex(RuntimeError,"ENCUT"):
                onsite.potcar_gate(f)
            config=Path(t)/"site.json"
            config.write_text(json.dumps({"scheduler":"pbs","queue":"CP_001",
                "ncpus":4,"mem_gb":8,"walltime":"02:00:00",
                "vasp_binary":"/not/a/real/vasp_std","mpi_command":"mpirun",
                "license_confirmed":False,"validated":False}))
            with self.assertRaisesRegex(RuntimeError,"Site not"):
                onsite.load_site(config)

if __name__=="__main__":
    unittest.main()
