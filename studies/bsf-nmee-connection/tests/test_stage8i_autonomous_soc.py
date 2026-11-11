"""Stage8i: offline-only source-anchored SOC job generator, parser and mock PBS tests.

Never runs VASP or submits real qsub. All sources are synthetic fixtures or
original GitHub Stage8g input geometry. Real Stage8h CONTCAR SHA assertions
are separately checked against the committed sanitized raw audit.
"""
import importlib.util,json,re,subprocess,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
STAGE8I=ROOT/"hpc-stage8i"
sys.path.insert(0,str(STAGE8I))

def load(name):
    target=STAGE8I/(name+".py")
    spec=importlib.util.spec_from_file_location("stage8i_test_"+name,target)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

prep=load("prepare_soc")
pbs=load("onsite_soc")
eig=load("eigenval_gate")

class Stage8iAutonomousSOC(unittest.TestCase):
    def _fixture(self,tmp):
        base=Path(tmp)/"stage8g"
        # Reference starting geometry as synthetic fixture, no DFT claims.
        source=ROOT/"material-screening/source-structures/WSe2_BiDB_row2685_3R_polar_zscan.vasp"
        for old in ("path_A_s0.00","path_A_s0.50"):
            folder=base/old;folder.mkdir(parents=True)
            (folder/"CONTCAR").write_bytes(source.read_bytes())
        return base

    def test_anchored_source_hashes_are_from_stage8h_raw_evidence(self):
        audit=json.loads((ROOT/"hpc-stage8g/STAGE8H_AUDIT_RAW_20261021.json").read_text())
        for name,j in prep.EXPECTED.items():
            match=next(x for x in audit["jobs"] if x["job"]==j["stage8g"])
            self.assertEqual(j["sha256"],match["raw_output_checks"]["CONTCAR_sha256"])
            self.assertTrue(match["symmetry_pass"])
        self.assertEqual(set(prep.EXPECTED),set(pbs.NAMES))

    def test_exactly_two_jobs_and_no_potcar_or_job_submission(self):
        with tempfile.TemporaryDirectory() as t:
            base=self._fixture(t)
            # synthetic fixtures have synthetic hashes, never accept as actual DFT.
            patched={k:{**v,"sha256":prep.digest(base/v["stage8g"]/"CONTCAR")}
                     for k,v in prep.EXPECTED.items()}
            with patch.object(prep,"EXPECTED",patched):
                m=prep.build(base,Path(t)/"prepared")
            self.assertEqual(len(m["jobs"]),2)
            self.assertFalse(list((Path(t)/"prepared").rglob("POTCAR")))
            _,dat=pbs.verify_inputs(Path(t)/"prepared")
            self.assertEqual(dat["max_new_jobs"],2)
            for key in pbs.NAMES:
                text=(Path(t)/"prepared"/key/"INCAR").read_text()
                self.assertIn("LSORBIT = .TRUE.",text)
                self.assertIn("NSW = 0",text)
                self.assertNotRegex(text,r"(?m)^\s*ISPIN\s*=")
            self.assertFalse(any(x.endswith("WAVECAR") for x in
                (y.name for y in (Path(t)/"prepared").rglob("*"))))
            with self.assertRaisesRegex(RuntimeError,"Refuse overwriting"):
                with patch.object(prep,"EXPECTED",patched):
                    prep.build(base,Path(t)/"prepared")

    def test_tampering_stage8h_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            base=self._fixture(t)
            with self.assertRaisesRegex(RuntimeError,"FAILED_STAGE8H_CONTCAR_SHA"):
                prep.build(base,Path(t)/"prepared")

    def test_two_jobs_mock_pbs_with_no_duplicate_submissions(self):
        with tempfile.TemporaryDirectory() as t:
            base=self._fixture(t)
            patched={k:{**v,"sha256":prep.digest(base/v["stage8g"]/"CONTCAR")}
                     for k,v in prep.EXPECTED.items()}
            with patch.object(prep,"EXPECTED",patched):
                root=Path(t)/"prepared";prep.build(base,root)
            # Patch both in-memory expected potcar checksum and VALIDATED site.
            vasp=Path(t)/"vasp_ncl";vasp.write_text("FAKE_BINARY_NEVER_EXECUTED")
            site=Path(t)/"site.json";site.write_text("{}")
            cfg={"scheduler":"pbs","validated":True,
                 "license_confirmed":True,"soc_binary_verified":True,
                 "queue":"CP_001","ncpus":16,"mem_gb":32,
                 "walltime":"12:00:00","vasp_binary":str(vasp),
                 "mpi_command":"mpirun"}
            pot=("TITEL = PAW_PBE W 08Apr2002\nENMAX = 223.057;\n"
                 "TITEL = PAW_PBE Se 06Sep2000\nENMAX = 211.555;\n")
            for name in pbs.NAMES:(root/name/"POTCAR").write_text(pot)
            ids=iter(["500001.cmspbs1","500002.cmspbs1"])
            def fake_qsub(*a,**kw):
                return subprocess.CompletedProcess(a[0],0,
                    stdout=next(ids)+"\n",stderr="")
            with patch.object(pbs,"site_gate",return_value=cfg), \
                 patch.object(pbs,"EXPECTED_POTCAR_SHA",pbs.digest(root/pbs.NAMES[0]/"POTCAR")), \
                 patch.object(pbs.subprocess,"run",side_effect=fake_qsub) as qsub:
                dry=pbs.launch(root,site,submit=False)
                self.assertEqual(dry["status"],"VALIDATED_DRY_RUN")
                first=pbs.launch(root,site,submit=True)
                self.assertEqual(len(first["pbs_ids"]),2)
                second=pbs.launch(root,site,submit=True)
                self.assertEqual(second["pbs_ids"],first["pbs_ids"])
                self.assertEqual(qsub.call_count,2)
            self.assertEqual(pbs.inspect(root)["material_static_NMEE"],None)

    def test_soc_site_must_confirm_real_noncollinear_executable(self):
        with tempfile.TemporaryDirectory() as t:
            site=Path(t)/"site.json"
            site.write_text(json.dumps({"scheduler":"pbs","validated":True,
                "license_confirmed":True,"soc_binary_verified":False,
                "queue":"CP_001","ncpus":16,"mem_gb":32,
                "walltime":"12:00:00","vasp_binary":"/no/vasp_ncl",
                "mpi_command":"mpirun"}))
            with self.assertRaisesRegex(RuntimeError,"noncollinear LSORBIT"):
                pbs.site_gate(site)

    def test_noncollinear_eigenval_parser_is_conservative(self):
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/"EIGENVAL"
            header=["header"]*5+["2 2 16"]
            data=[]
            for k in range(2):
                data.append("")
                data.append(f"{.5*k:.8f} 0.000000 0.000000 0.5")
                for j in range(16):
                    energy=((-1.0+.1*k) if j==1 else
                            (.2+.1*k) if j==2 else j-3)
                    occupancy=1. if j<2 else 0.
                    data.append(f"{j+1} {energy:.6f} {occupancy:.1f}")
            f.write_text("\n".join(header+data)+"\n")
            result=eig.analyze_eigenval(f,expected_kpoints=2)
            self.assertEqual(result["num_kpoints"],2)
            self.assertAlmostEqual(result["sampled_indirect_gap_eV"],1.1,places=9)
            self.assertFalse(result["fractional_occupation_warning"])
            self.assertIn("sampling", "k-grid sampling")

if __name__=="__main__":
    unittest.main()
