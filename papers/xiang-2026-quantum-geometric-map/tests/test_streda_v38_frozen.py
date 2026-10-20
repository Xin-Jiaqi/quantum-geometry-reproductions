"""Frozen numerical reference for v38 independent clean antisymmetric sea benchmark."""
import json
from pathlib import Path
import subprocess
import sys
import unittest
from qgmap.streda_clean_sea_v38 import bulk_sea,metallic_pole_audit

ROOT=Path(__file__).resolve().parents[1]
FILE=ROOT/'results/frozen/v38_clean_streda_reference.json'

class CleanSeaFrozen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref=json.loads(FILE.read_text())
        cls.to=cls.ref['topology']; cls.es=cls.ref['eta_scan']
    def test_insulator_phase_reference(self):
        for m,d in self.to.items():
            with self.subTest(m=m):
                r=bulk_sea(61,float(m),eta=.02)
                self.assertAlmostEqual(r['sea_sigma_over_e2h'],d['sea_eta_0p02'],delta=1e-11)
                self.assertAlmostEqual(r['occupied_berry_chern'],d['occupied_berry'],delta=1e-11)
    def test_eta_scan(self):
        for eta,d in self.es.items():
            with self.subTest(eta=eta):
                r=bulk_sea(61,-1,eta=float(eta))
                self.assertAlmostEqual(r['sea_sigma_over_e2h'],d['sea'],delta=1e-11)
    def test_metal_pole_cancellation(self):
        x=metallic_pole_audit(101,-1,1.1,.02)
        for k,v in self.ref['metal'].items():
            if isinstance(v,(float,int)):
                self.assertAlmostEqual(x[k],v,delta=2e-11)
    def test_cli_science_contract(self):
        command=[sys.executable,str(ROOT/'scripts/reproduce_streda_v38.py'),'--mesh','61','--mass','-1','--eta','.02']
        p=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,timeout=30)
        self.assertEqual(p.returncode,0,msg=p.stderr)
        a=json.loads(p.stdout)
        self.assertEqual(a['scientific_status'],'PARTIAL')
        self.assertEqual(a['scope'],'clean two-band antisymmetric Kubo-Streda Fermi-sea benchmark ONLY')
        self.assertAlmostEqual(a['fhs_chern'],-1.,delta=1e-12)
if __name__=='__main__':unittest.main()
