import sys
import unittest
from pathlib import Path
from dataclasses import replace
import numpy as np
here=Path(__file__).resolve().parents[1]
if str(here) not in sys.path:sys.path.insert(0,str(here))
from scripts.reproduce import evaluate
from nmee_lab.ptbackend.cumnas import Model,points,extrinsic_mu_scan
from nmee_lab.ptbackend.split_response import split

class TestReproductionContract(unittest.TestCase):
    def test_branch_provenance(self):
        x=evaluate('fig2-printed-vs-static',19,.025,.15,10,50,1)
        vals=x['protocols']
        self.assertAlmostEqual(vals['static_eq9_full']-vals['prl_printed_eq1'],vals['static_eq9_minus_printed'],places=12)
        self.assertNotAlmostEqual(vals['static_eq9_full'],vals['prl_printed_eq1'],places=5)

    def test_finite_temperature_spin_space_projector_formula(self):
        c=evaluate('fig1d-printed',39,.03,.15,100,50,1.0)
        self.assertTrue(np.isfinite(c['response']))
        self.assertIn('PRL_Eq2',c['protocol'])

    def test_tau_scaling_eq2(self):
        a=evaluate('fig1d-printed',41,.03,.15,100,50,1)['response']
        b=evaluate('fig1d-printed',41,.03,.15,100,100,1)['response']
        self.assertAlmostEqual(b,2*a,places=11)

    def test_length_scaling_eq2(self):
        a=evaluate('fig1d-printed',37,.03,.15,100,50,1)['response']
        b=evaluate('fig1d-printed',37,.03,.15,100,50,2)['response']
        self.assertAlmostEqual(b,a/2,places=11)

    def test_neel_reversal_eq2(self):
        a=evaluate('fig1d-printed',47,.03,.15,100,50,1)['response']
        b=evaluate('fig1d-printed',47,.03,.15+np.pi,100,50,1)['response']
        self.assertLess(abs(a+b),2e-9)

    def test_staggered_intrinsic_static_diagnostic(self):
        a=evaluate('fig1c-static-diagnostic',37,.03,0,100,50,1)
        self.assertLess(abs(a['A_z']+a['B_z']),1e-8)
        self.assertIn('diagnostic',a['protocol'])

    def test_finite_freq_split_import_and_shape(self):
        # Regression: previous R32 prototype had an erroneous absolute import
        # of cumnas.KB on finite-T projector route.
        pts=points(5)
        out=split(pts,.15,.03,eta_e=1e-3)
        self.assertEqual(out['full'].shape,(25,))
        self.assertTrue(np.isfinite(out['full'].real).all())

if __name__=='__main__':unittest.main()
