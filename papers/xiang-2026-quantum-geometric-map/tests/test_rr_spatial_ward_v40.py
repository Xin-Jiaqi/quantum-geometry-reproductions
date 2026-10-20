"""Independent real-vector-potential finite differences: test spatial RR Ward."""
import json
import subprocess
import sys
import unittest
import numpy as np
from pathlib import Path
from qgmap.rr_spatial_ward_v40 import audit,twisted_qwz

class SpatialWard40(unittest.TestCase):
    def test_both_axes_independent_twist(self):
        v=audit(mesh=21,W=.12,energy=1.1,step=2e-4)
        self.assertEqual(v['science_status'],'PARTIAL')
        for a in ('x','y'):
            with self.subTest(axis=a):
                d=v['directions'][a]
                self.assertLess(d['rr_ladder_residual'],1e-10)
                self.assertLess(d['ra_ladder_residual'],1e-10)
                self.assertLess(d['twisted_hamiltonian_velocity_error'],2e-7)
                self.assertLess(d['dSigma_twist_vs_RR_ladder_error'],8e-5)
                self.assertLess(d['dGR_twist_vs_RR_vertex_error'],2e-3)
                self.assertGreater(d['RA_current_correction_norm'],5*d['RR_current_correction_norm'])
                self.assertLess(d['selfenergy_residual_plus'],1e-9)

    def test_zero_disorder_twist(self):
        v=audit(mesh=17,W=0.,energy=.7,step=3e-4)
        for d in v['directions'].values():
            self.assertEqual(d['RR_current_correction_norm'],0.)
            self.assertEqual(d['RA_current_correction_norm'],0.)
            self.assertEqual(d['dSigma_twist_vs_RR_ladder_error'],0.)
            self.assertLess(d['dGR_twist_vs_RR_vertex_error'],5e-3)

    def test_translated_model_vs_direct_bloch(self):
        n=17;dx=.034;dy=-.027
        h=twisted_qwz(n,-1.,dx,dy)[0]
        shifted=twisted_qwz(n,-1.,offset=(n*dx/(2*np.pi),n*dy/(2*np.pi)))[0]
        np.testing.assert_allclose(h,shifted,atol=1e-14)

    def test_fails_closed_and_independent_cli(self):
        with self.assertRaises(ValueError):audit(step=0)
        with self.assertRaises(ValueError):audit(step=np.inf)
        with self.assertRaises(ValueError):twisted_qwz(21,-1.,ax=np.nan)
        root=Path(__file__).resolve().parents[1]
        run=subprocess.run([sys.executable,str(root/'scripts/reproduce_rr_ward_v40.py'),'--mesh','17','--W','0'],
                           cwd=root,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=45)
        self.assertEqual(run.returncode,0,msg=run.stderr)
        data=json.loads(run.stdout)
        self.assertEqual(data['science_status'],'PARTIAL')
        self.assertIn('x',data['directions'])
        self.assertIn('y',data['directions'])

if __name__=='__main__':unittest.main()
