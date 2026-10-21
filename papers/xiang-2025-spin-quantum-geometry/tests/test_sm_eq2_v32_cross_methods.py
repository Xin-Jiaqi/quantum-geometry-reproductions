"""V32 cross check: V30 orbital Liouville and immutable V14 scalar snapshots."""
import json,sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from methods.sm_eq2_v32 import Model,probe,eq2_four_sources,eta0_eq2_and_printed_eq6
from methods.finite_eta import second_order
class CrossV32(unittest.TestCase):
    def test_independent_v30_liouville_full_matrix(self):
        for n in (2,3,4,5):
            m=Model(n);d=m.data(.23)
            def H(x,y):
                q=m.data(x)
                return q['U']@np.diag(q['en'])@q['U'].conj().T
            eq2=d['U']@eq2_four_sources(m,.23,.02).sum(0)@d['U'].conj().T
            v30=second_order(H,.23,0.,np.eye(n),mu=.25,temperature=.45,
                             eta=.02,step=.005)
            self.assertLess(np.max(abs(eq2-v30['rho2'])),2e-6)
    def test_four_parts_hermitian(self):
        for n in (2,3,4,5):
            a=eq2_four_sources(Model(n),.23,.02)
            self.assertEqual(a.shape,(4,n,n))
            self.assertLess(np.max(abs(a.sum(0)-a.sum(0).conj().T)),2e-11)
    def test_frozen_v14_numbers(self):
        ref=json.loads((ROOT/'results/REFERENCE_V32.json').read_text())
        for case in ref['cases']:
            out=probe(case['n'])
            for key in ('extra_max','triple_max','delta_spin_trace',
                        'eq2_offdiag_spin_trace','eq6_offdiag_spin_trace'):
                self.assertAlmostEqual(out[key],case[key],delta=1e-8)
            self.assertLess(out['eq2_vs_orbital_max'],2e-8)
    def test_flat_f_zero(self):
        m=Model(4);m.slopes[:]=0.;m.curvature[:]=0.
        r=eta0_eq2_and_printed_eq6(m,.23)
        self.assertLess(np.max(abs(r['extra'])),1e-13)
    def test_source_mismatch_not_fitted(self):
        r=probe(4)
        self.assertGreater(r['extra_max'],.002)
        self.assertGreater(abs(r['delta_spin_trace']),6e-4)
