"""Xiang Fig1(a) Eq13 V29 matrix and frozen numeric gates."""
import json
import sys
import unittest
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P))
import paper_models as model
from paper_eq13 import intrinsic_reduced_eq13,kernel_matrix_at_k,kernel_analytic_at_k

VX=.6582119569*100.; ALPHA=.1*VX

class IntrinsicEq13PaperTests(unittest.TestCase):
    def test_analytic_vs_matrix_geometry_all_valleys_bands(self):
        for s in (-1,1):
            for b in (-1,1):
                for k in ((.12,-.07),(-.31,.26)):
                    a=kernel_analytic_at_k(*k,s,b,VX,VX,ALPHA,10.)
                    m=kernel_matrix_at_k(*k,s,b,VX,VX,ALPHA,10.)
                    for key in a:
                        self.assertAlmostEqual(a[key],m[key],delta=1e-12)
    def test_fermi_sea_surface_both_frozen(self):
        frozen=json.loads((P/'results/REFERENCE_V29.json').read_text())
        for p in frozen['intrinsic_eq13']:
            got=intrinsic_reduced_eq13(p['mu_meV'],VX,VX,ALPHA,10.,nr=48,ntheta=96)
            for key in ('sea','surface','total'):
                self.assertAlmostEqual(got[key],p[key],delta=1e-12)
    def test_mu20_numeric_exact_v27_reference(self):
        v=intrinsic_reduced_eq13(20.,VX,VX,ALPHA,10.,nr=48,ntheta=96)
        self.assertAlmostEqual(v['total'],1.6079096721715977e-05,delta=1e-13)
        self.assertAlmostEqual(v['sea']+v['surface'],v['total'],delta=1e-18)
    def test_valence_hole_and_electron_symmetry(self):
        for mu in (12.,20.,40.):
            pos=intrinsic_reduced_eq13(mu,VX,VX,ALPHA,10.,nr=32,ntheta=64)
            neg=intrinsic_reduced_eq13(-mu,VX,VX,ALPHA,10.,nr=32,ntheta=64)
            for key in pos:self.assertAlmostEqual(pos[key],-neg[key],delta=1e-12)
    def test_inside_indirect_gap_exact_zero(self):
        for mu in (0.,9.9,-9.9):
            self.assertEqual(intrinsic_reduced_eq13(mu,VX,VX,ALPHA,10.)['total'],0.)
    def test_quadrature_convergence(self):
        for mu in (12.,20.,40.):
            coarse=intrinsic_reduced_eq13(mu,VX,VX,ALPHA,10.,nr=32,ntheta=64)
            fine=intrinsic_reduced_eq13(mu,VX,VX,ALPHA,10.,nr=64,ntheta=128)
            for key in ('sea','surface','total'):
                self.assertAlmostEqual(coarse[key],fine[key],delta=2e-10)
    def test_no_type_ii_cone_allowed(self):
        with self.assertRaises(ValueError):
            intrinsic_reduced_eq13(20.,VX,VX,1.2*VX,10.)
    def test_reproduce_cli_explicit_protocol(self):
        d=model.reproduce()
        self.assertEqual(d['scientific_gates']['fig1a_intrinsic_eq13'],'PASS_REDUCED_FORMULA_ONLY')
        self.assertEqual(d['scientific_gates']['metal_DC_protocol'],'NOT_VALIDATED')
        self.assertEqual(len(d['fig1a']['intrinsic_eq13_reduced_mu_scan']),4)

if __name__=='__main__':unittest.main()
