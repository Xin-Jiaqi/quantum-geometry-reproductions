"""v53 independent original SM orbital Eq48/49/52/56 physical-model gates."""
import unittest
import numpy as np
from qgmap.orbital_sm_v53 import (metric_eq49_to_eq52,orbital_moment_three_band,
                            orbital_phe_eq56,full_audit,SM_PDF_SHA256)

class OriginalOrbitalSMV53(unittest.TestCase):
    def test_metric_source48_disagrees_with_printed_eq52(self):
        d=metric_eq49_to_eq52(48)
        self.assertGreater(abs(d['residual_original']),6e-6)
        self.assertGreater(abs(d['missing_velocity']),6e-6)

    def test_correct_coefficient_two_converges_to_eq48(self):
        ds=[metric_eq49_to_eq52(n) for n in (20,28,36,48)]
        errs=[abs(d['residual_correct']) for d in ds]
        self.assertTrue(all(a>b for a,b in zip(errs,errs[1:])))
        self.assertLess(errs[-1],1.5e-7)
        self.assertGreater(abs(ds[-1]['residual_original']),40*errs[-1])

    def test_missing_velocity_is_mesh_invariant_physical_bloch_metric(self):
        ds=[metric_eq49_to_eq52(n) for n in (20,28,48)]
        self.assertLess(max(x['missing_velocity'] for x in ds)-min(x['missing_velocity'] for x in ds),1e-10)
        self.assertAlmostEqual(ds[-1]['missing_velocity'],6.9295435748e-6,delta=1e-9)

    def test_when_band_sum_velocity_zero_original_eq52_metric_agrees(self):
        for n in (20,32,48):
            a=metric_eq49_to_eq52(n,tilt_scale=0.)
            self.assertAlmostEqual(a['missing_velocity'],0.,delta=1e-15)
            self.assertLess(abs(a['residual_original']),1e-16)

    def test_3band_crossed_orbital_moment_Berry_process_large(self):
        d=orbital_moment_three_band(56,mu=-1.1,T=.45)
        self.assertAlmostEqual(d['omitted_threeband'],.00074273209625,delta=1e-8)
        self.assertGreater(d['ratio_abs_to_complete'],.9)
        self.assertLess(d['ratio_abs_to_complete'],1.)
        self.assertLess(d['complete_moment_Berry'],0.)

    def test_cross_term_is_not_total_orbital_MNHE_fraction(self):
        d=orbital_moment_three_band(48,mu=-1.1,T=.45)
        self.assertAlmostEqual(d['omitted_threeband'],d['complete_moment_Berry']-d['printed_eq52_m_eq_l'],delta=1e-15)
        self.assertAlmostEqual(sum(row[2] for row in d['per_band']),d['omitted_threeband'],delta=1e-15)

    def test_threeband_near_exact_brillouin_mesh_convergence(self):
        a=orbital_moment_three_band(32,mu=-1.1,T=.45)
        b=orbital_moment_three_band(56,mu=-1.1,T=.45)
        self.assertLess(abs(a['omitted_threeband']-b['omitted_threeband']),1e-10)

    def test_threeband_interband_U1_gauge_invariance(self):
        a=orbital_moment_three_band(36,mu=-1.1,T=.45)
        b=orbital_moment_three_band(36,mu=-1.1,T=.45,band_phases=(.38,-1.1,2.4))
        for key in ('complete_moment_Berry','printed_eq52_m_eq_l','omitted_threeband'):
            self.assertAlmostEqual(a[key],b[key],delta=1e-14)

    def test_cross_term_changes_fermi_window_and_need_not_dominate(self):
        a=orbital_moment_three_band(36,mu=-1.1,T=.45)
        b=orbital_moment_three_band(36,mu=1.5,T=.45)
        self.assertGreater(a['ratio_abs_to_complete'],.9)
        self.assertLess(b['ratio_abs_to_complete'],.07)

    def test_eq56_both_printed_integral_representations_agree_in_limit(self):
        seq=[orbital_phe_eq56(N) for N in (24,32,48,64)]
        delta=[abs(v['absolute_delta']) for v in seq]
        self.assertTrue(all(a>b for a,b in zip(delta,delta[1:])))
        self.assertLess(delta[-1],3e-6)
        self.assertGreater(seq[-1]['first_eq56'],.007)

    def test_eq56_finite_difference_error_approaches_second_order(self):
        a=orbital_phe_eq56(32);b=orbital_phe_eq56(64)
        factor=abs(a['absolute_delta']/b['absolute_delta'])
        self.assertGreater(factor,3.6)
        self.assertLess(factor,4.3)

    def test_invalid_mesh_source_provenance_and_PARTIAL(self):
        for fun in (metric_eq49_to_eq52,orbital_moment_three_band,orbital_phe_eq56):
            with self.assertRaises(ValueError):fun(8)
        with self.assertRaises(ValueError):orbital_moment_three_band(32,band_phases=(0.,1.))
        self.assertEqual(len(SM_PDF_SHA256),64)
        audit=full_audit()
        self.assertEqual(audit['scientific_status'],'PARTIAL')
        self.assertIn('Not a total-MNHE percentage',audit['scope']['threeband'])

if __name__=='__main__':unittest.main()
