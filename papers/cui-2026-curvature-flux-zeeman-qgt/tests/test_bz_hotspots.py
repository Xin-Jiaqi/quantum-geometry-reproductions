"""Round20 spectral hotspot, occupation and cancellation physics checks."""
import unittest
import numpy as np
from flux_lab.hexagonal_bz import hexagonal_symmetrized_model,hexagonal_bz_response
from flux_lab.primitive_bz import primitive_bz_response
from flux_lab.nonabelian_response import cluster_interband_kubo
from flux_lab.bz_hotspots import (local_transition_diagnostic,field_diagnostics,
                                   bz_hotspot_audit)

K=np.array([.43,-.39])

class BZHotspots(unittest.TestCase):
    def test_pair_denominator_oracle_vs_original_causal_interband_kubo(self):
        for group in ("6","6m'm'"):
            H,V,S=hexagonal_symmetrized_model(K,group)
            for z in (.071+.004j,.051+.019j,-.043+.017j):
                d=local_transition_diagnostic(H,V,S,z=z)
                p=cluster_interband_kubo(H,V,S,z=z,mu=.13,temperature=.17)[0]
                m=cluster_interband_kubo(H,V,S,z=-z,mu=.13,temperature=.17)[0]
                np.testing.assert_allclose(d["even"],(p+m)/2,atol=3e-13,rtol=1e-12)
                np.testing.assert_allclose(d["odd"],(p-m)/2,atol=3e-13,rtol=1e-12)

    def test_occupation_reporting_threshold_does_not_truncate_response(self):
        H,V,S=hexagonal_symmetrized_model(K,"6")
        a=local_transition_diagnostic(H,V,S,occupation_cutoff=0.)
        b=local_transition_diagnostic(H,V,S,occupation_cutoff=2.)
        for name in ("even","odd"):
            np.testing.assert_allclose(a[name],b[name],atol=1e-14,rtol=1e-13)
        self.assertGreater(a["active_transition_count"],0)
        self.assertEqual(b["active_transition_count"],0)
        self.assertIsNone(b["closest_active_resonance_distance"])

    def test_energy_origin_plus_mu_shift_is_invariant(self):
        H,V,S=hexagonal_symmetrized_model(K,"6")
        d=local_transition_diagnostic(H,V,S,mu=.13)
        e=local_transition_diagnostic(H+4.7*np.eye(4),V,S,mu=4.83)
        for name in ("even","odd"):
            np.testing.assert_allclose(d[name],e[name],atol=4e-13,rtol=1e-12)
        self.assertAlmostEqual(d["min_intercluster_gap"],e["min_intercluster_gap"],places=11)

    def test_positive_intercluster_gap_and_nontrivial_occupation(self):
        H,V,S=hexagonal_symmetrized_model(K,"6")
        d=local_transition_diagnostic(H,V,S)
        self.assertGreater(d["min_intercluster_gap"],0.)
        self.assertGreater(d["active_transition_count"],0)
        self.assertGreater(d["max_pair_occupation_difference"],.1)
        self.assertGreaterEqual(d["closest_active_resonance_distance"],0.)

    def test_synthetic_positive_integrand_has_no_cancellation(self):
        field=np.ones((6,2,3),complex)*(.2+.7j)
        d=field_diagnostics(field,np.ones(6)/6,component=(0,1))
        self.assertAlmostEqual(d["selected_retention"],1.,places=13)
        self.assertAlmostEqual(d["selected_cancellation"],0.,places=13)
        self.assertAlmostEqual(d["top_share"],1/6,places=13)

    def test_synthetic_opposite_complex_hotspots_cancel(self):
        field=np.zeros((2,2,3),complex)
        field[0,1,2]=1+2j;field[1,1,2]=-1-2j
        d=field_diagnostics(field,[.5,.5],component=(1,2),top_fraction=.5)
        self.assertAlmostEqual(abs(d["selected_net"]),0.,places=14)
        self.assertAlmostEqual(d["selected_cancellation"],1.,places=14)
        self.assertAlmostEqual(d["top_share"],.5,places=14)

    def test_weighted_abs_triangle_inequality_all_tensor_components(self):
        H,V,S=hexagonal_symmetrized_model(K,"6")
        f=np.stack([local_transition_diagnostic(H,V,S)["even"]]*7)
        d=field_diagnostics(f,np.ones(7)/7)
        self.assertTrue(np.all(abs(d["net"])<=d["integrated_abs"]+1e-13))
        self.assertTrue(0<=d["selected_retention"]<=1+1e-13)

    def test_field_diagnostics_fail_closed(self):
        X=np.zeros((4,2,3),complex)
        with self.assertRaises(ValueError):field_diagnostics(X,np.ones(3)/3)
        with self.assertRaises(ValueError):field_diagnostics(X,[1,1,1,1])
        with self.assertRaises(ValueError):field_diagnostics(X,[.25]*4,component=(2,0))
        with self.assertRaises(ValueError):field_diagnostics(X,[.25]*4,top_fraction=0.)

    def test_hexagon_diagnostic_unprojected_response_equals_original(self):
        d=bz_hotspot_audit("6",domain="hexagon",order=3)
        ref=hexagonal_bz_response("6",order=3)
        np.testing.assert_allclose(d["even"]["net"],ref["even"],atol=3e-13,rtol=1e-11)
        np.testing.assert_allclose(d["odd"]["net"],ref["odd"],atol=3e-13,rtol=1e-11)

    def test_primitive_diagnostic_unprojected_response_equals_original(self):
        d=bz_hotspot_audit("6",domain="primitive",method="midpoint",
                            order=5,shift=(.23,-.17))
        ref=primitive_bz_response("6",order=5,method="midpoint",shift=(.23,-.17))
        np.testing.assert_allclose(d["even"]["net"],ref["even"],atol=3e-13,rtol=1e-11)
        np.testing.assert_allclose(d["odd"]["net"],ref["odd"],atol=3e-13,rtol=1e-11)

    def test_hotspot_concentration_weighted_not_unweighted(self):
        d=bz_hotspot_audit("6",domain="hexagon",order=3,hotspot_count=4)
        self.assertEqual(len(d["hotspots"]),4)
        self.assertGreaterEqual(d["even"]["top_share"],0.)
        self.assertLessEqual(d["even"]["top_share"],1.+1e-12)
        self.assertTrue(0<=d["even"]["selected_cancellation"]<=1.+1e-12)
        self.assertGreater(d["min_gap"],0.)
        self.assertTrue(np.all(np.isfinite(d["even"]["net"])))

    def test_fixed_model_temperature_changes_spectral_occupation(self):
        H,V,S=hexagonal_symmetrized_model(K,"6")
        a=local_transition_diagnostic(H,V,S,temperature=.1)
        b=local_transition_diagnostic(H,V,S,temperature=.5)
        self.assertGreater(np.max(abs(a["even"]-b["even"])),1e-12)

    def test_invalid_domain_and_hotspot_count(self):
        with self.assertRaises(ValueError):bz_hotspot_audit("6",domain="wrong",order=3)
        with self.assertRaises(ValueError):bz_hotspot_audit("6",domain="hexagon",order=3,hotspot_count=0)
        with self.assertRaises(ValueError):bz_hotspot_audit("4mm",domain="hexagon",order=3)

if __name__=="__main__":unittest.main()
