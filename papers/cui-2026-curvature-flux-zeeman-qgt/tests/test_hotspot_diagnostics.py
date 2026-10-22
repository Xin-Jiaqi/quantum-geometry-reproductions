"""Round20 microscopic occupation-sector / k-hotspot diagnostic acceptance."""
import unittest
import numpy as np
from flux_lab.hexagonal_bz import hexagonal_symmetrized_model,hexagonal_bz_response
from flux_lab.primitive_bz import primitive_bz_response
from flux_lab.group_microscopic import build_group
from flux_lab.nonabelian_response import cluster_interband_kubo
from flux_lab.hotspot_diagnostics import pair_resolved,weighted_k_audit

K=np.array([.5145,-1.0456]);Z=.071+.004j
def point(group="6",k=K):
    H,V,S=hexagonal_symmetrized_model(k,group)
    return H,V,S

class HotspotOccupations(unittest.TestCase):
    def test_pointwise_pair_decomposition_reproduces_kubo_both_complex_branches(self):
        for group in ("6","6m'm'"):
            H,V,S=point(group)
            d=pair_resolved(H,V,S,z=Z)
            plus=cluster_interband_kubo(H,V,S,z=Z,mu=.13,temperature=.17)[0]
            minus=cluster_interband_kubo(H,V,S,z=-Z,mu=.13,temperature=.17)[0]
            np.testing.assert_allclose(d["plus"],plus,atol=5e-13)
            np.testing.assert_allclose(d["minus"],minus,atol=5e-13)
            np.testing.assert_allclose(d["even"],(plus+minus)/2,atol=5e-13)
            np.testing.assert_allclose(d["odd"],(plus-minus)/2,atol=5e-13)

    def test_cross_mu_and_same_side_partition_sums_to_total(self):
        H,V,S=point()
        d=pair_resolved(H,V,S)
        np.testing.assert_allclose(d["components"]["cross_mu"]+
                                   d["components"]["same_side"],d["even"],atol=5e-13)

    def test_near_pole_small_occupation_difference_not_assumed_large_response(self):
        # At these k points a tiny spin splitting can nearly match omega,
        # but same-side Fermi occupations differ by only ~1e-6.
        H,V,S=point(k=np.array([-.77308044,2.03298929]))
        d=pair_resolved(H,V,S)
        self.assertGreater(d["resonant_transition_count"],0)
        self.assertLess(d["max_occupation_difference_resonant"],2e-5)
        self.assertLess(d["min_pole_distance"],.0045)

    def test_hotspot_group6_has_severe_phase_cancellation(self):
        d=weighted_k_audit("6",order=6,method="hexagon",
                           components=((0,0),),top_n=4)
        comp=d["components"]["00"]
        self.assertGreater(comp["cancellation_ratio"],100.)
        self.assertGreater(comp["absolute_integral"],.002)
        self.assertLess(comp["net_abs"],2e-5)
        self.assertEqual(len(comp["top_nodes"]),4)

    def test_sectors_and_integrated_response_at_identical_hex_grid(self):
        d=weighted_k_audit("6",order=4,method="hexagon",components=((0,0),(0,1)))
        ref=hexagonal_bz_response("6",order=4)
        for key in ("00","01"):
            a=int(key[0]);b=int(key[1])
            want=ref["even"][a,b];got=complex(*d["components"][key]["net"])
            self.assertAlmostEqual(abs(got-want),0,delta=3e-12)
        self.assertLess(d["max_total_partition_residual"],1e-12)

    def test_unprojected_primitive_response_matches_partition_at_fixed_grid(self):
        d=weighted_k_audit("6",order=5,method="gauss",components=((0,0),))
        ref=primitive_bz_response("6",order=5,method="gauss")
        val=complex(*d["components"]["00"]["net"])
        self.assertLess(abs(val-ref["even"][0,0]),3e-12)

    def test_thermal_occupation_sensitivity_is_not_grid_gap_disappearance(self):
        # Same H(k), only temperature changes the thermal factors.
        H,V,S=point()
        a=pair_resolved(H,V,S,temperature=.08)
        b=pair_resolved(H,V,S,temperature=.25)
        self.assertAlmostEqual(a["min_intercluster_gap"],b["min_intercluster_gap"],places=12)
        self.assertGreater(np.max(abs(a["components"]["same_side"]-
                                      b["components"]["same_side"])),1e-6)

    def test_antiunitary_mirror_forbidden_xx_and_allowed_xy(self):
        d=weighted_k_audit("6m'm'",order=4,method="hexagon",
                           components=((0,0),(0,1)))
        self.assertLess(d["components"]["00"]["net_abs"],1e-10)
        self.assertGreater(d["components"]["01"]["net_abs"],1e-8)

    def test_invalid_models_parameters_and_components_fail_closed(self):
        H,V,S=point()
        with self.assertRaises(ValueError):pair_resolved(H,V,S,temperature=0)
        with self.assertRaises(ValueError):pair_resolved(H,V,S,z=.03-.01j)
        with self.assertRaises(ValueError):pair_resolved(H,V,S,resonance_widths=0)
        with self.assertRaises(ValueError):weighted_k_audit("6",method="bogus")
        with self.assertRaises(ValueError):weighted_k_audit("6",components=((2,1),),order=2)

    def test_diagnostics_are_not_raw_probability_or_convergence_tests(self):
        d=weighted_k_audit("6",order=3,method="midpoint",components=((0,0),))
        self.assertIn("diagnostic only",d["status"])
        self.assertAlmostEqual(d["weight_sum"],1,places=12)
        self.assertGreaterEqual(d["components"]["00"]["cancellation_ratio"],1.)
        self.assertGreaterEqual(d["resonant_transition_visits"],0)

if __name__=="__main__":unittest.main()
