"""Round21: complement to accepted Round20 hotspot module — Fermi sectors."""
import unittest
import numpy as np
from flux_lab.occupation_partition import pair_resolved,integrate_fermi_sectors
from flux_lab.hexagonal_bz import hexagonal_symmetrized_model,hexagonal_bz_response
from flux_lab.primitive_bz import primitive_bz_response
from flux_lab.nonabelian_response import cluster_interband_kubo
from flux_lab.bz_hotspots import bz_hotspot_audit

K=np.array([.5145,-1.0456]);Z=.071+.004j
class FermiSectors(unittest.TestCase):
    def test_direct_pair_sum_equals_existing_full_complex_kubo(self):
        for group in ("6","6m'm'"):
            H,V,S=hexagonal_symmetrized_model(K,group)
            d=pair_resolved(H,V,S)
            plus=cluster_interband_kubo(H,V,S,z=Z,mu=.13,temperature=.17)[0]
            minus=cluster_interband_kubo(H,V,S,z=-Z,mu=.13,temperature=.17)[0]
            np.testing.assert_allclose(d["plus"],plus,atol=3e-13)
            np.testing.assert_allclose(d["minus"],minus,atol=3e-13)
            np.testing.assert_allclose(d["even"],(plus+minus)/2,atol=3e-13)

    def test_cross_mu_same_side_reconstruct_exact_local_kubo(self):
        H,V,S=hexagonal_symmetrized_model(K,"6")
        d=pair_resolved(H,V,S)
        np.testing.assert_allclose(d["components"]["cross_mu"]+
                                   d["components"]["same_side"],d["even"],atol=3e-13)

    def test_near_resonance_denominator_is_not_occupancy_weight(self):
        H,V,S=hexagonal_symmetrized_model([-.77308044,2.03298929],"6")
        d=pair_resolved(H,V,S)
        self.assertGreater(d["resonant_transition_count"],0)
        self.assertLess(d["max_occupation_difference_resonant"],2e-5)
        self.assertLess(d["min_pole_distance"],.0045)

    def test_fixed_H_thermal_scan_changes_same_side_not_gap(self):
        H,V,S=hexagonal_symmetrized_model(K,"6")
        a=pair_resolved(H,V,S,temperature=.08)
        b=pair_resolved(H,V,S,temperature=.25)
        self.assertAlmostEqual(a["min_intercluster_gap"],
                               b["min_intercluster_gap"],places=12)
        self.assertGreater(np.max(abs(a["components"]["same_side"]-
                                      b["components"]["same_side"])),1e-6)

    def test_hex_finite_grid_sector_sum_matches_independent_hotspot_evaluator(self):
        d=integrate_fermi_sectors("6",order=3,component=(0,0))
        ref=bz_hotspot_audit("6",domain="hexagon",order=3,component=(0,0))
        val=complex(*d["sectors"]["total"]["net"])
        self.assertLess(abs(val-ref["even"]["net"][0,0]),3e-12)
        self.assertLess(d["complete_partition_residual"],3e-13)

    def test_primitive_sector_sum_matches_previous_primitive_BZ_solver(self):
        d=integrate_fermi_sectors("6",order=5,domain="primitive",
                                 method="gauss",component=(0,1))
        ref=primitive_bz_response("6",order=5,method="gauss")
        self.assertLess(abs(complex(*d["sectors"]["total"]["net"])-ref["even"][0,1]),3e-12)

    def test_high_cancellation_and_opposite_group_behaviour(self):
        a=integrate_fermi_sectors("6",order=6,component=(0,0))
        b=integrate_fermi_sectors("6m'm'",order=4,component=(0,1))
        self.assertGreater(a["cancellation_ratio"],100)
        self.assertLess(b["cancellation_ratio"],10)
        self.assertGreater(a["sectors"]["cross_mu"]["weighted_absolute"],
                           100*a["sectors"]["same_side"]["weighted_absolute"])
        self.assertGreater(b["sectors"]["same_side"]["net_abs"],
                           10*b["sectors"]["cross_mu"]["net_abs"])

    def test_bad_inputs_fail_closed(self):
        H,V,S=hexagonal_symmetrized_model(K,"6")
        for kwargs in ({"temperature":0},{"z":.01-.1j},{"resonance_widths":0}):
            with self.assertRaises(ValueError):pair_resolved(H,V,S,**kwargs)
        with self.assertRaises(ValueError):integrate_fermi_sectors("6",component=(2,1))
        with self.assertRaises(ValueError):integrate_fermi_sectors("6",domain="other")
if __name__=="__main__":unittest.main()
