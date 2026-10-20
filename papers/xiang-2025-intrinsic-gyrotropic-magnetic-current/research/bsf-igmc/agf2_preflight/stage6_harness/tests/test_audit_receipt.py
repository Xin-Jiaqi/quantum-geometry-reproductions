"""Synthetic receipts test only the honesty/consistency gate, not VASP physics."""
from __future__ import annotations
import copy
import importlib.util
import json
import math
import pathlib
import unittest

BASE = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("igmc_stage6_audit", BASE / "audit_receipt.py")
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class GateTests(unittest.TestCase):
    def template(self):
        return json.loads((BASE / "stage6_receipt.template.json").read_text(encoding="utf-8"))

    def complete_synthetic(self):
        x = self.template()
        x["data_origin"] = "real_dft"  # synthetic fixture shape only, never published as evidence
        c = x["campaign"]
        c.update({"source_structure_url":"https://example.org/structures", "method_source_url":"https://example.org/methods",
                  "source_structure_sha256":"a"*64, "protocol_sha256":"b"*64, "potcar_lock_sha256":"c"*64,
                  "functional":"EXAMPLE", "vdw_method":"EXAMPLE_NONE", "ag_hubbard_u_ev":0,
                  "physics_approved_by_user":True})
        for n in ("B1","B2"):
            s=x["states"][n]
            s.update({"poscar_path":f"structures/{n}/POSCAR", "poscar_sha256":"a"*64,
                      "relaxed_outcar_sha256":"d"*64, "num_fu":4, "total_energy_ev":-80.0,
                      "max_force_ev_per_ang":0.001, "scf_converged":True,"ionic_converged":True,
                      "magnetic_group_verified":True,"magnetic_axis":"x", "direct_gap_ev":0.5,
                      "indirect_gap_ev":0.3,"operators_gauge_validated":True})
        x["calculation_evidence"].update({"software_and_version":"EXAMPLE","kmesh_and_convergence_record":"synthetic",
                                             "energy_cutoff_convergence_record":"synthetic",
                                             "vacuum_convergence_record":"synthetic","magnetic_order_comparison_record":"synthetic",
                                             "wannier_operator_audit_record":"synthetic","field_frequency_thz":1,"zeeman_g_factor":2})
        return x

    def test_empty_source_never_passes(self):
        r = module.audit(self.template())
        self.assertEqual(r["status"], "BLOCKED_SOURCE_PROTOCOL")
        self.assertIsNone(r["delta_E_B2_minus_B1_meV_per_fu"])

    def test_synthetic_without_operator_is_blocked(self):
        x = self.complete_synthetic()
        x["states"]["B1"]["operators_gauge_validated"] = False
        self.assertEqual(module.audit(x)["status"], "BLOCKED_OPERATOR_AUDIT")

    def test_typeI_odd_dxx_dyy(self):
        x = self.complete_synthetic()
        x["states"]["B1"]["igmc"]={"D_xx":0.1,"D_yy":-0.12,"D_xx_abs_error":1e-4,"D_yy_abs_error":1e-4,"units":"sheet_alpha_P_per_B"}
        x["states"]["B2"]["igmc"]={"D_xx":-0.10001,"D_yy":0.12001,"D_xx_abs_error":1e-4,"D_yy_abs_error":1e-4,"units":"sheet_alpha_P_per_B"}
        r=module.audit(x)
        self.assertEqual(r["status"],"REVIEW_REQUIRED_NUMERICALLY_CONSISTENT")
        self.assertEqual(r["symmetry_parity_diagnostic"]["D_xx"]["status"],"CONSISTENT_NUMERICALLY")

    def test_typeII_even_dxx_dyy(self):
        x = self.complete_synthetic()
        x["campaign"].update({"magnetic_branch":"type_II","strain_percent":4})
        for n in ("B1","B2"):
            x["states"][n]["igmc"]={"D_xx":0.1,"D_yy":0.2,"D_xx_abs_error":1e-4,"D_yy_abs_error":1e-4,"units":"sheet_alpha_P_per_B"}
        self.assertEqual(module.audit(x)["status"],"REVIEW_REQUIRED_NUMERICALLY_CONSISTENT")

    def test_zero_signal_is_not_a_switch(self):
        x = module.parity_diagnostic(1e-8,-1e-8,-1,1e-6,1e-6)
        self.assertEqual(x["status"], "UNRESOLVED_ZERO_SIGNAL")

    def test_wrong_sign_requires_review(self):
        x = module.parity_diagnostic(0.1,0.1,-1,1e-4,1e-4)
        self.assertEqual(x["status"],"REQUIRES_PHYSICS_REVIEW")

    def test_bad_fu_and_nan(self):
        x = self.complete_synthetic()
        x["states"]["B1"]["num_fu"]=0
        x["states"]["B2"]["direct_gap_ev"]=float("nan")
        r = module.audit(x)
        self.assertEqual(r["status"],"BLOCKED_DFT_MAGNETIC")
        self.assertIsNone(r["delta_E_B2_minus_B1_meV_per_fu"])

    def test_no_real_dft_claim_from_synthetic_shape(self):
        x = self.complete_synthetic()
        x["data_origin"]="NOT_RUN"
        x["states"]["B1"]["igmc"]={"D_xx":0.1,"D_yy":0.2,"D_xx_abs_error":1e-4,"D_yy_abs_error":1e-4,"units":"sheet_alpha_P_per_B"}
        x["states"]["B2"]["igmc"]={"D_xx":-0.1,"D_yy":-0.2,"D_xx_abs_error":1e-4,"D_yy_abs_error":1e-4,"units":"sheet_alpha_P_per_B"}
        self.assertEqual(module.audit(x)["status"],"BLOCKED_REAL_DFT_PROVENANCE")

if __name__=='__main__':unittest.main()
