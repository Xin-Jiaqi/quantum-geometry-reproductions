"""Scope: registration, file evidence and fail-closed design contract only."""
import copy
import unittest
from pathlib import Path
from integration.physics_interface import (load_registry, validate_registry,
                                            get_module, plan_module)
from integration.physics_interface.registry import InvalidRegistry
from integration.physics_interface.contracts import ResponseRequestV1, evaluate

ROOT=Path(__file__).resolve().parents[3]

class RegistryTests(unittest.TestCase):
    def setUp(self):
        self.reg=load_registry()
    def test_four_known_papers(self):
        self.assertEqual({x["id"] for x in self.reg["modules"]},
           {"jia2026-nmee","xiang2025-igmc","xiang2025-spin-qgt","xiang2026-qgmap"})
    def test_safe_live_files_and_evidence(self):
        result=validate_registry(self.reg,repo_root=ROOT)
        self.assertEqual(result["module_count"],4)
        self.assertEqual(result["status"],"METADATA_VALID_ONLY")
    def test_each_paper_is_partial(self):
        self.assertTrue(all(x["science_status"]=="PARTIAL" for x in self.reg["modules"]))
    def test_all_protocol_ids_unique_and_scoped(self):
        ids=[p["id"] for m in self.reg["modules"] for p in m["protocols"]]
        self.assertEqual(len(ids),len(set(ids)))
        self.assertGreaterEqual(len(ids),10)
    def test_no_physics_dispatch_claim(self):
        with self.assertRaises(NotImplementedError):
            evaluate(None,None)
    def test_unknown_protocol_refused(self):
        with self.assertRaises(InvalidRegistry):
            plan_module("xiang2025-spin-qgt","nmee.printed",self.reg)
    def test_known_plan_is_nonexecuting(self):
        result=plan_module("xiang2025-spin-qgt","spinqgt.finite_eta",self.reg)
        self.assertEqual(result["execution"],"NOT_EXECUTED")
        self.assertEqual(result["cwd"],"papers/xiang-2025-spin-quantum-geometry")
    def test_mutated_duplicate_rejected(self):
        reg=copy.deepcopy(self.reg);reg["modules"][1]["id"]=reg["modules"][0]["id"]
        with self.assertRaises(InvalidRegistry):
            validate_registry(reg,repo_root=ROOT)
    def test_parent_path_rejected(self):
        reg=copy.deepcopy(self.reg);reg["modules"][0]["directory"]="papers/../../tmp"
        with self.assertRaises(InvalidRegistry):
            validate_registry(reg,repo_root=ROOT,check_files=False)
    def test_source_revision_required_for_study(self):
        from dataclasses import replace
        r=ResponseRequestV1(method_id="spinqgt.finite_eta",protocol_id="spinqgt.finite_eta",
           observable_id="Jx_sy", observable_basis_id="orbital-A",hamiltonian_basis_id="orbital-A",
           occupation="fermi-dirac",temperature_value=.1,temperature_unit="eV",k_convention="cartesian",
           k_unit="1/angstrom",energy_unit="eV",eta_value=.01,eta_kind="adiabatic_regulator",
           source_revision="UNPINNED")
        with self.assertRaisesRegex(ValueError,"pin"):
            r.validate()
        replace(r,source_revision="0123456789abcdef").validate()
    def test_operator_basis_guard(self):
        from dataclasses import replace
        r=ResponseRequestV1(method_id="igmc.low_frequency",protocol_id="igmc.low_frequency",
            observable_id="current",observable_basis_id="A",hamiltonian_basis_id="B",
            occupation="fermi-dirac",temperature_value=0.,temperature_unit="K",
            k_convention="fractional_reciprocal",k_unit="fractional",
            energy_unit="eV",source_revision="0123456789abcdef")
        with self.assertRaisesRegex(ValueError,"same"):
            r.validate()
        replace(r,observable_basis_id="B").validate()
    def test_laurent_is_not_generic_intrinsic(self):
        r=ResponseRequestV1(method_id="spinqgt.laurent_C",protocol_id="spinqgt.intrinsic_tau0",
            observable_id="spin",observable_basis_id="B",hamiltonian_basis_id="B",
            occupation="fermi-dirac",temperature_value=.1,temperature_unit="K",
            k_convention="cartesian",k_unit="1/angstrom",energy_unit="eV",
            eta_value=.02,eta_kind="adiabatic_regulator",source_revision="abcdef123456")
        with self.assertRaisesRegex(ValueError,"not validated"):
            r.validate()

if __name__=="__main__":
    unittest.main()
