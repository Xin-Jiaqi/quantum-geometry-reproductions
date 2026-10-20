"""Finite-group algebra and input-gate safety regression."""
from pathlib import Path
import json
import sys
import unittest
from tempfile import TemporaryDirectory
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from agf2_preflight import (symmetry_result,allowed_components,partner_sign,
                           validate_dft_contract,ahe_z_sign,R_MZ)
ROOT=Path(__file__).resolve().parents[1]

class TestAgf2Symmetry(unittest.TestCase):
 def test_mxprime_tensor_basis(self):
    self.assertEqual(allowed_components('C'),['xy','xz','yx','zx'])
    self.assertEqual(allowed_components('D'),['xx','yy','yz','zy','zz'])
 def test_paired_phase_selectivity(self):
    s=symmetry_result()['switches']
    self.assertEqual((s['type_I']['C_xy'],s['type_I']['D_xx'],s['type_I']['AHE_xy']),(-1,-1,1))
    self.assertEqual((s['type_II_strained']['C_xy'],s['type_II_strained']['D_xx'],s['type_II_strained']['AHE_xy']),(-1,1,-1))
 def test_connectors_reverse_polarization(self):
    self.assertEqual(symmetry_result()['switches']['type_I']['P_z'],-1)
    self.assertEqual(symmetry_result()['switches']['type_II_strained']['P_z'],-1)
 def test_ahe_internal_mxprime_allows_axial_z(self):
    from agf2_preflight import R_MX
    self.assertEqual(ahe_z_sign(R_MX,True),1)
 def test_original_model_axis_separation(self):
    self.assertEqual(partner_sign('m_z','D','xx'),-1)
    self.assertEqual(partner_sign("m_z'",'D','xx'),1)
 def test_antisymmetric_hall_tensor_not_conflated_with_dxx(self):
    self.assertNotEqual(symmetry_result()['switches']['type_I']['AHE_xy'],symmetry_result()['switches']['type_I']['D_xx'])
 def test_example_fails_closed_without_atomic_coordinates(self):
    data=json.loads((ROOT/'inputs/agf2_dft_contract.example.json').read_text())
    response=validate_dft_contract(data,ROOT/'inputs')
    self.assertEqual(response['status'],'BLOCKED_INPUTS')
    self.assertIn('MISSING_ATOMIC_STRUCTURE_B1',response['reasons'])
    self.assertIn('MISSING_ATOMIC_STRUCTURE_B2',response['reasons'])
    self.assertIn('MISSING_PHYSICS_hub_u_ag_ev',response['reasons'])
 def test_no_path_traversal(self):
    data=json.loads((ROOT/'inputs/agf2_dft_contract.example.json').read_text())
    data['structures']['B1']={'poscar':'../../sensitive/POSCAR','sha256':'a'*64}
    self.assertIn('INVALID_STRUCTURE_PATH_B1',validate_dft_contract(data,ROOT/'inputs')['reasons'])
 def test_strain_specific_campaign_must_not_mix_branches(self):
    data=json.loads((ROOT/'inputs/agf2_dft_contract.example.json').read_text())
    data['physics']['magnetic_order']='type_II'
    self.assertIn('TYPE_II_SOURCE_STRAIN_SHOULD_BE_FOUR_PERCENT',validate_dft_contract(data,ROOT/'inputs')['reasons'])
 def test_no_fake_poscar_passes_shape_gate(self):
    with TemporaryDirectory() as temp:
        root=Path(temp);path=root/'B1.vasp';path.write_text('nonsense\n')
        from agf2_preflight import validate_poscar_minimal
        self.assertIn('POSCAR_TOO_SHORT',validate_poscar_minimal(path))
 def test_duplicate_atomic_structure_rejected(self):
    with TemporaryDirectory() as temp:
     root=Path(temp)
     (root/'POSCAR').write_text('synthetic-only\n1.0\n3 0 0\n0 4 0\n0 0 25\nAg F\n2 4\nDirect\n0 0 0.4\n0.5 0.5 0.4\n0.25 0 0.4\n0.75 0 0.4\n0.25 0.5 0.4\n0.75 0.5 0.4\n')
     import hashlib
     sha=hashlib.sha256((root/'POSCAR').read_bytes()).hexdigest()
     payload={'material':'AgF2','task':'bsf_igmc_stage4','physics':{
      'functional':'PBE','vdw':'D3','hub_u_ag_ev':4,'soc':True,'vacuum_axis':'z',
      'magnetic_axis':'x','magnetic_order':'type_I','strain_percent':0,'kmesh':[9,9,1],
      'ecut_ev':500,'pseudopotential_lock_sha256':'a'*64},
      'structures':{'B1':{'poscar':'POSCAR','sha256':sha},'B2':{'poscar':'POSCAR','sha256':sha}}}
     res=validate_dft_contract(payload,root)
     self.assertEqual(res['status'],'BLOCKED_INPUTS')
     self.assertIn('STRUCTURAL_PARTNERS_IDENTICAL',res['reasons'])

if __name__=='__main__':unittest.main()
