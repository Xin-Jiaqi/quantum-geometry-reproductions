"""Independent magnetic group closure, exact tensor nullspace and microscopic controls."""
import unittest
import numpy as np
from flux_lab.core import model,qgt_eigen,PAULI
from flux_lab.magnetic_symmetry import (load_groups,SECTORS,matrices,closure,
  projector_basis,independent_nullity,allowed,audit_table,transverse_symmetric_allowed)

class MagneticSymmetry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.groups=load_groups()

    def test_31_groups_and_group_orders(self):
        expected={1:1,2:2,6:2,7:4,8:2,9:2,10:4,11:2,20:4,21:8,22:4,23:4,
          29:4,30:8,31:4,44:8,45:16,46:8,47:8,60:3,61:6,68:6,69:12,
          70:6,76:6,77:12,78:6,91:12,92:24,93:12,94:12}
        self.assertEqual({g['number']:len(g['operations']) for g in self.groups.values()},expected)

    def test_closed_orthogonal_groups_and_time_reversal_xor(self):
        for name,g in self.groups.items():
            for R,t in g['operations']:
                np.testing.assert_allclose(R@R.T,np.eye(2),atol=1e-11)
                self.assertAlmostEqual(abs(np.linalg.det(R)),1.,places=8)
                for S,u in g['operations']:
                    self.assertTrue(any(v==(t^u) and np.linalg.norm(Q-R@S)<1e-8 for Q,v in g['operations']),name)

    def test_reynolds_projector_and_independent_nullspaces_all_124(self):
        for name,g in self.groups.items():
            for c,(kind,eta) in SECTORS.items():
                pr,_=projector_basis(g['operations'],kind,eta)
                nc=independent_nullity(g['operations'],kind,eta)
                self.assertEqual(pr,nc,(name,c))

    def test_magnetic_time_reversal_gray_group_candidate(self):
        self.assertEqual(tuple(allowed(self.groups["11'"],x) for x in SECTORS),(0,3,1,0))

    def test_nonmagnetic_rotation_and_mirror_candidate(self):
        for name in ('3m','4mm','6mm'):
            self.assertEqual(tuple(allowed(self.groups[name],x) for x in SECTORS),(0,0,1,1),name)

    def test_invariant_symmetric_metric_not_transverse_symmetric_response(self):
        for name in ['3','3m\'','4','4m\'m\'','6','6m\'m\'']:
            self.assertGreater(allowed(self.groups[name],'gN'),0,name)
            self.assertFalse(transverse_symmetric_allowed(self.groups[name]),name)

    def test_source_table_discrepancies_are_recorded_not_overridden(self):
        differences={(r['group'],s) for r in audit_table() for s in r['disagreements']}
        expected={(name,sector) for name in ["2'","m'","m'm2'"] for sector in ('gN','omegaA')}
        self.assertEqual(differences,expected)

    def test_all_normal_curvature_candidate_allowed(self):
        self.assertTrue(all(allowed(g,'omegaN')==1 for g in self.groups.values()))

    def test_explicit_c2z_time_reversal_spinor_geometry(self):
        # Physical spinor A=sigma_x*K with A H(k) A^-1=H(k), k fixed.
        x,y=.29,-.21;d,dk=model(x,y,m=0.,vx=1.2,vy=.85)
        h=np.einsum('a,aij->ij',d,PAULI)
        U=PAULI[0]
        np.testing.assert_allclose(U@h.conj()@U.conj().T,h,atol=1e-13)
        q=qgt_eigen(d,dk,-1)
        self.assertLess(np.max(abs(q['T'][:2,:2].real)),1e-12)
        self.assertLess(abs(q['omegaA'][2]),1e-12)
        self.assertGreater(np.linalg.norm(q['omegaA'][:2]),.1)

    def test_closure_generator_order_independent(self):
        a=closure([('C4',0),('Mx',1)])
        b=closure([('Mx',1),('C4',0)])
        self.assertEqual(len(a),len(b))
        for R,t in a:
            self.assertTrue(any(u==t and np.linalg.norm(S-R)<1e-8 for S,u in b))

    def test_groups_are_finite_and_bad_generator_fails(self):
        with self.assertRaises(ValueError):closure([('C7',0)])
        with self.assertRaises(ValueError):closure([('Mx',2)])

    def test_source_has_no_implied_finite_frequency_microscopic_oracle(self):
        # A group invariant tensor is a necessary condition only, not proof
        # that a material has a nonzero dynamical response or even that the
        # *assumed* antiunitary Jahn character is correct for every channel.
        z=self.groups['4mm']
        self.assertEqual(allowed(z,'gN'),0)
        self.assertEqual(allowed(z,'omegaA'),1)
        # Projector dimension alone does not produce a conductivity value.
        self.assertIsNone(projector_basis(z['operations'],'antisymmetric',-1)[1])

if __name__=='__main__':unittest.main()
