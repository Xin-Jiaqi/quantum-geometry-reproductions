"""Independent symmetry, charge and derivative checks for the static protocol."""
from dataclasses import replace
import unittest
import numpy as np
from bsf_nmee.bilayer_audit import AuditSpec, AuditBilayer, integration_grid, symmetry_residuals
from bsf_nmee.static_observables import observables, flatten_hamiltonian, LAYER
from bsf_nmee.sliding_tight_binding import SlidingBilayer, SlidingSpec
from bsf_nmee.symmetry_scaling import static_inplane_matrix


class BilayerAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec = AuditSpec(shift=(.023,.041),tunneling_eV=.005)
        cls.k, cls.w, _ = integration_grid(AuditBilayer(cls.spec),29,25,48)

    def test_default_matches_frozen_hamiltonian_and_response(self):
        new = AuditBilayer(self.spec)
        oldspec = SlidingSpec(shift=self.spec.shift,tunneling_eV=.005)
        old = SlidingBilayer(oldspec)
        for axis in (None,0,1):
            np.testing.assert_allclose(new.matrices(self.k[:7],axis),
                                       old.matrices(self.k[:7],axis),atol=1e-14)
        actual = observables(new,self.k,self.w)['C']
        expected = static_inplane_matrix(oldspec,self.k,self.w,step_T=.125)['gamma']
        np.testing.assert_allclose(actual,expected,atol=2e-12)

    def test_physical_symmetries_and_intentional_breaking(self):
        residual = symmetry_residuals(self.spec)
        for name in ('C3_eV','mirror_pair_eV','time_reversal_eV'):
            self.assertLess(residual[name],1e-12)
        self.assertGreater(residual['graphene_inversion_eV'],1e-3)
        broken = symmetry_residuals(replace(self.spec,bond_anisotropy_eV=.008))
        self.assertGreater(broken['C3_eV'],1e-3)
        self.assertLess(broken['mirror_pair_eV'],1e-12)

    def test_graphene_direct_inversion_forbids_net_response(self):
        s = replace(self.spec,stacking='direct',staggered_eV=0.,
                    rashba_soc_eV=0.,mu_eV=0.)
        self.assertLess(symmetry_residuals(s)['graphene_inversion_eV'],1e-12)
        result = observables(AuditBilayer(s),self.k,self.w)
        self.assertLess(abs(result['delta_N_per_cell']),1e-13)
        self.assertLess(np.max(abs(np.array(result['C']))),1e-12)

    def test_electronic_dipole_sign_and_mirror(self):
        a = observables(AuditBilayer(self.spec),self.k,self.w)
        b = observables(AuditBilayer(replace(self.spec,shift=tuple(-np.array(self.spec.shift)))),self.k,self.w)
        np.testing.assert_allclose(a['C'],-np.array(b['C']),atol=1e-12)
        self.assertAlmostEqual(a['dipole_over_e_d_per_cell'],-b['dipole_over_e_d_per_cell'],delta=1e-14)
        self.assertAlmostEqual(a['dipole_over_e_d_per_cell'],-a['delta_N_per_cell']/2,delta=1e-16)

    def test_charge_hellmann_feynman_independent_energy_derivative(self):
        s = replace(self.spec,shift=(.23,.41),tunneling_eV=.028)
        model = AuditBilayer(s)
        result = observables(model,self.k,self.w,compute_response=False)
        # H(lambda)=H+lambda*(P_top-P_bottom); dE/dlambda = layer imbalance.
        step = 2e-6; energy_derivative = 0.
        for start in range(0,len(self.k),1024):
            h = model.matrices(self.k[start:start+1024])
            ep = np.linalg.eigvalsh(h+step*LAYER)[:,:4].sum(axis=1)
            em = np.linalg.eigvalsh(h-step*LAYER)[:,:4].sum(axis=1)
            energy_derivative += self.w[start:start+1024]@((ep-em)/(2*step))
        energy_derivative *= abs(np.linalg.det(model.a))
        self.assertAlmostEqual(energy_derivative,result['delta_N_per_cell'],delta=2e-10)

    def test_flattening_has_fixed_spectrum_correct_velocity_and_same_charge(self):
        model = AuditBilayer(self.spec); k=self.k[:9]; step=1e-5; gap=.03
        h=model.matrices(k); v=[model.matrices(k,axis=j) for j in range(2)]
        hf,vf=flatten_hamiltonian(h,v,gap)
        np.testing.assert_allclose(np.linalg.eigvalsh(hf),np.tile([-gap/2]*4+[gap/2]*4,(len(k),1)),atol=1e-13)
        for j in range(2):
            shift=np.eye(2)[j]*step
            hp,_=flatten_hamiltonian(model.matrices(k+shift),[],gap)
            hm,_=flatten_hamiltonian(model.matrices(k-shift),[],gap)
            np.testing.assert_allclose(vf[j],(hp-hm)/(2*step),atol=1e-9)
        _,u=np.linalg.eigh(h); _,uf=np.linalg.eigh(hf)
        np.testing.assert_allclose(u[:,:,:4]@u[:,:,:4].conj().swapaxes(1,2),
                                   uf[:,:,:4]@uf[:,:,:4].conj().swapaxes(1,2),atol=1e-13)

    def test_global_hopping_sign_is_a_layer_gauge(self):
        a = observables(AuditBilayer(self.spec),self.k,self.w)
        b = observables(AuditBilayer(replace(self.spec,same_scale=-1.,cross_scale=-1.)),self.k,self.w)
        np.testing.assert_allclose(a['C'],b['C'],atol=1e-12)
        self.assertAlmostEqual(a['delta_N_per_cell'],b['delta_N_per_cell'],delta=1e-14)


if __name__ == '__main__':
    unittest.main()
