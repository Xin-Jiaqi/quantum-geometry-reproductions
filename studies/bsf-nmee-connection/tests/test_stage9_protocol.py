"""Physical coordinates, microscopic selection rules and response cross-checks."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.bilayer_audit import AuditSpec, AuditBilayer, integration_grid, symmetry_residuals
from bsf_nmee.bilayer_protocol import fraction120_to60, ab_ba_shift, half_filled_spec, prepared_B_response
from bsf_nmee.static_observables import observables, LAYER


class ProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.k, cls.w, _ = integration_grid(AuditBilayer(), 29, 25, 48)
        cls.spec = half_filled_spec(AuditSpec(shift=(.023, .041), tunneling_eV=.005), cls.k)

    def test_coordinates_and_lattice_periodicity(self):
        model = AuditBilayer(self.spec)
        basis120 = np.array([model.a[0], model.a[1]-model.a[0]])
        uv = np.array([1/3, 2/3])
        np.testing.assert_allclose(fraction120_to60(uv)@model.a, uv@basis120, atol=1e-15)
        length = np.linalg.norm((np.array(ab_ba_shift(1))-ab_ba_shift(0))@model.a)
        self.assertAlmostEqual(length, self.spec.bond_nm)
        # Embedded-gauge tunneling sum relabels real-space translations.
        periodic = AuditBilayer(replace(self.spec, shift=(1.023, .041)))
        np.testing.assert_allclose(model.matrices(self.k[:9]), periodic.matrices(self.k[:9]), atol=1e-12)

    def test_no_staggering_restores_inversion_even_with_rashba(self):
        spec = half_filled_spec(replace(self.spec, staggered_eV=0.), self.k)
        self.assertLess(symmetry_residuals(spec)['graphene_inversion_eV'], 1e-12)
        result = observables(AuditBilayer(spec), self.k, self.w)
        self.assertLess(np.max(np.abs(result['C'])), 1e-12)
        self.assertLess(abs(result['delta_N_per_cell']), 1e-13)

    def test_spin_conservation_and_rashba_reversal(self):
        reference = observables(AuditBilayer(self.spec), self.k, self.w)
        reversed_r = observables(AuditBilayer(replace(self.spec, rashba_soc_eV=-.02)), self.k, self.w)
        np.testing.assert_allclose(reference['C'], -np.array(reversed_r['C']), atol=1e-12)
        self.assertAlmostEqual(reference['delta_N_per_cell'], reversed_r['delta_N_per_cell'], delta=1e-13)
        no_r = half_filled_spec(replace(self.spec, rashba_soc_eV=0.), self.k)
        result = observables(AuditBilayer(no_r), self.k, self.w)
        self.assertLess(np.max(np.abs(result['C'])), 1e-12)

    def test_same_sublattice_tunneling_has_pointwise_equal_layer_charge(self):
        model = AuditBilayer(replace(self.spec, cross_scale=0.))
        _, states = np.linalg.eigh(model.matrices(self.k[::53]))
        occupied = states[:, :, :4]
        imbalance = np.einsum('kav,ab,kbv->k', occupied.conj(), LAYER, occupied).real
        np.testing.assert_allclose(imbalance, 0., atol=2e-13)

    def test_dynamic_static_limit_and_layer_additivity(self):
        model = AuditBilayer(self.spec)
        dynamic = prepared_B_response(model, self.k, self.w, [0.], eta_eV=1e-7)
        static = observables(model, self.k, self.w, step_T=.0625)
        np.testing.assert_allclose(dynamic['real'][0], static['C'], rtol=1e-7, atol=1e-13)
        layers = [prepared_B_response(model, self.k, self.w, [0.], eta_eV=1e-7,
                                     observation=s)['real'][0] for s in ('top', 'bottom')]
        np.testing.assert_allclose(np.sum(layers, axis=0), dynamic['real'][0], atol=2e-13)

    def test_ab_ba_isotropy_and_local_quadratic_scalar(self):
        ab = np.array(ab_ba_shift(0))
        results = {}
        for q in (0., -.01, .01, -.02, .02):
            spec = half_filled_spec(replace(self.spec, shift=tuple(ab+q*np.array([.23,.41]))), self.k)
            results[q] = observables(AuditBilayer(spec), self.k, self.w)
        c = np.array(results[0.]['C'])
        np.testing.assert_allclose(c, np.eye(2)*np.trace(c)/2, atol=1e-12)
        for key in ('dipole_over_e_d_per_cell', 'trace'):
            even = [results[q][key]+results[-q][key]-2*results[0.][key] for q in (.01,.02)]
            self.assertAlmostEqual(np.log2(abs(even[1]/even[0])), 2., delta=.01)
        odd = [results[q]['difference']-results[-q]['difference'] for q in (.01,.02)]
        self.assertAlmostEqual(np.log2(abs(odd[1]/odd[0])), 1., delta=.01)


if __name__ == '__main__':
    unittest.main()
