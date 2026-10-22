"""v70 multi-frequency ordered L/eta limit gates; high-resolution values are regenerated."""
import json,pathlib,unittest
import numpy as np
from qgmap.finite_lead_lehmann_v69 import Spec,lehmann_response
from scripts.reproduce_v70_ordered_limits import create,production_value

ROOT=pathlib.Path(__file__).resolve().parents[1]
FROZEN=json.loads((ROOT/'results/frozen/v70_ordered_limits.json').read_text())
FREQ=('0.20','0.29','0.40')

class V70OrderedLimits(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.quick=create(False)
        cls.medium={};cls.production={}
        for w in (.20,.29,.40):
            key=f'{w:.2f}'
            cls.medium[key]=lehmann_response(Spec(omega=w,length=500,eta=.01))['sigma_device']
            if key=='0.29':
                cls.production[key]=production_value(w,length=1000,eta=.005)
            else:
                z=FROZEN['production_sequence'][key]['sigma_device']
                cls.production[key]=complex(z['real'],z['imag'])
        cls.endcap={}
        for cap in (-1.1,0.,1.1):
            cls.endcap[f'{cap:+.1f}']=lehmann_response(Spec(omega=.29,length=500,eta=.01,end_cap=cap))['sigma_device']

    def test_01_three_frequencies_present(self):
        self.assertEqual(tuple(self.quick['quick_sequence'].keys()),FREQ)

    def test_02_three_continuum_references_present(self):
        self.assertEqual(tuple(self.quick['keldysh_reference'].keys()),FREQ)

    def test_03_all_reference_values_finite(self):
        for k in FREQ:
            z=self.quick['keldysh_reference'][k]['sigma_device']
            self.assertTrue(np.isfinite([z['real'],z['imag']]).all())

    def test_04_quick_sequence_starts_at_L250(self):
        for k in FREQ:self.assertEqual(self.quick['quick_sequence'][k][0]['length'],250)

    def test_05_quick_sequence_second_point_L350(self):
        for k in FREQ:self.assertEqual(self.quick['quick_sequence'][k][1]['length'],350)

    def test_06_joint_limit_uses_correct_eta(self):
        for k in FREQ:
            self.assertEqual(self.quick['quick_sequence'][k][0]['eta'],.02)
            self.assertEqual(self.quick['quick_sequence'][k][1]['eta'],.02)

    def test_07_L500_eta001_values_are_finite(self):
        for k in FREQ:self.assertTrue(np.isfinite([self.medium[k].real,self.medium[k].imag]).all())

    def test_08_L1000_eta0005_values_are_finite(self):
        for k in FREQ:self.assertTrue(np.isfinite([self.production[k].real,self.production[k].imag]).all())

    def test_09_L500_reduces_error_at_omega020(self):
        ref=self.quick['keldysh_reference']['0.20']['sigma_device'];z=complex(ref['real'],ref['imag'])
        self.assertLess(abs(self.medium['0.20']-z),self.quick['quick_sequence']['0.20'][0]['absolute_error'])

    def test_10_L500_reduces_error_at_omega029(self):
        self.assertLess(abs(self.medium['0.29']-complex(-.005292649615898412,-.001194257399094567)),
                        self.quick['quick_sequence']['0.29'][0]['absolute_error'])

    def test_11_L500_reduces_error_at_omega040(self):
        ref=self.quick['keldysh_reference']['0.40']['sigma_device'];z=complex(ref['real'],ref['imag'])
        self.assertLess(abs(self.medium['0.40']-z),self.quick['quick_sequence']['0.40'][0]['absolute_error'])

    def test_12_L1000_reduces_error_at_omega020(self):
        ref=self.quick['keldysh_reference']['0.20']['sigma_device'];z=complex(ref['real'],ref['imag'])
        self.assertLess(abs(self.production['0.20']-z),abs(self.medium['0.20']-z))

    def test_13_L1000_reduces_error_at_omega029(self):
        ref=self.quick['keldysh_reference']['0.29']['sigma_device'];z=complex(ref['real'],ref['imag'])
        self.assertLess(abs(self.production['0.29']-z),abs(self.medium['0.29']-z))

    def test_14_L1000_reduces_error_at_omega040(self):
        ref=self.quick['keldysh_reference']['0.40']['sigma_device'];z=complex(ref['real'],ref['imag'])
        self.assertLess(abs(self.production['0.40']-z),abs(self.medium['0.40']-z))

    def test_15_production_error_bounded_omega020(self):
        self.assertLess(FROZEN['production_sequence']['0.20']['absolute_error'],1.5e-4)

    def test_16_production_error_bounded_omega029(self):
        self.assertLess(FROZEN['production_sequence']['0.29']['absolute_error'],2.5e-4)

    def test_17_production_error_bounded_omega040(self):
        self.assertLess(FROZEN['production_sequence']['0.40']['absolute_error'],2.5e-4)

    def test_18_production_residual_not_zero_omega020(self):
        self.assertGreater(FROZEN['production_sequence']['0.20']['absolute_error'],1e-5)

    def test_19_production_residual_not_zero_omega029(self):
        self.assertGreater(FROZEN['production_sequence']['0.29']['absolute_error'],1e-5)

    def test_20_production_residual_not_zero_omega040(self):
        self.assertGreater(FROZEN['production_sequence']['0.40']['absolute_error'],1e-5)

    def test_21_relative_error_is_under_five_percent_omega020(self):
        self.assertLess(FROZEN['production_sequence']['0.20']['relative_error'],.05)

    def test_22_relative_error_is_under_five_percent_omega029(self):
        self.assertLess(FROZEN['production_sequence']['0.29']['relative_error'],.05)

    def test_23_relative_error_is_under_five_percent_omega040(self):
        self.assertLess(FROZEN['production_sequence']['0.40']['relative_error'],.05)

    def test_24_intermediate_L350_is_not_assumed_monotonic(self):
        e250=self.quick['quick_sequence']['0.20'][0]['absolute_error']
        e350=self.quick['quick_sequence']['0.20'][1]['absolute_error']
        self.assertGreater(e350,e250)

    def test_25_intermediate_L350_improves_at_omega029(self):
        e=self.quick['quick_sequence']['0.29']
        self.assertLess(e[1]['absolute_error'],e[0]['absolute_error'])

    def test_26_intermediate_L350_improves_at_omega040(self):
        e=self.quick['quick_sequence']['0.40']
        self.assertLess(e[1]['absolute_error'],e[0]['absolute_error'])

    def test_27_endcap_minus_plus_symmetry(self):
        self.assertLess(abs(self.endcap['-1.1']-self.endcap['+1.1']),2e-8)

    def test_28_endcap_spread_is_small(self):
        center=self.endcap['+0.0']
        self.assertLess(max(abs(v-center) for v in self.endcap.values()),3e-5)

    def test_29_endcap_error_smaller_than_continuum_residual(self):
        ref=self.quick['keldysh_reference']['0.29']['sigma_device'];z=complex(ref['real'],ref['imag'])
        spread=max(abs(v-self.endcap['+0.0']) for v in self.endcap.values())
        residual=abs(self.medium['0.29']-z)
        self.assertLess(spread,.1*residual)

    def test_30_medium_values_regenerate_frozen_omega029(self):
        z=FROZEN['medium_sequence']['0.29']['sigma_device']
        self.assertLess(abs(self.medium['0.29']-complex(z['real'],z['imag'])),1e-7)

    def test_31_omega020_isolated_production_record_finite(self):
        z=FROZEN['production_sequence']['0.20']['sigma_device']
        self.assertTrue(np.isfinite([z['real'],z['imag']]).all())
        self.assertLess(FROZEN['production_sequence']['0.20']['absolute_error'],1.5e-4)

    def test_32_production_values_regenerate_frozen_omega029(self):
        z=FROZEN['production_sequence']['0.29']['sigma_device']
        self.assertLess(abs(self.production['0.29']-complex(z['real'],z['imag'])),1e-7)

    def test_33_omega040_isolated_production_record_finite(self):
        z=FROZEN['production_sequence']['0.40']['sigma_device']
        self.assertTrue(np.isfinite([z['real'],z['imag']]).all())
        self.assertLess(FROZEN['production_sequence']['0.40']['absolute_error'],2.5e-4)

    def test_34_spectral_eta_is_not_labeled_as_reservoir_gamma(self):
        self.assertIn('spectral resolution',FROZEN['limit_path']['eta_role'])

    def test_35_original_paper_status_remains_partial(self):
        self.assertIn('PARTIAL',FROZEN['scientific_status'])

if __name__=='__main__':unittest.main()
