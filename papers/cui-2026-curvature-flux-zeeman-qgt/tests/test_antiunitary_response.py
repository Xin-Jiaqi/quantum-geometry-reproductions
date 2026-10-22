"""Round15 independent antiunitary spinor, clean Kubo and Liouville tests."""
import json,unittest
from pathlib import Path
import numpy as np
from flux_lab.core import PAULI
from flux_lab.causal_response import interband_kubo,interband_liouville,source_four_sector_expansion
from flux_lab.antiunitary_response import (examples,toy_bloch,operator_residual,
 qgt_band_pair,qgt_symmetry_residual,kubo_symmetry_residual,
 pair_frequency_parities,wrong_retarded_branch_error,toy_spin_orbital,
 fourband_antiunitary_residual)
K=np.array([.43,-.38]);Z=.07+.014j
HERE=Path(__file__).resolve().parents[1]

def bz(kind,n=16,w=.035):
    ks=-np.pi+(np.arange(n)+.5)*2*np.pi/n
    p=np.zeros((2,3),complex);m=p.copy()
    for x in ks:
        for y in ks:
            h,v,s=toy_bloch([x,y],kind)
            p+=interband_kubo(h,v,s,z=w,mu=.1,temperature=.16)[0]
            m+=interband_kubo(h,v,s,z=-w,mu=.1,temperature=.16)[0]
    odd=(p-m)/(2*n*n);even=(p+m)/(2*n*n)
    return odd,even

class MicroscopicAntiunitary(unittest.TestCase):
    def test_three_exact_spinor_groups(self):
        for name in examples():
            self.assertLess(max(operator_residual(name,K).values()),2e-14)
    def test_tr_has_k_to_minus_k_all_spin_odd(self):
        s=examples()['TR']
        np.testing.assert_allclose(s.R,-np.eye(2))
        np.testing.assert_array_equal(s.spin_parity,[-1,-1,-1])
    def test_c2zt_has_local_k_even_inplane_spin_and_odd_z(self):
        s=examples()['C2zT']
        np.testing.assert_allclose(s.R,np.eye(2))
        np.testing.assert_array_equal(s.spin_parity,[1,1,-1])
    def test_mxt_momentum_mirror_spin_signs(self):
        s=examples()['mxT']
        np.testing.assert_allclose(s.R,np.diag([1,-1]))
        np.testing.assert_array_equal(s.spin_parity,[-1,1,1])
    def test_explicit_velocity_derivative_not_guessed(self):
        step=2e-6
        for name in examples():
            _,v,_=toy_bloch(K,name)
            for a in range(2):
                dx=np.zeros(2);dx[a]=step
                df=(toy_bloch(K+dx,name)[0]-toy_bloch(K-dx,name)[0])/(2*step)
                np.testing.assert_allclose(v[a],df,atol=2e-10)
    def test_exact_qgt_pair_transformation_at_generic_points(self):
        for name in examples():
            for k in (K,np.array([.6,.2])):
                self.assertLess(qgt_symmetry_residual(name,k),4e-14)
    def test_tr_real_zeeman_t_odd_imag_even(self):
        g=examples()['TR'];h,v,s=toy_bloch(K,'TR')
        H,V,S=toy_bloch(g.R@K,'TR')
        t=qgt_band_pair(h,v,s);u=qgt_band_pair(H,V,S)
        np.testing.assert_allclose(u.real,-t.real,atol=1e-13)
        np.testing.assert_allclose(u.imag,t.imag,atol=1e-13)
    def test_c2zt_local_real_inplane_zeeman_t_vanishes(self):
        h,v,s=toy_bloch(K,'C2zT');t=qgt_band_pair(h,v,s)
        self.assertLess(np.max(abs(t[:,:2].real)),2e-14)
        self.assertLess(np.max(abs(t[:,2].imag)),2e-14)
        self.assertGreater(np.max(abs(t[:,:2].imag)),1e-3)
        self.assertGreater(np.max(abs(t[:,2].real)),1e-3)
    def test_full_complex_frequency_kubo_symmetry_not_approximation(self):
        for name in examples():
            for z in (.07,-.07,.07+.014j,-.05+.017j):
                q=kubo_symmetry_residual(name,K,z=z)
                self.assertLess(q['chi'],1e-12)
                self.assertLess(q['alpha'],1e-12)
    def test_orbital_basis_liouville_independent_of_spectral_oracle(self):
        for name in examples():
            q=kubo_symmetry_residual(name,K,z=Z,method='liouville')
            self.assertLess(max(q['chi'],q['alpha']),2e-12)
            h,v,s=toy_bloch(K,name)
            a=interband_kubo(h,v,s,z=Z,mu=.1,temperature=.16)
            b=interband_liouville(h,v,s,z=Z,mu=.1,temperature=.16)
            for x,y in zip(a,b):np.testing.assert_allclose(x,y,atol=2e-12)
    def test_retarded_at_negative_frequency_is_not_advanced_minus_z(self):
        for name in examples():
            self.assertGreater(wrong_retarded_branch_error(name,K),1e-4)
    def test_tr_clean_pair_forbids_odd_frequency_but_not_even(self):
        q=pair_frequency_parities('TR',K)
        self.assertLess(np.max(abs(q['odd'])),1e-12)
        self.assertGreater(np.max(abs(q['even'])),.2)
    def test_c2zt_clean_inplane_even_outofplane_odd(self):
        q=pair_frequency_parities('C2zT',K)
        self.assertLess(np.max(abs(q['odd'][:,:2])),1e-13)
        self.assertLess(np.max(abs(q['even'][:,2])),1e-13)
        self.assertGreater(np.max(abs(q['odd'][:,2])),1e-3)
        self.assertGreater(np.max(abs(q['even'][:,:2])),.1)
    def test_mirror_tr_clean_parity_matrix(self):
        q=pair_frequency_parities('mxT',K)
        self.assertLess(np.max(abs(q['even'][:,:2][np.diag_indices(2)])),1e-13)
        self.assertLess(max(abs(q['odd'][0,1]),abs(q['odd'][1,0])),1e-13)
        self.assertGreater(abs(q['odd'][1,1]),1e-3)
        self.assertGreater(abs(q['even'][1,0]),.5)
    def test_odd_even_low_frequency_correspond_to_source_four_sectors(self):
        for name in examples():
            k=toy_bloch(K,name);q=source_four_sector_expansion(*k,omega_source=.008,mu=.1,temperature=.16)
            p=source_four_sector_expansion(*toy_bloch(examples()[name].R@K,name),
                                      omega_source=.008,mu=.1,temperature=.16)
            linear=q['linear_coefficient']+p['linear_coefficient']
            quadratic=q['quadratic_coefficient']+p['quadratic_coefficient']
            sign=np.diag(examples()[name].R)[:,None]*examples()[name].spin_parity[None,:]
            self.assertLess(np.max(abs(linear[sign==1])),3e-13)
            if np.any(sign==-1):self.assertLess(np.max(abs(quadratic[sign==-1])),3e-13)
    def test_symmetry_breaking_negative_control(self):
        g=examples()['C2zT'];h,_,_=toy_bloch(K,'C2zT')
        broken=h+.2*PAULI[2]
        self.assertGreater(np.max(abs(broken-g.U@broken.conj()@g.U.conj().T)),.3)
    def test_bz_unweighted_pair_cancellations_and_nonzero_allowed(self):
        for name in examples():
            odd,even=bz(name)
            sign=np.diag(examples()[name].R)[:,None]*examples()[name].spin_parity[None,:]
            self.assertLess(np.max(abs(odd[sign==1])),1e-12)
            if np.any(sign==-1):self.assertLess(np.max(abs(even[sign==-1])),1e-12)
            self.assertGreater(np.max(abs(even)),1e-3)
        self.assertGreater(np.max(abs(bz('C2zT')[0])),5e-4)
    def test_round15_frozen_reference_reproducible(self):
        q=json.loads((HERE/'results/frozen/round15_antiunitary_microscopic.json').read_text())
        for name in examples():
            d=q['models'][name]
            self.assertAlmostEqual(qgt_symmetry_residual(name,K),d['qgt_pair_max_residual'],delta=2e-14)
            odd,even=bz(name)
            self.assertAlmostEqual(float(np.max(abs(odd))),d['BZ16']['max_odd'],delta=1e-12)
            self.assertAlmostEqual(float(np.max(abs(even))),d['BZ16']['max_even'],delta=1e-12)
    def test_four_band_genuine_spin_orbital_mixing_antiunitary(self):
        for name in examples():
            d=fourband_antiunitary_residual(name,K)
            self.assertLess(max(d.values()),2e-12)
    def test_four_band_independent_liouville(self):
        for name in examples():
            h,v,s=toy_spin_orbital(K,name)
            a=interband_kubo(h,v,s,z=Z,mu=.1,temperature=.16)
            b=interband_liouville(h,v,s,z=Z,mu=.1,temperature=.16)
            for x,y in zip(a,b):np.testing.assert_allclose(x,y,atol=2e-12)

if __name__=='__main__':unittest.main()
