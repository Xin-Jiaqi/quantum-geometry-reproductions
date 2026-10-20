"""Physical and numeric contract for displacement-generated bilayer tunnelings."""
import unittest
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingBilayer, SlidingSpec, M_GLOBAL
from bsf_nmee.valley_quadrature import valley_quadrature
from scripts.sliding_static_scan import response_weighted


class TestMicroscopicSliding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec=SlidingSpec()
        cls.m=SlidingBilayer(cls.spec)
        cls.k=np.array([[0.21,0.54],[2.17,-1.92],[8.4,9.3]])

    def test_exact_hermiticity_velocity_and_mirror(self):
        m=self.m;p=SlidingBilayer(SlidingSpec(shift=(-.23,-.41)))
        h=m.matrices(self.k)
        self.assertLess(np.max(abs(h-h.conj().transpose(0,2,1))),1e-12)
        for axis in (0,1):
            eps=1e-5;dk=np.zeros_like(self.k);dk[:,axis]=eps
            finite=(m.matrices(self.k+dk)-m.matrices(self.k-dk))/(2*eps)
            self.assertLess(np.max(abs(m.matrices(self.k,axis)-finite)),1e-9)
            self.assertLess(np.max(abs(p.matrices(self.k,axis)-M_GLOBAL@m.matrices(self.k,axis)@M_GLOBAL.conj().T)),1e-11)
        self.assertLess(np.max(abs(p.matrices(self.k)-M_GLOBAL@h@M_GLOBAL.conj().T)),1e-11)

    def test_integer_bravais_shift_and_reciprocal_gauge(self):
        m=self.m
        shifted=SlidingBilayer(SlidingSpec(shift=(1.23,.41)))
        for axis in (None,0,1):
            self.assertLess(np.max(abs(shifted.matrices(self.k,axis)-m.matrices(self.k,axis))),1e-12)
        g=m.model().reciprocal[0]
        positions=np.concatenate((m.pos,m.pos+np.asarray(m.spec.shift)@m.a))
        U=np.diag(np.repeat(np.exp(1j*positions@g),2))
        for axis in (None,0,1):
            self.assertLess(np.max(abs(m.matrices(self.k+g,axis)-U.conj().T@m.matrices(self.k,axis)@U)),1e-11)

    def test_weighted_valley_quadrature_and_mirror_response(self):
        k,w,meta=valley_quadrature(self.m.model().reciprocal,exterior_n=41,radius=1.5,radial_n=35,angular_n=64)
        # Small-grid comparisons test symmetry, not amplitude convergence.
        self.assertLess(abs(meta['normalization_error']),.05)
        src=response_weighted(self.spec,k,w)
        mirror=response_weighted(SlidingSpec(shift=(-.23,-.41)),k,w)
        self.assertGreater(abs(src['response']['total_x']),1e-8)
        self.assertLess(abs(mirror['response']['total_x']+src['response']['total_x']),1e-10)
        self.assertLess(abs(mirror['response']['top_x']+src['response']['bottom_x']),1e-10)
        self.assertLess(abs(mirror['response']['bottom_x']+src['response']['top_x']),1e-10)
        self.assertLess(abs(src['total_vs_layers_error']),1e-10)

    def test_actual_toolkit_connection_feeds_independent_hamiltonian(self):
        try:
            from group_theory_operations import point_group_operations
            from bsf_nmee.connection import connect_interface
        except ImportError:
            self.skipTest('Install actual group-theory-operations-toolkit')
        mirror=next(op for op in point_group_operations('6/mmm') if op.name=='m_001')
        conn=connect_interface((.23,.41),mirror,expected_target=(.77,.59))
        self.assertTrue(conn.layer_exchanging)
        source=self.m
        target=SlidingBilayer(SlidingSpec(shift=conn.target_shift))
        # Since the model is periodic in fractional shift, (0.77,0.59)
        # is equivalent to (-0.23,-0.41), not a hand-assigned response.
        diff=np.max(abs(target.matrices(self.k)-M_GLOBAL@source.matrices(self.k)@M_GLOBAL.conj().T))
        self.assertLess(diff,1e-12)

    def test_decoupled_and_rashba_free_controls(self):
        k,w,_=valley_quadrature(self.m.model().reciprocal,41,1.5,35,64)
        zero_t=response_weighted(SlidingSpec(tunneling_eV=0),k,w)['response']['total_x']
        no_rashba=response_weighted(SlidingSpec(rashba_soc_eV=0),k,w)['response']['total_x']
        self.assertLess(abs(zero_t),1e-10)
        self.assertLess(abs(no_rashba),1e-10)

    def test_orbital_selectivity_without_claiming_material_prediction(self):
        k,w,_=valley_quadrature(self.m.model().reciprocal,41,1.5,35,64)
        full=response_weighted(SlidingSpec(),k,w)['response']['total_x']
        same=response_weighted(SlidingSpec(cross_sublattice=False),k,w)['response']['total_x']
        self.assertGreater(abs(full),100*abs(same))

if __name__=='__main__':unittest.main()
