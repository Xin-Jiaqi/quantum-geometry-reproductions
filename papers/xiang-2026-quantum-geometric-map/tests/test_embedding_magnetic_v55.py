"""v55 15 independent physical frame, embedding and magnetic Peierls tests."""
import unittest
import numpy as np
from scipy.linalg import expm
from qgmap.orbital_geometry_v54 import periodic_three_band_model,band_geometry
from qgmap.embedding_magnetic_v55 import (
 covariant_velocity,frame_transform,connected_band_geometry,
 orbital_center_connection,peierls_field,field_gauge_unitary,
 magnetic_vertex_frame_law,plaquette_wilson)

class EmbeddingAndPeierlsV55(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        k,H,V=periodic_three_band_model(12)
        cls.H=H[4,7];cls.V=V[:,4,7];cls.k=(float(k[4]),float(k[7]))
        cls.A=orbital_center_connection([[0.,0.],[.19,-.25],[-.32,.47]])
        cls.K=np.array([[.2,.4+.1j,0],[.4-.1j,-.13,.2j],[0,-.2j,.36]])
        x,y=cls.k
        cls.W=expm(1j*(.35*np.sin(x)-.21*np.cos(y))*cls.K)
        cls.dW=np.stack([1j*a*cls.K@cls.W for a in (.35*np.cos(x),.21*np.sin(y))])
        cls.site=np.array([[0,0],[1.25,0],[1.25,.8],[0,.8]],float)
        cls.H0=np.array([[.2,1.,.3,.8],[1.,-.1,1.2,.1],
                         [.3,1.2,.05,.9],[.8,.1,.9,.3]],complex)
        cls.B=.32;cls.q=-.87

    def test_01_zero_connection_gives_matched_analytic_velocity(self):
        np.testing.assert_allclose(covariant_velocity(self.H,self.V,np.zeros_like(self.V)),self.V)

    def test_02_embedded_point_orbitals_modify_analytic_interorbital_velocity(self):
        a=covariant_velocity(self.H,self.V,self.A)
        self.assertGreater(np.max(abs(a-self.V)),.1)
        np.testing.assert_allclose(a,np.swapaxes(a.conj(),-1,-2),atol=1e-12)

    def test_03_nontrivial_orbital_W_derivative_yields_exact_frame_covariance(self):
        o=frame_transform(self.H,self.V,self.A,self.W,self.dW)
        a=covariant_velocity(o['H'],o['derivatives'],o['connection'])
        b=np.stack([self.W.conj().T@v@self.W for v in covariant_velocity(self.H,self.V,self.A)])
        np.testing.assert_allclose(a,b,atol=1e-12)
        self.assertGreater(np.max(abs(o['derivatives']-self.V)),.05)

    def test_04_all_three_geometric_tensors_match_connected_frame(self):
        t=frame_transform(self.H,self.V,self.A,self.W,self.dW)
        a=connected_band_geometry(self.H,self.V,self.A)
        b=connected_band_geometry(t['H'],t['derivatives'],t['connection'])
        for name in ('metric','curvature','moment_energy','metric_pairs','curvature_pairs'):
            np.testing.assert_allclose(getattr(a,name),getattr(b,name),atol=1e-12)
        self.assertEqual(b.acceptance,'VERIFIED_NONDEGENERATE_COMPLETE_FINITE_ORBITAL_CONNECTION')
        wrong=band_geometry(t['H'],t['derivatives'])
        self.assertGreater(np.max(abs(wrong.curvature-a.curvature)),.005)

    def test_05_independent_centered_FD_of_rotated_entire_H(self):
        x,y=self.k
        def h(px,py):
            a=np.zeros((3,3),complex)
            a[0,0]=-1.3+.23*np.cos(px);a[1,1]=.3+.3*np.cos(py)
            a[2,2]=2.1+.2*np.cos(px+py)
            a[0,1]=.43*(np.sin(px)+1j*np.sin(py))+.13j*np.cos(px)
            a[0,2]=.38*(np.sin(py)+1j*np.sin(px))
            a[1,2]=.49*(np.sin(px+py)+1j*np.cos(px))
            for i,j in ((0,1),(0,2),(1,2)):a[j,i]=a[i,j].conj()
            w=expm(1j*(.35*np.sin(px)-.21*np.cos(py))*self.K)
            return w.conj().T@a@w
        o=frame_transform(self.H,self.V,self.A,self.W,self.dW)
        for i in (0,1):
            p=[x,y];q=[x,y];p[i]+=1e-5;q[i]-=1e-5
            np.testing.assert_allclose((h(*p)-h(*q))/2e-5,o['derivatives'][i],atol=3e-10)

    def test_06_B_independent_W_transforms_given_magnetic_operator(self):
        M=np.array([[.08,.2+.1j,-.04j],[.2-.1j,.18,.05],[.04j,.05,-.12]])
        t=frame_transform(self.H,self.V,self.A,self.W,self.dW,magnetic_vertex=M)
        np.testing.assert_allclose(t['magnetic_vertex'],self.W.conj().T@M@self.W,atol=1e-13)

    def test_07_fail_closed_invalid_bloch_embedding(self):
        with self.assertRaises(ValueError):orbital_center_connection([[0,0]])
        with self.assertRaises(ValueError):orbital_center_connection([[0,0],[np.nan,1]])
        with self.assertRaises(ValueError):covariant_velocity(self.H,self.V,self.A[:1])
        with self.assertRaises(ValueError):frame_transform(self.H,self.V,self.A,2*np.eye(3),self.dW)
        with self.assertRaises(ValueError):frame_transform(self.H,self.V,self.A,np.eye(3),self.dW)

    def test_08_symmetric_and_landau_graph_magnetic_spectra_identical(self):
        s=peierls_field(self.H0,self.site,self.B,self.q)
        l=peierls_field(self.H0,self.site,self.B,self.q,gauge='landau_y')
        g,_=field_gauge_unitary(self.site,self.B,self.q)
        np.testing.assert_allclose(l.hamiltonian,g@s.hamiltonian@g.conj().T,atol=1e-14)
        np.testing.assert_allclose(np.linalg.eigvalsh(l.hamiltonian),np.linalg.eigvalsh(s.hamiltonian),atol=1e-14)

    def test_09_Wilson_loop_matches_exact_qBarea_over_hbar(self):
        expected=np.exp(1j*self.q*self.B*1.25*.8)
        for gauge in ('symmetric','landau_y'):
            h=peierls_field(self.H0,self.site,self.B,self.q,gauge=gauge).hamiltonian
            self.assertAlmostEqual(plaquette_wilson(h,[0,1,2,3]),expected,delta=1e-14)
            self.assertAlmostEqual(plaquette_wilson(h,[3,2,1,0]),expected.conjugate(),delta=1e-14)

    def test_10_magnetic_vertex_matches_full_field_difference_in_both_gauges(self):
        step=1e-5
        for gauge in ('symmetric','landau_y'):
            m=peierls_field(self.H0,self.site,self.B,self.q,gauge=gauge).magnetic_vertex
            p=peierls_field(self.H0,self.site,self.B+step,self.q,gauge=gauge).hamiltonian
            n=peierls_field(self.H0,self.site,self.B-step,self.q,gauge=gauge).hamiltonian
            np.testing.assert_allclose(m,(p-n)/(2*step),atol=3e-11)

    def test_11_field_dependent_gauge_B_derivative_contact_is_essential(self):
        s=peierls_field(self.H0,self.site,self.B,self.q)
        l=peierls_field(self.H0,self.site,self.B,self.q,gauge='landau_y')
        g,dg=field_gauge_unitary(self.site,self.B,self.q)
        full=magnetic_vertex_frame_law(s.hamiltonian,s.magnetic_vertex,g,dg)
        np.testing.assert_allclose(full,l.magnetic_vertex,atol=1e-14)
        self.assertGreater(np.max(abs(g@s.magnetic_vertex@g.conj().T-l.magnetic_vertex)),.1)

    def test_12_B_field_reversal_and_real_zero_field_graph(self):
        for gauge in ('symmetric','landau_y'):
            pos=peierls_field(self.H0,self.site,self.B,self.q,gauge=gauge)
            neg=peierls_field(self.H0,self.site,-self.B,self.q,gauge=gauge)
            np.testing.assert_allclose(pos.hamiltonian.conj(),neg.hamiltonian,atol=1e-14)
            np.testing.assert_allclose(-pos.magnetic_vertex.conj(),neg.magnetic_vertex,atol=1e-14)

    def test_13_origin_change_only_changes_link_gauge(self):
        translation=[-.55,1.42]
        for gauge in ('symmetric','landau_y'):
            p=peierls_field(self.H0,self.site,self.B,self.q,gauge=gauge).hamiltonian
            q=peierls_field(self.H0,self.site+translation,self.B,self.q,gauge=gauge).hamiltonian
            np.testing.assert_allclose(np.linalg.eigvalsh(p),np.linalg.eigvalsh(q),atol=5e-14)
            self.assertAlmostEqual(plaquette_wilson(p,[0,1,2,3]),plaquette_wilson(q,[0,1,2,3]),delta=1e-14)

    def test_14_bad_hermitian_path_and_polygon_fail_closed(self):
        with self.assertRaises(ValueError):peierls_field(self.H0,self.site,self.B,self.q,gauge='bad')
        with self.assertRaises(ValueError):peierls_field(self.H0,self.site[:3],self.B,self.q)
        with self.assertRaises(ValueError):peierls_field(self.H0,self.site,self.B,0)
        with self.assertRaises(ValueError):peierls_field(self.H0+1j*np.eye(4),self.site,self.B,self.q)
        with self.assertRaises(ValueError):plaquette_wilson(self.H0,[0,1,2,2])
        with self.assertRaises(ValueError):plaquette_wilson(np.eye(4),[0,1,2,3])

    def test_15_magnetic_operator_not_fake_total_magnetotransport(self):
        out=peierls_field(self.H0,self.site,self.B,self.q)
        self.assertEqual(out.acceptance,'VERIFIED_FINITE_GRAPH_STRAIGHT_PEIERLS_MAGNETIC_VERTEX')
        self.assertFalse(hasattr(out,'orbital_mnhe'))

if __name__=='__main__':unittest.main()
