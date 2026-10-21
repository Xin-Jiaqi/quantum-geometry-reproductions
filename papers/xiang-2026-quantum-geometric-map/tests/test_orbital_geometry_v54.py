"""v54 VERIFIED core: independent projector physics, multiband covariance and decomposition."""
import unittest
import numpy as np
from qgmap.orbital_geometry_v54 import (
 band_geometry,curvature_moment_product,moment_curvature_kernel,
 periodic_three_band_model,finite_difference_projector_reference)

sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex)
sz=np.diag([1.,-1.]).astype(complex)

def two(x,y,m=1.):
    return np.sin(x)*sx+np.sin(y)*sy+(m+.5*(np.cos(x)+np.cos(y)))*sz

def vtwo(x,y):
    return np.array([np.cos(x)*sx-.5*np.sin(x)*sz,np.cos(y)*sy-.5*np.sin(y)*sz])

def three(x,y):
    H=np.zeros((3,3),complex)
    H[0,0]=-1.3+.23*np.cos(x);H[1,1]=.3+.3*np.cos(y)
    H[2,2]=2.1+.2*np.cos(x+y)
    H[0,1]=.43*(np.sin(x)+1j*np.sin(y))+.13j*np.cos(x)
    H[0,2]=.38*(np.sin(y)+1j*np.sin(x))
    H[1,2]=.49*(np.sin(x+y)+1j*np.cos(x))
    for i,j in ((0,1),(0,2),(1,2)):H[j,i]=H[i,j].conj()
    return H

class OrbitalGeometryV54(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.k,cls.H,cls.V=periodic_three_band_model(30)
        cls.geo=band_geometry(cls.H,cls.V)

    def test_two_band_analytic_velocities_vs_independent_projector(self):
        for x,y in ((.29,.48),(-.42,.27),(.8,-1.1)):
            a=band_geometry(two(x,y),vtwo(x,y))
            b=finite_difference_projector_reference(two,x,y,8e-5)
            np.testing.assert_allclose(a.curvature,b['curvature'],atol=2e-8)
            np.testing.assert_allclose(a.metric[...,0,1],b['metric_xy'],atol=2e-8)
            np.testing.assert_allclose(a.moment_energy,b['moment_energy'],atol=3e-8)

    def test_three_band_projector_is_independent_oracle(self):
        k=self.k
        a=band_geometry(self.H[4,8],self.V[:,4,8])
        b=finite_difference_projector_reference(three,float(k[4]),float(k[8]),8e-5)
        np.testing.assert_allclose(a.curvature,b['curvature'],atol=2e-8)
        np.testing.assert_allclose(a.metric[...,0,1],b['metric_xy'],atol=2e-8)
        np.testing.assert_allclose(a.moment_energy,b['moment_energy'],atol=3e-8)

    def test_band_pair_curvature_trace_and_metric_positivity(self):
        g=self.geo
        np.testing.assert_allclose(g.curvature.sum(axis=-1),0,atol=3e-14)
        np.testing.assert_allclose(g.curvature_pairs,-np.swapaxes(g.curvature_pairs,-1,-2),atol=3e-14)
        np.testing.assert_allclose(g.metric_pairs,
            np.swapaxes(g.metric_pairs,-4,-3).swapaxes(-2,-1),atol=3e-14)
        self.assertGreater(np.min(g.metric[...,0,0]),-1e-12)
        self.assertGreater(np.min(g.metric[...,1,1]),-1e-12)

    def test_arbitrary_independent_band_U1_phases_at_each_k(self):
        E,U=np.linalg.eigh(self.H)
        U2=U*np.exp(1j*np.random.default_rng(54).normal(size=E.shape))[...,None,:]
        b=band_geometry(self.H,self.V,eigenvectors=U2)
        for name in ('metric','curvature','moment_energy','curvature_pairs'):
            np.testing.assert_allclose(getattr(self.geo,name),getattr(b,name),atol=2e-12)

    def test_random_constant_orbital_basis_unitary(self):
        rng=np.random.default_rng(23)
        A=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        Q,_=np.linalg.qr(A)
        b=band_geometry(Q@self.H@Q.conj().T,Q@self.V@Q.conj().T)
        for name in ('metric','curvature','moment_energy'):
            np.testing.assert_allclose(getattr(self.geo,name),getattr(b,name),atol=2e-12)

    def test_energy_shift_and_uniform_rescaling(self):
        b=band_geometry(self.H+12*np.eye(3),self.V)
        np.testing.assert_allclose(self.geo.curvature,b.curvature,atol=2e-12)
        np.testing.assert_allclose(self.geo.moment_energy,b.moment_energy,atol=2e-12)
        a=3.2
        c=band_geometry(self.H*a,self.V*a)
        np.testing.assert_allclose(self.geo.metric,c.metric,atol=2e-12)
        np.testing.assert_allclose(self.geo.curvature,c.curvature,atol=2e-12)
        np.testing.assert_allclose(self.geo.moment_energy*a,c.moment_energy,atol=2e-12)

    def test_two_band_moment_equal_upper_lower(self):
        g=band_geometry(two(.3,.5),vtwo(.3,.5))
        self.assertAlmostEqual(g.moment_energy[0],g.moment_energy[1],delta=2e-14)
        self.assertAlmostEqual(g.curvature[0],-g.curvature[1],delta=2e-14)

    def test_three_band_explicit_m_not_equal_l_double_sum(self):
        g=self.geo
        out=curvature_moment_product(g)
        crossed=np.zeros_like(g.curvature)
        for n in range(3):
            for m in range(3):
                for l in range(3):
                    if m!=l:crossed[...,n]+=g.curvature_pairs[...,n,m]*g.moment_pairs[...,n,l]
        np.testing.assert_allclose(crossed,out['cross_m_ne_l'],atol=4e-15)
        np.testing.assert_allclose(out['full'],out['diagonal_m_eq_l']+crossed,atol=4e-15)

    def test_independent_v53_three_band_geometric_subsector(self):
        _,H,V=periodic_three_band_model(48)
        out=moment_curvature_kernel(band_geometry(H,V),-1.1,.45)
        self.assertAlmostEqual(out['full'],-.000799451053562,delta=7e-12)
        self.assertAlmostEqual(out['diagonal_m_eq_l'],-.001542183149815,delta=8e-12)
        self.assertAlmostEqual(out['cross_m_ne_l'],+.000742732096253,delta=8e-12)

    def test_fermi_window_changes_three_band_fraction(self):
        a=moment_curvature_kernel(self.geo,-1.1,.45)
        b=moment_curvature_kernel(self.geo,1.5,.45)
        self.assertGreater(abs(a['cross_m_ne_l']/a['full']),.9)
        self.assertLess(abs(b['cross_m_ne_l']/b['full']),.07)

    def test_analytically_derived_velocity_matches_finite_k_H(self):
        k,H,V=periodic_three_band_model(8)
        x,y=float(k[1]),float(k[4]);step=1e-5
        for i in (0,1):
            p=[x,y];q=[x,y];p[i]+=step;q[i]-=step
            np.testing.assert_allclose(V[i,1,4],(three(*p)-three(*q))/(2*step),atol=2e-10)

    def test_reject_degenerate_nonsquare_nonhermitian_and_bad_eigenvectors(self):
        with self.assertRaisesRegex(ValueError,'degeneracy'):
            band_geometry(np.zeros((3,3)),np.zeros((2,3,3)))
        with self.assertRaisesRegex(ValueError,'velocities shape'):
            band_geometry(np.diag([0.,1.,2.]),np.zeros((3,3)))
        bad=np.diag([0.,1.,2.]).astype(complex);bad[0,1]=1
        with self.assertRaisesRegex(ValueError,'Hermitian'):
            band_geometry(bad,np.zeros((2,3,3)))
        with self.assertRaisesRegex(ValueError,'not matched'):
            band_geometry(self.H[0,0],self.V[:,0,0],eigenvectors=np.eye(3))

    def test_verified_geometry_but_not_claimed_total_response(self):
        g=self.geo
        parts=curvature_moment_product(g)
        self.assertEqual(g.acceptance,'VERIFIED_NONDEGENERATE_FIXED_ORBITAL_BASIS')
        self.assertEqual(parts['acceptance'],'VERIFIED_PRODUCT_DECOMPOSITION_NOT_COMPLETE_RESPONSE')
        self.assertNotIn('conductivity',parts)

if __name__=='__main__':unittest.main()
