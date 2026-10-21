"""v58 fixed-density planar coherent Hall bars: independent device/gauge/size checks."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.planar_hall_scaling_v58 import (
    PlanarSpec,planar_bar,gamma_profiles,caroli_multisite,planar_response,
    full_B_derivative,compare_planar_family,adaptive_planar_response)
from qgmap.lead_transport_v57 import LeadSpec,_transmission_at_energy
from qgmap.embedding_magnetic_v55 import peierls_field

class PlanarHallV58(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lead=LeadSpec(mu=.65,temperature=.12,field=.12,charge_over_hbar=-.87,energy_points=80)
        cls.p5=PlanarSpec(nx=5,ny=5,coverage=1.,gamma_per_site=.6)
        cls.p75=replace(cls.p5,nx=7)
        cls.p77=replace(cls.p5,nx=7,ny=7)
        cls.a=planar_response(cls.p5,cls.lead)
        cls.b=planar_response(cls.p75,cls.lead)

    def test_01_fixed_density_local_hopping_onsite_and_edge_count(self):
        for p in (self.p5,self.p75,self.p77):
            positions,H,sites=planar_bar(p)
            self.assertEqual(len(positions),p.nx*p.ny)
            self.assertEqual(sum(len(s) for s in sites),2*(p.nx-2)+2*(p.ny-2))
            self.assertAlmostEqual(self.a['lattice_density'],1.)
            for i,j in np.argwhere(abs(H)>1e-12):
                if i!=j:
                    self.assertAlmostEqual(abs(H[i,j]),1.)
                    self.assertAlmostEqual(np.linalg.norm(positions[i]-positions[j]),1.)
            self.assertTrue(all(abs(H[i,i]-p.onsite)<1e-12 for i in range(len(H))))

    def test_02_identical_flux_per_plaquette_at_different_sample_size(self):
        for p in (self.p5,self.p75,self.p77):
            a=planar_response(p,self.lead)
            self.assertAlmostEqual(a['plaquette_flux'],-.1044,delta=1e-13)
            self.assertAlmostEqual(a['lattice_density'],1.)

    def test_03_disjoint_symmetric_contacts_on_opposite_edges(self):
        for cov in (.2,.4,.6,.8,1.):
            for p in (replace(self.p5,coverage=cov),replace(self.p77,coverage=cov)):
                positions,H,leads=planar_bar(p)
                used=[i for edge in leads for i in edge]
                self.assertEqual(len(used),len(set(used)))
                self.assertEqual(len(leads[0]),len(leads[2]))
                self.assertEqual(len(leads[1]),len(leads[3]))
                for alpha,edge in enumerate(leads):
                    coords=positions[list(edge)]
                    axis=1 if alpha in (0,2) else 0
                    self.assertAlmostEqual(coords[:,axis].mean(),
                       (p.ny-1 if axis==1 else p.nx-1)/2,delta=1e-12)

    def test_04_multisite_trace_matches_independent_Caroli_matrices(self):
        pos,H,leadsites=planar_bar(self.p5)
        gam=gamma_profiles(leadsites,len(H),.6)
        pb=peierls_field(H,pos,.12,-.87)
        E=.31
        T,dT=caroli_multisite(E,pb.hamiltonian,pb.magnetic_vertex,gam)
        G=np.linalg.inv(E*np.eye(len(H))-pb.hamiltonian+.5j*np.diag(gam.sum(axis=0)))
        dG=G@pb.magnetic_vertex@G
        for a in range(4):
            for b in range(4):
                if a==b:continue
                A=np.diag(gam[a]);B=np.diag(gam[b]);Ga=G.conj().T
                exact=np.trace(A@G@B@Ga).real
                derivative=2*np.real(np.trace(A@dG@B@Ga))
                self.assertAlmostEqual(T[a,b],exact,delta=2e-13)
                self.assertAlmostEqual(dT[a,b],derivative,delta=2e-12)

    def test_05_single_site_contact_reduces_to_independent_v57(self):
        p=replace(self.p5,coverage=.4)
        pos,H,leads=planar_bar(p)
        self.assertEqual(tuple(map(len,leads)),(1,1,1,1))
        gam=gamma_profiles(leads,len(H),.6)
        pb=peierls_field(H,pos,self.lead.field,self.lead.charge_over_hbar)
        T,dT=caroli_multisite(.27,pb.hamiltonian,pb.magnetic_vertex,gam)
        single,dsingle=_transmission_at_energy(.27,pb.hamiltonian,pb.magnetic_vertex,
                         [lead[0] for lead in leads],[.6]*4)
        np.testing.assert_allclose(T,single,atol=5e-13)
        np.testing.assert_allclose(dT,dsingle,atol=5e-12)

    def test_06_multi_contact_T_and_derivative_B_gauge_covariance(self):
        for p in (self.p5,self.p75,self.p77):
            s=planar_response(p,self.lead,gauge='landau_y')
            o=planar_response(p,self.lead)
            np.testing.assert_allclose(s['T'],o['T'],atol=3e-12)
            np.testing.assert_allclose(s['dT_dB'],o['dT_dB'],atol=3e-11)
            self.assertAlmostEqual(s['observables']['hall_ratio'],o['observables']['hall_ratio'],delta=3e-12)

    def test_07_optical_theorem_and_terminal_current_conservation(self):
        for p in (self.p5,self.p75,self.p77):
            r=planar_response(p,self.lead)
            np.testing.assert_allclose(r['T'].sum(axis=1),r['T'].sum(axis=0),atol=2e-11)
            self.assertLess(r['observables']['max_floating_probe_current'],1e-11)
            self.assertLess(r['observables']['max_current_conservation_residual'],1e-11)

    def test_08_onsager_casimir_for_all_contacts(self):
        plus=self.a
        minus=planar_response(self.p5,replace(self.lead,field=-self.lead.field))
        np.testing.assert_allclose(plus['T'],minus['T'].T,atol=5e-12)
        np.testing.assert_allclose(plus['dT_dB'],-minus['dT_dB'].T,atol=5e-11)

    def test_09_symmetric_rectangles_have_odd_Hall_and_zero_B0(self):
        for p in (self.p5,self.p75,self.p77):
            plus=planar_response(p,self.lead)
            minus=planar_response(p,replace(self.lead,field=-self.lead.field))
            zero=planar_response(p,replace(self.lead,field=0))
            self.assertAlmostEqual(plus['observables']['hall_ratio'],
                                   -minus['observables']['hall_ratio'],delta=2e-11)
            self.assertLess(abs(zero['observables']['hall_ratio']),2e-11)
            self.assertAlmostEqual(plus['observables']['hall_resistance_h_over_q2'],
                                  -minus['observables']['hall_resistance_h_over_q2'],delta=3e-11)

    def test_10_true_current_normalized_Hall_resistance_and_analytic_slope(self):
        obs=self.a['observables']
        self.assertAlmostEqual(obs['hall_resistance_h_over_q2'],
            obs['hall_ratio']/obs['longitudinal_conductance_q2_over_h'],delta=1e-14)
        p=obs['longitudinal_conductance_q2_over_h']
        target=(obs['dHall_dB']*p-obs['hall_ratio']*obs['dLongitudinal_dB'])/p**2
        self.assertAlmostEqual(obs['dHall_resistance_dB'],target,delta=1e-13)

    def test_11_analytic_field_slope_vs_whole_contact_voltage_fd(self):
        for p in (self.p5,self.p75):
            a=planar_response(p,self.lead)
            fd=full_B_derivative(p,self.lead,step=1e-4)
            np.testing.assert_allclose(a['dT_dB'],fd['dT_dB'],atol=5e-8)
            self.assertAlmostEqual(a['observables']['dHall_dB'],fd['dHall_dB'],delta=2e-8)
            self.assertAlmostEqual(a['observables']['dHall_resistance_dB'],
                                   fd['dHall_resistance_dB'],delta=3e-8)

    def test_12_fixed_total_reservoir_gamma_vs_per_site_gamma(self):
        p=replace(self.p77,coverage=1.)
        pos,h,leads=planar_bar(p)
        per=gamma_profiles(leads,len(h),.6,'per_site')
        total=gamma_profiles(leads,len(h),.6,'fixed_reservoir_sum')
        np.testing.assert_allclose(per.sum(axis=1),[3.,3.,3.,3.])
        np.testing.assert_allclose(total.sum(axis=1),[.6]*4)
        a=planar_response(p,self.lead)
        b=planar_response(replace(p,gamma_mode='fixed_reservoir_sum'),self.lead)
        self.assertGreater(abs(a['observables']['hall_ratio']-b['observables']['hall_ratio']),.01)

    def test_13_coverage_changes_contact_width_and_Hall_even_at_fixed_size(self):
        p=replace(self.p77,coverage=.4)
        a=planar_response(p,self.lead)
        b=planar_response(self.p77,self.lead)
        self.assertEqual(a['sample_sites'],b['sample_sites'])
        self.assertEqual(a['contact_site_counts'],(1,1,1,1))
        self.assertEqual(b['contact_site_counts'],(5,5,5,5))
        self.assertGreater(abs(a['observables']['hall_resistance_h_over_q2']-
                               b['observables']['hall_resistance_h_over_q2']),.1)

    def test_14_width_and_length_resize_keep_density_but_do_not_force_convergence(self):
        result=compare_planar_family([self.p5,self.p75,self.p77],self.lead)
        self.assertEqual(len(result['rows']),3)
        self.assertEqual(len({r['density'] for r in result['rows']}),1)
        self.assertEqual(len({r['plaquette_flux'] for r in result['rows']}),1)
        self.assertFalse(result['thermodynamic_claim'])
        assert len({round(r['dHall_dB'],4) for r in result['rows']})==3
        self.assertLess(result['rows'][0]['sites'],result['rows'][1]['sites'])

    def test_15_rotation_of_rectangular_sample_preserves_current_normalized_Hall(self):
        a=planar_response(replace(self.p5,nx=7,ny=5),self.lead)
        b=planar_response(replace(self.p5,nx=5,ny=7),self.lead)
        self.assertAlmostEqual(a['observables']['hall_resistance_h_over_q2'],
                               b['observables']['hall_resistance_h_over_q2'],delta=3e-12)
        self.assertAlmostEqual(a['observables']['dHall_resistance_dB'],
                               b['observables']['dHall_resistance_dB'],delta=2e-11)

    def test_16_contact_variation_separate_from_density_change(self):
        p1=replace(self.p77,coverage=.6)
        p2=replace(self.p77,coverage=1.)
        result=compare_planar_family([p1,p2],self.lead)
        self.assertEqual(result['rows'][0]['sites'],result['rows'][1]['sites'])
        self.assertEqual(result['rows'][0]['density'],result['rows'][1]['density'])
        self.assertGreater(abs(result['rows'][0]['hall_ratio']-result['rows'][1]['hall_ratio']),.005)

    def test_17_detect_gaussian_underconvergence_then_adaptive_reference(self):
        p=self.p5
        a=planar_response(p,replace(self.lead,energy_points=120))
        b=planar_response(p,replace(self.lead,energy_points=180))
        # A too-small finite quadrature MUST NOT be used as an accepted source.
        self.assertGreater(abs(a['observables']['hall_ratio']-b['observables']['hall_ratio']),5e-5)
        reference=adaptive_planar_response(p,self.lead)
        strict=adaptive_planar_response(p,self.lead,epsrel=7e-9,epsabs=3e-11)
        self.assertTrue(reference['quadrature_success'])
        self.assertAlmostEqual(reference['observables']['hall_ratio'],strict['observables']['hall_ratio'],delta=2e-8)
        self.assertAlmostEqual(reference['observables']['dHall_dB'],strict['observables']['dHall_dB'],delta=2e-7)
        self.assertAlmostEqual(reference['observables']['hall_ratio'],.014841643736028873,delta=3e-8)

    def test_18_same_lattice_geometry_reproduces_change_of_gauge_and_field_scale(self):
        p=replace(self.p5,lattice_spacing=2.)
        b=replace(self.lead,field=self.lead.field/4)
        # Peierls phase depends on q B area and total spectrum on fixed hopping.
        a=planar_response(self.p5,self.lead)
        c=planar_response(p,b)
        np.testing.assert_allclose(a['T'],c['T'],atol=2e-12)
        self.assertAlmostEqual(c['observables']['dHall_dB'],4*a['observables']['dHall_dB'],delta=5e-11)

    def test_19_reject_overlap_invalid_density_and_invalid_gamma_protocol(self):
        for vals in ({'nx':2},{'ny':2},{'coverage':0},{'coverage':1.01},
                     {'gamma_per_site':0},{'gamma_mode':'unknown'},{'lattice_spacing':0}):
            with self.assertRaises(ValueError):PlanarSpec(**vals)
        with self.assertRaises(ValueError):gamma_profiles(((1,),(1,),(2,),(3,)),10,.6)
        with self.assertRaises(ValueError):gamma_profiles(((1,),(),(2,),(3,)),10,.6)
        with self.assertRaises(ValueError):caroli_multisite(.2,np.eye(4),np.eye(4),np.ones((3,4)))
        with self.assertRaises(ValueError):full_B_derivative(self.p5,self.lead,step=-1)

    def test_21_adaptive_two_magnetic_gauges_resolve_tiny_resonances(self):
        a=adaptive_planar_response(self.p5,self.lead)
        b=adaptive_planar_response(self.p5,self.lead,gauge='landau_y')
        np.testing.assert_allclose(a['T'],b['T'],atol=3e-10)
        np.testing.assert_allclose(a['dT_dB'],b['dT_dB'],atol=3e-9)
        self.assertAlmostEqual(a['observables']['dHall_dB'],b['observables']['dHall_dB'],delta=2e-8)

    def test_22_adaptive_current_contacts_and_true_Hall_resistance(self):
        a=adaptive_planar_response(self.p75,self.lead)
        self.assertLess(a['observables']['max_floating_probe_current'],3e-11)
        self.assertLess(a['observables']['max_current_conservation_residual'],3e-11)
        self.assertAlmostEqual(a['observables']['hall_resistance_h_over_q2'],
          a['observables']['hall_ratio']/a['observables']['longitudinal_conductance_q2_over_h'],delta=1e-13)

    def test_23_adaptive_full_field_derivative_uses_full_voltage_probe_resolve(self):
        s=self.lead;p=self.p75;step=1e-4
        source=adaptive_planar_response(p,s)
        parts=[adaptive_planar_response(p,replace(s,field=s.field+d*step)) for d in (2,1,-1,-2)]
        phi=lambda k: (-parts[0]['observables'][k]+8*parts[1]['observables'][k]
                        -8*parts[2]['observables'][k]+parts[3]['observables'][k])/(12*step)
        self.assertAlmostEqual(source['observables']['dHall_dB'],phi('hall_ratio'),delta=4e-7)
        self.assertAlmostEqual(source['observables']['dHall_resistance_dB'],
                               phi('hall_resistance_h_over_q2'),delta=6e-7)

    def test_24_adaptive_contact_coverage_is_not_an_energy_mesh_artifact(self):
        full=adaptive_planar_response(self.p5,self.lead)
        narrow=adaptive_planar_response(replace(self.p5,coverage=.4),self.lead)
        self.assertGreater(abs(narrow['observables']['hall_ratio']-full['observables']['hall_ratio']),.1)
        self.assertEqual(narrow['contact_site_counts'],(1,1,1,1))
        self.assertNotEqual(narrow['integrand_evaluations'],full['integrand_evaluations'])

    def test_25_adaptive_rejects_unverified_numerical_tolerances(self):
        for e in (0,-1,float('nan')):
            with self.assertRaises(ValueError):adaptive_planar_response(self.p5,self.lead,epsrel=e)

    def test_20_scope_has_no_bulk_material_tensor_or_fake_size_convergence(self):
        a=compare_planar_family([self.p5,self.p77],self.lead)
        self.assertIn('NOT_BULK',self.a['status'])
        self.assertFalse(a['thermodynamic_claim'])
        self.assertTrue(any('not a bulk conductivity' in s for s in a['caveats']))

if __name__=='__main__':unittest.main()
