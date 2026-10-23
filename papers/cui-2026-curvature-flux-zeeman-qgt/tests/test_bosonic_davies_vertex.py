"""Round31 finite Fock/Ohmic-boson spectral weak-coupling and dynamic Ward gates."""
import unittest
import numpy as np
from flux_lab.fermionic_fock_bath import (
    fock_soc_hubbard,grand_gibbs,unitary_phase,local_charge_balance,
    thermal_sector_jumps)
from flux_lab.thermal_lindblad_ward import dissipator
from flux_lab.bosonic_davies_vertex import (
    bose_occupation,ohmic_spectral_density,oscillator_fourier_rate,
    finite_oscillator_correlator,bose_davies_fock_jumps,
    bath_magnetic_derivative,linear_response,bath_contact)

BETA=.83
G2=.045
CUTOFF=2.4
MU=.16
def sample():
    d=fock_soc_hubbard()
    b=bose_davies_fock_jumps(d['H'],d['N'],d['coupling'],beta=BETA,g2=G2,cutoff=CUTOFF)
    rho=grand_gibbs(d['H'],d['N'],beta=BETA,mu=MU)
    return d,b,rho

class OhmicBosonicFockWard(unittest.TestCase):
    def test_bose_occupations_and_ohmic_cutoff_positive(self):
        self.assertGreater(bose_occupation(.7,BETA),0)
        self.assertGreater(ohmic_spectral_density(1.2,g2=G2,cutoff=CUTOFF),0)
        self.assertEqual(ohmic_spectral_density(0,g2=G2,cutoff=CUTOFF),0)
        self.assertGreater(oscillator_fourier_rate(-1.2),oscillator_fourier_rate(1.2))

    def test_derived_bath_spectrum_kms_for_multiple_positive_frequencies(self):
        for w in (.04,.27,1.2,4.5):
            absorption=oscillator_fourier_rate(w,beta=BETA,g2=G2,cutoff=CUTOFF)
            emission=oscillator_fourier_rate(-w,beta=BETA,g2=G2,cutoff=CUTOFF)
            self.assertAlmostEqual(absorption/emission,np.exp(-BETA*w),delta=2e-14)
        self.assertAlmostEqual(oscillator_fourier_rate(0,beta=BETA,g2=G2),
                               2*np.pi*G2/BETA,places=13)

    def test_independent_finite_harmonic_oscillator_bath_trace(self):
        w=1.25
        d=finite_oscillator_correlator(w,beta=BETA,levels=48)
        self.assertAlmostEqual(d["absorption"]/d["emission"],np.exp(-BETA*w),delta=3e-14)
        self.assertAlmostEqual(d["emission"],d["occupation"]+1,delta=5e-14)
        self.assertAlmostEqual(d["absorption"],d["occupation"],delta=5e-14)
        self.assertAlmostEqual(np.trace(d["rho"]).real,1,places=14)

    def test_exact_oscillator_correlation_matches_thermal_matrix_trace(self):
        w=1.25;t=.37;d=finite_oscillator_correlator(w,beta=BETA,levels=48)
        E=np.diag(d["H"]).real
        U=np.diag(np.exp(1j*t*E))
        B=d["b"]+d["b"].conj().T
        Bt=U@B@U.conj().T
        direct=np.trace(d["rho"]@Bt@B)
        spectral=d["emission"]*np.exp(-1j*w*t)+d["absorption"]*np.exp(1j*w*t)
        np.testing.assert_allclose(direct,spectral,atol=4e-14)

    def test_interacting_fermion_sectors_kms_from_bose_spectral_function(self):
        d,b,rho=sample()
        self.assertEqual([int(np.sum(b["numbers"]==i)) for i in range(5)],[1,4,6,4,1])
        self.assertGreater(b["transition_count"],20)
        self.assertEqual(len(b["jumps"]),b["transition_count"]+1)
        p=np.exp(-BETA*(b["energies"]-MU*b["numbers"]))
        p/=p.sum()
        np.testing.assert_allclose(b["rates"]*p[None,:],
                                   b["rates"].T*p[:,None],atol=4e-14)
        self.assertEqual(np.count_nonzero(b["rates"][
            b["numbers"][:,None]!=b["numbers"][None,:]]),0)

    def test_grand_gibbs_stationary_including_zero_frequency_bath(self):
        d,b,rho=sample()
        np.testing.assert_allclose(dissipator(rho,b["jumps"]),0,atol=3e-13)
        self.assertGreater(np.linalg.norm(b["zero_frequency_jump"]),.02)
        np.testing.assert_allclose(b["zero_frequency_jump"],
                                   b["zero_frequency_jump"].conj().T,atol=2e-14)

    def test_not_merely_round30_exponential_postulated_rates(self):
        d,b,_=sample()
        old=thermal_sector_jumps(d["H"],d["N"],d["coupling"],beta=BETA,gamma=.14)
        nz=(b["rates"]>.000001)&(old["rates"]>.000001)
        self.assertGreater(np.count_nonzero(nz),10)
        ratio=b["rates"][nz]/old["rates"][nz]
        self.assertGreater(np.max(ratio)-np.min(ratio),.05)

    def test_total_electron_number_conserved_for_general_fock_density(self):
        d,b,_=sample()
        rng=np.random.default_rng(31)
        Q=rng.normal(size=(16,16))+1j*rng.normal(size=(16,16))
        X=Q@Q.conj().T;X/=np.trace(X)
        source=dissipator(X,b["jumps"])
        for n in (d["N"],d["N"]@d["N"]):
            self.assertLess(abs(np.trace(n@source)),3e-13)
        for J in b["jumps"]:
            np.testing.assert_allclose(J@d["N"],d["N"]@J,atol=5e-13)

    def test_local_bath_source_required_in_exact_fock_continuity(self):
        d,b,_=sample()
        X=np.zeros((16,16),complex);X[3,3]=1
        result=local_charge_balance(d,X,b["jumps"])
        self.assertGreater(np.max(abs(result["bath_source"])),.005)
        np.testing.assert_allclose(result["continuity"],0,atol=3e-13)
        self.assertLess(abs(result["total_source"]),3e-13)

    def test_joint_gauge_covar_stale_jump_failure_and_contact(self):
        d,b,rho=sample();G=d["gauge"];phi=.39
        R=lambda X:unitary_phase(X,G,phi)
        js=np.array([R(j) for j in b["jumps"]])
        X=np.zeros((16,16),complex);X[3,3]=.7;X+=.3*rho
        np.testing.assert_allclose(dissipator(R(X),js),
                                   R(dissipator(X,b["jumps"])),atol=4e-13)
        self.assertGreater(np.max(abs(dissipator(R(X),b["jumps"])-
                                  R(dissipator(X,b["jumps"])))),.0001)
        gauge_dr=1j*(G@rho-rho@G)
        full=dissipator(gauge_dr,b["jumps"])+bath_contact(rho,b["jumps"],G)
        np.testing.assert_allclose(full,0,atol=5e-13)
        self.assertGreater(np.max(abs(dissipator(gauge_dr,b["jumps"]))),.0001)

    def test_gauge_bath_vertex_matches_independent_finite_amplitude(self):
        d,b,rho=sample();G=d["gauge"];eps=2e-5
        def at(s):
            jumps=np.array([unitary_phase(j,G,s*eps) for j in b["jumps"]])
            return dissipator(rho,jumps)
        numeric=(at(1)-at(-1))/(2*eps)
        np.testing.assert_allclose(numeric,bath_contact(rho,b["jumps"],G),
                                   atol=2e-9)

    def test_zeeman_bath_vertex_includes_bosonic_energy_dependence(self):
        d,b,rho=sample()
        x=bath_magnetic_derivative(d,beta=BETA,g2=G2,cutoff=CUTOFF,mu=MU)
        self.assertGreater(np.max(abs(x["equilibrium_collision_state"])),.0001)
        np.testing.assert_allclose(x["residual"],0,atol=2e-8)
        self.assertGreater(np.max(abs(x["collision_contact"])),.0001)

    def test_finite_frequency_pure_gauge_ward_with_bath_vertex(self):
        d,b,rho=sample()
        r=linear_response(d,mode="gauge",beta=BETA,mu=MU,g2=G2,cutoff=CUTOFF)
        exact=1j*(d["gauge"]@rho-rho@d["gauge"])
        np.testing.assert_allclose(r["dr"],exact,atol=2e-11)
        np.testing.assert_allclose(r["density"],0,atol=2e-11)
        self.assertLess(abs(r["current"]),2e-11)
        np.testing.assert_allclose(r["ward"],0,atol=2e-11)

    def test_omitted_dynamic_bath_vertex_fails_pure_gauge(self):
        d,b,rho=sample()
        good=linear_response(d,mode="gauge",beta=BETA,mu=MU,g2=G2,cutoff=CUTOFF)
        bad=linear_response(d,mode="gauge",beta=BETA,mu=MU,g2=G2,cutoff=CUTOFF,
                            include_bath_contact=False)
        self.assertGreater(np.max(abs(good["dr"]-bad["dr"])),1e-5)
        gauge_dr=1j*(d["gauge"]@rho-rho@d["gauge"])
        self.assertGreater(np.max(abs(bad["dr"]-gauge_dr)),1e-5)

    def test_dynamic_physical_zeeman_bath_response(self):
        d,b,_=sample()
        x=linear_response(d,mode="magnetic",beta=BETA,mu=MU,g2=G2,cutoff=CUTOFF)
        self.assertGreater(np.max(abs(x["dr"])),1e-3)
        np.testing.assert_allclose(x["ward"],0,atol=3e-10)
        bad=linear_response(d,mode="magnetic",beta=BETA,mu=MU,g2=G2,cutoff=CUTOFF,
                            include_bath_contact=False)
        self.assertGreater(np.max(abs(x["dr"]-bad["dr"])),1e-5)

    def test_fail_closed_beta_frequency_bath_inputs(self):
        d,_,_=sample()
        with self.assertRaises(ValueError):bose_occupation(0,BETA)
        with self.assertRaises(ValueError):ohmic_spectral_density(-.1)
        with self.assertRaises(ValueError):oscillator_fourier_rate(1,beta=0)
        with self.assertRaises(ValueError):oscillator_fourier_rate(0,g2=-1)
        with self.assertRaises(ValueError):finite_oscillator_correlator(1,levels=1)
        with self.assertRaises(ValueError):bose_davies_fock_jumps(d["H"],d["N"],d["coupling"],beta=-1)
        with self.assertRaises(ValueError):bose_davies_fock_jumps(d["H"],d["N"],d["coupling"],g2=0)
        with self.assertRaises(ValueError):linear_response(d,mode="gauge",z=.4)
        with self.assertRaises(ValueError):linear_response(d,mode="bad")

if __name__=="__main__":unittest.main()
