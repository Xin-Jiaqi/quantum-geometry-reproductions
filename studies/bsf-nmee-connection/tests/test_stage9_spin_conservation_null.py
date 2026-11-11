"""Stage9 falsification: Ising Sz conservation kills transverse static NMEE."""
import unittest
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingBilayer,SlidingSpec
from bsf_nmee.valley_kp import static_kernel

class SpinConservationNoGoTests(unittest.TestCase):
    def test_spin_conserving_limit_kills_even_and_odd_spin_response(self):
        spec=SlidingSpec(shift=(.23*.4,.41*.4),
                         tunneling_eV=.005,rashba_soc_eV=0.)
        obj=SlidingBilayer(spec)
        centers=np.array([[1/3,2/3],[2/3,1/3]])@obj.model().reciprocal
        offsets=np.array([[.07,.01],[.13,.08],[-.04,.10]])
        momenta=np.concatenate([center+offsets for center in centers])
        h=obj.matrices(momenta)
        v=np.stack([obj.matrices(momenta,axis=j) for j in (0,1)])
        sz=obj.spins["total_z"]
        sx_even=obj.spins["total_x"]
        sx_odd=obj.spins["top_x"]-obj.spins["bottom_x"]
        for operator in [*h,*v[0],*v[1]]:
            self.assertLess(np.max(abs(operator@sz-sz@operator)),1e-11)
        self.assertLess(np.max(abs(sx_even@sz+sz@sx_even)),1e-12)
        self.assertLess(np.max(abs(sx_odd@sz+sz@sx_odd)),1e-12)
        response=static_kernel(h,v,[sx_even,sx_odd],
             sz,spec.mu_eV,g_factor=obj.mono.g_factor,step_T=.25)
        self.assertLess(float(np.max(abs(response))),1e-10)

    def test_full_spin_mixing_is_not_symmetry_forbidden_in_model(self):
        spec=SlidingSpec(shift=(.23*.4,.41*.4),
                         tunneling_eV=.005,rashba_soc_eV=.020)
        obj=SlidingBilayer(spec)
        momenta=(np.array([[1/3,2/3]])@obj.model().reciprocal)+np.array([.08,.04])
        h=obj.matrices(momenta)
        sz=obj.spins["total_z"]
        self.assertGreater(float(np.max(abs(h[0]@sz-sz@h[0]))),1e-4)
        # This verifies only the required nonconservation gate, NOT
        # finite integrated response, which depends on other physics.

if __name__=="__main__":
    unittest.main()
