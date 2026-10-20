"""Independent field-transformation checks for the Stage-5 physics interpretation.

These checks derive inter-domain coefficients from field transformation, not
from the original `response_sign` code. No numerical AgF2 IGMC is claimed.
"""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from agf2_preflight import allowed_components, partner_sign, ahe_z_sign, R_MZ

R = (1, 1, -1)  # polar vector action for mirror z

def induced_polarization_sign(i, j, antiunitary=False):
    # P polar / T even. B axial / T odd => field B_j has factor det(R)*R_j*(-1)^theta.
    parity_p = R[i]
    parity_b = (-1) * R[j] * (-1 if antiunitary else 1)
    return parity_p * parity_b

class TestAdiabaticPhysics(unittest.TestCase):
    def test_in_plane_longitudinal(self):
        for index,label in [(0,'xx'),(1,'yy')]:
            self.assertEqual(induced_polarization_sign(index,index), -1)
            self.assertEqual(induced_polarization_sign(index,index,True), +1)
            self.assertEqual(partner_sign('m_z','D',label), -1)
            self.assertEqual(partner_sign("m_z'",'D',label), +1)
    def test_mixed_yz_channel(self):
        self.assertEqual(induced_polarization_sign(1,2),+1)
        self.assertEqual(induced_polarization_sign(1,2,True),-1)
        self.assertIn('yz',allowed_components('D'))
    def test_ahe_parity_complementary(self):
        self.assertEqual(ahe_z_sign(R_MZ,False),+1)
        self.assertEqual(ahe_z_sign(R_MZ,True),-1)
        self.assertEqual(-ahe_z_sign(R_MZ,False),induced_polarization_sign(0,0))
        self.assertEqual(-ahe_z_sign(R_MZ,True),induced_polarization_sign(0,0,True))
    def test_ac_displacement_integrated_charge_zero(self):
        import math
        # Integral of J(t)=-omega*alpha*B0*sin(omega*t) over a full period.
        n=1000
        omega=0.17
        dt=2*math.pi/(omega*n)
        integral=sum(-omega*math.sin(omega*(k+0.5)*dt)*dt for k in range(n))
        self.assertAlmostEqual(integral,0,places=12)

if __name__=='__main__': unittest.main()