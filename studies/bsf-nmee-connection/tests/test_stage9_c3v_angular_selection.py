"""Stage9: C3v transverse ordinary Edelstein vs longitudinal Bz-assisted static symmetry.

These are point-group tensor covariance identities, not predictions of finite
material response amplitudes or driven-DC/static protocol equivalence.
"""
import unittest
import numpy as np

def c3v():
    reflection=np.diag((1.,-1.))
    out=[]
    for j in range(3):
        angle=2*np.pi*j/3
        rot=np.array([[np.cos(angle),-np.sin(angle)],
                      [np.sin(angle),np.cos(angle)]])
        out.extend(((rot,1.),(rot@reflection,-1.)))
    return out

def projected_generators(ordinary=False):
    out=[]
    for i in range(2):
        for j in range(2):
            seed=np.zeros((2,2));seed[i,j]=1.
            avg=sum((det if ordinary else 1.)*(R@seed@R.T)
                    for R,det in c3v())/6
            out.append(avg)
    return out

class C3vAngularSpinSelectionTests(unittest.TestCase):
    def test_conventional_transverse_one_dimensional(self):
        P=projected_generators(ordinary=True)
        self.assertEqual(np.linalg.matrix_rank(np.stack([x.ravel() for x in P]),tol=1e-8),1)
        for x in P:
            self.assertLess(abs(np.trace(x)),1e-12)
            self.assertLess(abs(x[0,1]+x[1,0]),1e-12)
            self.assertLess(max(abs(x[0,0]),abs(x[1,1])),1e-12)
        self.assertGreater(max(abs(x[0,1]) for x in P),.2)

    def test_Bz_assisted_longitudinal_one_dimensional(self):
        P=projected_generators(ordinary=False)
        self.assertEqual(np.linalg.matrix_rank(np.stack([x.ravel() for x in P]),tol=1e-8),1)
        for x in P:
            self.assertLess(abs(x[0,0]-x[1,1]),1e-12)
            self.assertLess(max(abs(x[0,1]),abs(x[1,0])),1e-12)
        self.assertGreater(max(abs(x[0,0]) for x in P),.2)

if __name__=="__main__":
    unittest.main()
