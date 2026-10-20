import unittest
import numpy as np
from igmc_lab.tensor_symmetry import (
    IDENTITY,C4Z,M_XZ,transform,allowed_basis_dim,example_generators)

class TensorClassificationTests(unittest.TestCase):
    def test_pure_T_forbids_D_but_not_C(self):
        z=np.zeros((3,3));z[0,1]=1;z[1,0]=-1
        self.assertTrue(np.array_equal(transform(z,IDENTITY,channel='C',antiunitary=True),z))
        self.assertTrue(np.array_equal(transform(z,IDENTITY,channel='D',antiunitary=True),-z))
        self.assertTrue(np.array_equal(transform(z,IDENTITY,channel='C',antiunitary=True,protocol='paper_printed_eq14'),-z))

    def test_SM_reference_generators_4mm1_prime(self):
        gens=example_generators()['4mm1_prime']
        self.assertEqual(allowed_basis_dim(gens,channel='C'),1)
        self.assertEqual(allowed_basis_dim(gens,channel='D'),0)
        # Printed Eq14 incorrectly forces C=0 in this nonmagnetic group.
        self.assertEqual(allowed_basis_dim(gens,channel='C',protocol='paper_printed_eq14'),0)

    def test_SM_reference_generators_4m_prime_m_prime(self):
        gens=example_generators()['4m_prime_m_prime']
        self.assertEqual(allowed_basis_dim(gens,channel='C'),1)
        self.assertEqual(allowed_basis_dim(gens,channel='D'),2)
        # SM Table II allows independent xx=yy and zz in 4m'm'.
        Hall=np.zeros((3,3));Hall[0,1]=1;Hall[1,0]=-1
        displacement=np.diag([1,1,0])
        for r,anti in gens:
            self.assertTrue(np.allclose(transform(Hall,r,antiunitary=anti,channel='C'),Hall))
            self.assertTrue(np.allclose(transform(displacement,r,antiunitary=anti,channel='D'),displacement))

    def test_inversion_forbids_both(self):
        gens=example_generators()['inversion']
        self.assertEqual(allowed_basis_dim(gens,channel='C'),0)
        self.assertEqual(allowed_basis_dim(gens,channel='D'),0)

if __name__=='__main__':unittest.main()
