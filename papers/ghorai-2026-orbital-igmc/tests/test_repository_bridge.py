"""Cross-paper contracts. These require the complete repository checkout."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
from numpy.testing import assert_allclose
from orbital_igmc import CuMnAs, geometry, integrate, mesh

PAPERS=Path(__file__).resolve().parents[2]

def load(name, relative):
    path=PAPERS/relative
    if not path.is_file():raise FileNotFoundError(f'Full repository checkout required: {path}')
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module

class RepositoryBridge(unittest.TestCase):
    def test_cumnas_matches_existing_nmee_with_explicit_parameters(self):
        old=load('_orbital_check_nmee','jia-2026-nonlin-magnetoelectric-edelstein/nmee_lab/ptbackend/cumnas.py')
        new=CuMnAs(a_angstrom=3.82,phi=.37)
        q=mesh(8);h,v,_=new.operators(q)
        # Old 'velocity' is dH/dk in eV*nm. a_nm=1 makes it dH/dq here.
        p=old.Model(t=new.t1,tp_ratio=new.t0/new.t1,lam_ratio=new.alpha/new.t1,
                    j_ratio=new.exchange/new.t1,a_nm=1.)
        assert_allclose(h,old.hamiltonian(q,new.phi,p),atol=1e-15)
        for a in (0,1):assert_allclose(v[a],old.velocity(q,new.phi,a,p),atol=1e-15)

    def test_two_band_spin_matches_existing_igmc_after_declared_normalization(self):
        old=load('_orbital_check_igmc','xiang-2025-intrinsic-gyrotropic-magnetic-current/src/igmc_lab/core.py')
        q=mesh(31);h,der=old.dirac_hamiltonian(q[:,0],q[:,1],a=1.)
        pair=old.pair_kernels(h,der,'x','y')
        old_result=old.kernel_response(pair,0,25)
        # In this test M=-sigma_y/2 per unit b; D11 uses B0/2 harmonics.
        g=geometry(h,np.array([der['x'],der['y']]),-.5*old.SY,group_sizes=(1,1))
        for protocol in ('printed_D10','density_A9'):
            new=integrate(g,0,25,protocol=protocol)
            assert_allclose(new['spin'][0],-.25*old_result['sigma_D_per_omega'],atol=2e-14)

if __name__=='__main__':unittest.main()
