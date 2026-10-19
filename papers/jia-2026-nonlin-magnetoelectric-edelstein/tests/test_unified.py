import os
import tempfile
import unittest
from dataclasses import replace
import numpy as np
from nmee_lab import (CleanSpec,CuMnAsBGKSpec,WannierTB2D,
                      read_wannier90_hr,to_clean_model,projector_from_indices,
                      evaluate_clean,evaluate_cumnas_bgk,mesh,transform_layer_tensors,connection_residual)
from nmee_lab.clean_core.models import kane_mele,SX,SY,SZ,I2
from nmee_lab.ptbackend import conserving_bgk,cumnas


def make_tb():
    rng=np.random.default_rng(314)
    A=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    H0=np.diag([-1.1,-.15,.65,1.45])+.055*(A+A.conj().T)
    Tx=.035*(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    Ty=.025*(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    r=np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1]])
    h=np.array([H0,Tx,Tx.conj().T,Ty,Ty.conj().T])
    lattice=np.array([[.45,0],[.08,.39]])
    return WannierTB2D(r,h,lattice,'artificial-four-band')


def onsite_spins():
    return dict(x=np.kron(I2,SX),y=np.kron(I2,SY),z=np.kron(I2,SZ))


class TestUnified(unittest.TestCase):
    def test_hr_velocity_central_difference(self):
        tb=make_tb();k=np.array([[.45,1.19],[3.,-2.3]])
        for j in (0,1):
            dk=np.zeros_like(k);dk[:,j]=1e-5
            vfd=(tb.matrices(k+dk)-tb.matrices(k-dk))/(2e-5)
            self.assertLess(np.max(abs(vfd-tb.matrices(k,j))),1e-9)
        self.assertTrue(np.allclose(tb.matrices(k),tb.matrices(k).conj().transpose(0,2,1)))

    def test_hr_roundtrip_with_degeneracy(self):
        tb=make_tb();n=tb.nwann
        with tempfile.TemporaryDirectory() as d:
            path=os.path.join(d,'hr.dat')
            deg=[2,1,1,1,1]
            with open(path,'w') as f:
                f.write('synthetic test, no real material\n4\n5\n')
                f.write(' '.join(map(str,deg))+'\n')
                for R,H,weight in zip(tb.translations,tb.hopping,deg):
                    for j in range(n):
                        for i in range(n):
                            v=H[i,j]*weight
                            f.write(f'{R[0]} {R[1]} 0 {i+1} {j+1} {v.real:.16e} {v.imag:.16e}\n')
            tb2=read_wannier90_hr(path,tb.lattice_nm)
            self.assertLess(np.max(abs(tb2.hopping-tb.hopping)),1e-13)
            self.assertLess(np.max(abs(tb2.matrices(np.array([[.3,2.]]))-tb.matrices(np.array([[.3,2.]])))),1e-13)
            # Never silently accept non-2D Wannier hoppings.
            from pathlib import Path
            edited=Path(path).read_text().replace('1 0 0 1 1','1 0 1 1 1',1)
            Path(path).write_text(edited)
            with self.assertRaises(ValueError):read_wannier90_hr(path,tb.lattice_nm)

    def test_hermiticity_missing_pair_and_duplicate_R(self):
        tb=make_tb()
        with self.assertRaises(ValueError):WannierTB2D(tb.translations[:-1],tb.hopping[:-1],tb.lattice_nm)
        with self.assertRaises(ValueError):WannierTB2D(np.array([[0,0],[0,0]]),tb.hopping[:2],tb.lattice_nm)

    def test_explicit_wannier_physics_gate(self):
        tb=make_tb()
        with self.assertRaises(ValueError):to_clean_model(tb,spin_onsite=onsite_spins())
        with self.assertRaises(ValueError):to_clean_model(tb,spin_onsite={'x':onsite_spins()['x']},acknowledge_zero_connection=True)
        with self.assertRaises(ValueError):to_clean_model(tb,spin_onsite=onsite_spins(),
                      layer_indices={'top':[0,1],'bottom':[1,2]},acknowledge_zero_connection=True)

    def test_layer_sum_for_actual_integrated_response(self):
        tb=make_tb()
        m=to_clean_model(tb,spin_onsite=onsite_spins(),layer_indices={'top':[0,1],'bottom':[2,3]},
                         acknowledge_zero_connection=True)
        k=mesh(m.reciprocal,14)
        value={a:evaluate_clean(m,k,CleanSpec(alpha=a,gamma='x',beta=0,mu_eV=0.1))['result']
               for a in ('z','top_z','bottom_z')}
        self.assertLess(abs(value['top_z']+value['bottom_z']-value['z']),1e-10)
        self.assertTrue(np.isfinite(value['z']))

    def test_unitary_basis_transform_of_wannier_H_and_spin(self):
        tb=make_tb();ops=onsite_spins(); rng=np.random.default_rng(55)
        z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        W,_=np.linalg.qr(z)
        tro=lambda x:W.conj().T@x@W
        rot=WannierTB2D(tb.translations,np.array([tro(h) for h in tb.hopping]),tb.lattice_nm)
        m=to_clean_model(tb,spin_onsite=ops,acknowledge_zero_connection=True)
        mr=to_clean_model(rot,spin_onsite={x:tro(y) for x,y in ops.items()},acknowledge_zero_connection=True)
        k=mesh(m.reciprocal,12)
        for branch in ('printed','full_static'):
            s=CleanSpec(alpha='z',gamma='x',mu_eV=.1,branch=branch)
            x=evaluate_clean(m,k,s)['result']; y=evaluate_clean(mr,k,s)['result']
            self.assertLess(abs(x-y),1e-10)

    def test_non_degenerate_gate_inherited(self):
        tb=make_tb(); h=tb.hopping.copy();h[0]=np.zeros((4,4));h[1:]=0
        model=to_clean_model(WannierTB2D(tb.translations,h,tb.lattice_nm),
                       spin_onsite=onsite_spins(),acknowledge_zero_connection=True)
        with self.assertRaises(ValueError):evaluate_clean(model,mesh(model.reciprocal,3),CleanSpec())

    def test_prior_KM_clean_benchmark(self):
        model=kane_mele();p=mesh(model.reciprocal,15)
        full=evaluate_clean(model,p,CleanSpec(alpha='x',gamma='z',mu_eV=.025))
        pub=evaluate_clean(model,p,CleanSpec(alpha='x',gamma='z',mu_eV=.025,branch='printed'))
        self.assertAlmostEqual(full['parts']['printed'],pub['result'],places=12)
        self.assertAlmostEqual(full['parts']['full_static'],full['parts']['printed']+full['parts']['diagonal_disputed'],places=12)
        self.assertEqual(full['protocol'],'nondegenerate_clean_band_kernel')

    def test_CuMnAs_BGK_adapter_equivalent_to_prior(self):
        spec=CuMnAsBGKSpec(theta_rad=.15,mu_eV=.03,tau_anisotropy=.4,
                          observable='total_x',magnetic_component='y',electric_beta=1)
        n=27
        wrapped=evaluate_cumnas_bgk(spec,n=n)
        p=cumnas.Model()
        direct=conserving_bgk.response(cumnas.points(n),.15,.03,p,alpha=.4,
                   spin=cumnas.IY,beta=1,obs=cumnas.IX)
        self.assertLess(abs(wrapped['result']-direct['gamma_bgk']),1e-12)
        self.assertEqual(wrapped['protocol'],'cumnas_PT_BGK_static_B_equilibrated_then_DC_E')

    def test_CuMnAs_layer_sublattice_sum(self):
        sa=CuMnAsBGKSpec(); sb=replace(sa,observable='B_z')
        a=evaluate_cumnas_bgk(sa,n=39)['result'];b=evaluate_cumnas_bgk(sb,n=39)['result']
        self.assertLess(abs(a+b),1e-8)

    def test_layer_connection_operator_inversion_switch(self):
        t=np.zeros((3,3,3));t[0,2,0]=1.25
        a={'top':t,'bottom':-t}
        transformed=transform_layer_tensors(a,-np.eye(3),{'top':'bottom','bottom':'top'})
        self.assertLess(np.max(abs(transformed['top']-a['top'])),1e-14)
        self.assertLess(connection_residual(a,-np.eye(3),{'top':'bottom','bottom':'top'}),1e-14)
        self.assertLess(np.max(abs(a['top']+a['bottom'])),1e-14)
        with self.assertRaises(ValueError):transform_layer_tensors(a,np.eye(3),{'top':'top','bottom':'top'})

    def test_projector_base_index(self):
        a=projector_from_indices(4,[0,2]);b=projector_from_indices(4,[1,3]);
        self.assertTrue(np.allclose(a+b,np.eye(4)))
        with self.assertRaises(ValueError):projector_from_indices(3,[0,0])


if __name__=='__main__':unittest.main()
