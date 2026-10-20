"""Microscopic displacement-generated interlayer hopping, for hypothesis tests only.

The layer-parity pair has opposite Rashba SOC, connected through physical
horizontal mirror; sublattice-dependent real-space interlayer tunneling is
summed over neighboring in-plane Bravais translations. Neither a DFT Wannier
fit nor a prediction of a stable ferroelectric switching path.

2D k and orbital positions in nm, energies in eV. Bloch gauge includes
absolute orbital displacements; interlayer phase is exp(i k dot d), not
exp(i k dot lattice_vector) alone. The same embedding enters velocities.
"""
from dataclasses import dataclass
import numpy as np
from nmee_lab.clean_core.core import BlochModel
from nmee_lab.clean_core.models import kane_mele, I2, SX, SY, SZ

LAYER_X=np.array([[0,1],[1,0]],complex)
LAYER_Z=np.diag([1,-1]).astype(complex)
SUB_Z=np.diag([1,-1]).astype(complex)
M_SPIN=np.kron(I2,SZ)
M_GLOBAL=np.kron(LAYER_X,M_SPIN)
SA_X=np.kron(np.eye(4),SX)
SA_Y=np.kron(np.eye(4),SY)
SA_Z=np.kron(np.eye(4),SZ)
P_TOP=np.kron(np.diag([1,0]),np.eye(4))
P_BOTTOM=np.kron(np.diag([0,1]),np.eye(4))

@dataclass(frozen=True)
class SlidingSpec:
    shift:tuple[float,float]=(0.23,0.41)   # primitive fractional translations
    tunneling_eV:float=0.028
    decay_nm:float=0.16
    cross_sublattice:bool=True
    lattice_shell:int=5
    mu_eV:float=0.02541381265
    bond_nm:float=0.142
    rashba_soc_eV:float=0.020
    intrinsic_soc_eV:float=0.010


def lattice_geometry(bond_nm=.142):
    nn=np.array([[0,1],[np.sqrt(3)/2,-.5],[-np.sqrt(3)/2,-.5]])*bond_nm
    a=np.stack((nn[0]-nn[1],nn[0]-nn[2]))
    basis=np.stack((np.zeros(2),nn[0]))
    return a,basis


class SlidingBilayer:
    def __init__(self,spec=SlidingSpec()):
        self.spec=spec
        self.mono=kane_mele(bond_nm=spec.bond_nm, rashba_soc=spec.rashba_soc_eV, intrinsic_soc=spec.intrinsic_soc_eV)
        self.a,self.pos=lattice_geometry(spec.bond_nm)
        l=spec.lattice_shell
        lattice_points=np.stack(np.meshgrid(np.arange(-l,l+1),np.arange(-l,l+1),indexing='ij'),-1).reshape(-1,2)
        t_cart=np.asarray(spec.shift)@self.a
        # Upper/lower orbital displacements including shift of lower layer.
        hops=[]
        for ia in range(2):
            for ib in range(2):
                if not spec.cross_sublattice and ia!=ib:continue
                disp=lattice_points@self.a + self.pos[ib]-self.pos[ia]+t_cart
                amp=spec.tunneling_eV*np.exp(-np.einsum('ij,ij->i',disp,disp)/spec.decay_nm**2)
                hops.append((ia,ib,disp,amp))
        self.hops=hops
        self.spins={
          'total_x':SA_X, 'top_x':P_TOP@SA_X, 'bottom_x':P_BOTTOM@SA_X,
          'total_y':SA_Y, 'top_y':P_TOP@SA_Y, 'bottom_y':P_BOTTOM@SA_Y,
          'total_z':SA_Z, 'top_z':P_TOP@SA_Z, 'bottom_z':P_BOTTOM@SA_Z}

    def tunneling(self,k,axis=None):
        k=np.atleast_2d(np.asarray(k,float))
        result=np.zeros((len(k),4,4),complex)
        for ia,ib,disp,amp in self.hops:
            phase=np.exp(1j*k@disp.T)
            if axis is None:
                val=phase@amp
            else:
                val=phase@(1j*disp[:,axis]*amp)
            result[:,2*ia:2*ia+2,2*ib:2*ib+2] += val[:,None,None]*I2
        return result

    def matrices(self,k,axis=None):
        k=np.atleast_2d(np.asarray(k,float))
        h=self.mono.hamiltonian(k) if axis is None else self.mono.velocity(k,axis)
        other=M_SPIN@h@M_SPIN
        t=self.tunneling(k,axis)
        r=np.zeros((len(k),8,8),complex)
        r[:,:4,:4]=h
        r[:,4:,4:]=other
        r[:,:4,4:]=t
        r[:,4:,:4]=t.conj().swapaxes(1,2) if axis is None else t.conj().swapaxes(1,2)
        # Crucial: derivative of T(k)^dagger is (dT/dk)^dagger; yes.
        return r

    def model(self):
        spins={key:(lambda k,mat=mat:np.broadcast_to(mat,(len(np.atleast_2d(k)),8,8)).copy()) for key,mat in self.spins.items()}
        return BlochModel(lambda k:self.matrices(k),lambda k,axis:self.matrices(k,axis),
                          spins,self.mono.reciprocal,g_factor=self.mono.g_factor,
                          name='8-band-sliding-microscopic-hopping')
