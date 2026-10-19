"""Insulating static interband E*B spin response with a nonzero Wannier Berry connection.

This is a *thermodynamic/interband static* diagnostic, not a derivation of the
finite-temperature driven-dc intrinsic response of Jia et al. (2026).

The model uses a smooth common orbital gauge with r_beta=i∂_{k_beta}+A_beta(k),
D_beta O = ∂_{k_beta} O - i[A_beta,O], and H_B=H-g*mu_B*B*S_gamma.
Here A has units of nm, H eV, k nm^-1, O dimensionless Pauli spin.

A(k) and its gauge must be supplied explicitly.  In particular,
Wannier90 hr.dat alone is insufficient: it does not define A(k) or S(k).
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Mapping
import numpy as np
from .clean_core.core import MU_B,chi_k

MatrixCallback=Callable[[np.ndarray],np.ndarray]
DerivativeCallback=Callable[[np.ndarray,int],np.ndarray]

@dataclass(frozen=True)
class CovariantBloch2D:
    hamiltonian:MatrixCallback
    derivative_h:DerivativeCallback
    connection:DerivativeCallback
    operators:Mapping[str,MatrixCallback]
    derivative_operators:Mapping[str,DerivativeCallback]
    reciprocal_nm:np.ndarray
    g_factor:float=10.
    label:str='covariant-synthetic'


def _herm(a,name):
    a=np.asarray(a,dtype=complex)
    if a.ndim!=3 or a.shape[1]!=a.shape[2] or not np.isfinite(a).all() or not np.allclose(a,a.conj().transpose(0,2,1),atol=1e-9):
        raise ValueError(f'{name} must be finite Hermitian (Nk,N,N)')
    return a


def covariant_derivative(op,d_op,conn):
    """D O = partial O - i[A,O], with r=i partial + A convention."""
    return np.asarray(d_op,complex)-1j*(conn@op-op@conn)


def evaluate_insulating_static(model:CovariantBloch2D, points:np.ndarray,*,
                               observable='x', magnetic='z', beta=0,
                               n_occupied:int, B_step_T:float=.02,
                               min_direct_gap_eV:float=1e-6):
    """Central-field derivative of occupied/unoccupied static interband susceptibility.

    Requires a *gapped, T=0 insulator*: all occupied bands are filled, all
    unoccupied empty. Degeneracy *within the same occupied group* is permitted.
    This approach is k-dependent orbital-gauge covariant if all position
    connections, observables and their derivatives share the same gauge.
    NOT a microscopic finite-T dc collision, extrinsic, or Fermi-surface solver.
    """
    pts=np.asarray(points,float)
    if pts.ndim!=2 or pts.shape[1]!=2 or len(pts)==0:
        raise ValueError('points must have shape (Nk,2)')
    if beta not in (0,1):raise ValueError('beta must be 0/1')
    if B_step_T<=0 or min_direct_gap_eV<=0:raise ValueError('positive B step / gap required')
    rec=np.asarray(model.reciprocal_nm,float)
    if rec.shape!=(2,2) or not abs(np.linalg.det(rec))>1e-12:raise ValueError('invalid reciprocal lattice')
    H=_herm(model.hamiltonian(pts),'H')
    dH=_herm(model.derivative_h(pts,beta),'dH')
    A=_herm(model.connection(pts,beta),'A')
    if not all(x in model.operators for x in (observable,magnetic)):
        raise ValueError('missing observable or magnetic operator')
    if magnetic not in model.derivative_operators:
        raise ValueError('provide dSgamma/dk explicitly (zero allowed only if physically constant)')
    obs=_herm(model.operators[observable](pts),'spin observable')
    mag=_herm(model.operators[magnetic](pts),'Zeeman spin')
    dmag=_herm(model.derivative_operators[magnetic](pts,beta),'dSgamma')
    if H.shape!=dH.shape or any(m.shape!=H.shape for m in (A,obs,mag,dmag)):
        raise ValueError('operator shapes inconsistent')
    bands=H.shape[-1]
    if not 0<n_occupied<bands:raise ValueError('n_occupied outside (0,N)')
    V0=covariant_derivative(H,dH,A)
    dmagcov=covariant_derivative(mag,dmag,A)
    _herm(V0,'covariant velocity');_herm(dmagcov,'covariant derivative spin')
    measure=abs(float(np.linalg.det(rec)))/(2*np.pi)**2/len(pts)
    values=[];gaps=[]
    for sgn in (+1.,-1.):
        energyB=H-sgn*B_step_T*model.g_factor*MU_B*mag
        velocityB=V0-sgn*B_step_T*model.g_factor*MU_B*dmagcov
        e,U=np.linalg.eigh(energyB)
        direct=np.min(e[:,n_occupied]-e[:,n_occupied-1])
        global_gap=float(np.min(e[:,n_occupied])-np.max(e[:,n_occupied-1]))
        if direct<min_direct_gap_eV or global_gap<min_direct_gap_eV:
            raise ValueError(f'T=0 insulating method unsafe: direct={direct:.5g}, global={global_gap:.5g} eV')
        occ=np.zeros_like(e);occ[:,:n_occupied]=1.
        values.append(float(measure*np.sum(chi_k(e,U,velocityB,obs,occ,degeneracy_tol=1e-12))))
        gaps.append((float(direct),global_gap))
    return dict(protocol='T0_gapped_static_interband_B_derivative_with_position_connection',
                response=(values[0]-values[1])/(2*B_step_T),
                chi_plus=values[0],chi_minus=values[1],
                min_gap_eV=min(g[0] for g in gaps),global_gap_eV=min(g[1] for g in gaps),
                mesh_n=len(pts),magnetic_step_T=float(B_step_T),
                normalization='legacy Pauli magnetic spin / electric field and B, (nm T V)^-1 for 2D',
                warnings=['static T=0 interband only, not experimental driven-DC intrinsic',
                          'Berry connection and dS/dk are caller supplied, with consistent gauge',
                          'finite magnetic-field step must be converged',
                          'Wannier hr.dat alone lacks projected position/spin operator matrices'])


def covariant_from_wannier(tb, *,position_connections:Mapping[int,DerivativeCallback],
                           spin_operators:Mapping[str,MatrixCallback],
                           spin_derivatives:Mapping[str,DerivativeCallback],
                           g_factor=10.,label=None):
    """Adapter requiring physical projected A_beta(k) in nm and spin callbacks.

    No default A=0; no extrapolation from hr.dat. No automatic r.dat parser.
    A callback must take (k,beta) and provide Hermitian (Nk,N,N) arrays.
    """
    if set(position_connections)!={0,1}:
        raise ValueError('two position connections A_x,A_y required')
    if not spin_operators or not spin_derivatives:
        raise ValueError('spin and spin-derivative callbacks required')
    def connection(k,b):
        return position_connections[b](k,b)
    return CovariantBloch2D(hamiltonian=tb.matrices,
         derivative_h=lambda k,b:tb.matrices(k,b),connection=connection,
         operators=spin_operators,derivative_operators=spin_derivatives,
         reciprocal_nm=tb.reciprocal_nm,g_factor=g_factor,
         label=label or tb.label)
