"""v55 physics-first orbital embedding, k-dependent frame and finite Peierls B.

r_i=i*d_ki+A_i, energy-velocity d_iH-i[A_i,H]. Under W(k):
H'=W†HW, A_i'=W†A_iW+iW†d_iW, v_i'=W†v_iW.
Finite OPEN graph has straight Peierls links H_ij=t_ij exp(i*q/hbar*int_j^i A.dr)
and analytic magnetic vertex d_B H. Not a full E×B conductivity.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
import numpy as np
from .orbital_geometry_v54 import band_geometry

def _adj(x):
    return np.swapaxes(np.asarray(x).conj(),-1,-2)

def _check_hermitian(M,label):
    if not np.isfinite(M).all() or not np.allclose(M,_adj(M),atol=2e-10,rtol=0):
        raise ValueError(label+' must be finite Hermitian matrices')

def covariant_velocity(H,derivatives,connection):
    """Energy-velocity dH-i[A,H] (physical velocity divides by hbar)."""
    H=np.asarray(H,complex);dH=np.asarray(derivatives,complex)
    A=np.asarray(connection,complex)
    if H.ndim<2 or H.shape[-1]!=H.shape[-2] or H.shape[-1]<2:
        raise ValueError('H dimensions (...,N,N), N>=2 required')
    if dH.shape!=(2,)+H.shape or A.shape!=(2,)+H.shape:
        raise ValueError('dH and orbital connection must each have shape (2,...,N,N)')
    _check_hermitian(H,'H');_check_hermitian(dH,'dH');_check_hermitian(A,'A')
    return dH-1j*(A@H-H@A)

def frame_transform(H,derivatives,connection,W,dW,*,magnetic_vertex=None):
    """Change smooth orbital basis; magnetic vertex W†MW assumes B-independent W."""
    H=np.asarray(H,complex);dH=np.asarray(derivatives,complex)
    A=np.asarray(connection,complex);W=np.asarray(W,complex);dW=np.asarray(dW,complex)
    if H.ndim<2 or H.shape[-1]!=H.shape[-2] or W.shape!=H.shape:
        raise ValueError('H/W shape')
    if dH.shape!=(2,)+H.shape or A.shape!=(2,)+H.shape or dW.shape!=(2,)+H.shape:
        raise ValueError('dH/A/dW shape')
    eye=np.eye(H.shape[-1],dtype=complex)
    if not np.allclose(_adj(W)@W,eye,atol=2e-9,rtol=0):
        raise ValueError('W must be unitary')
    if not np.allclose(_adj(dW)@W+_adj(W)@dW,0,atol=2e-9,rtol=0):
        raise ValueError('dW must differentiate the unitary W')
    _check_hermitian(H,'H');_check_hermitian(dH,'dH');_check_hermitian(A,'A')
    wd=_adj(W);H1=wd@H@W
    d1=np.stack([_adj(dW[i])@H@W+wd@dH[i]@W+wd@H@dW[i] for i in range(2)])
    A1=np.stack([wd@A[i]@W+1j*wd@dW[i] for i in range(2)])
    _check_hermitian(A1,'transformed A')
    M1=None
    if magnetic_vertex is not None:
        M=np.asarray(magnetic_vertex,complex)
        if M.shape!=H.shape:raise ValueError('magnetic vertex shape')
        _check_hermitian(M,'magnetic vertex')
        M1=wd@M@W
    return {'H':H1,'derivatives':d1,'connection':A1,'magnetic_vertex':M1,
            'acceptance':'VERIFIED_FRAME_TRANSFORMATION_FOR_B_INDEPENDENT_W'}

def connected_band_geometry(H,dH,A,*,gap_tol=1e-8):
    """Physically supplied COMPLETE finite orbital connection; no degeneracies."""
    geometry=band_geometry(H,covariant_velocity(H,dH,A),gap_tol=gap_tol)
    return replace(geometry,acceptance='VERIFIED_NONDEGENERATE_COMPLETE_FINITE_ORBITAL_CONNECTION')

def orbital_center_connection(positions):
    """Point orbital embedding only, NOT general Wannier dipole matrix."""
    r=np.asarray(positions,float)
    if r.ndim!=2 or r.shape[1]!=2 or r.shape[0]<2 or not np.isfinite(r).all():
        raise ValueError('orbital positions must be finite (N,2), N>=2')
    return np.stack([np.diag(r[:,axis]).astype(complex) for axis in range(2)])

def _line_phase_derivative(positions,gauge):
    r=np.asarray(positions,float);x=r[:,0];y=r[:,1]
    if gauge=='symmetric':
        return .5*(x[None,:]*y[:,None]-y[None,:]*x[:,None])
    if gauge=='landau_y':
        return .5*(x[:,None]+x[None,:])*(y[:,None]-y[None,:])
    raise ValueError('gauge must be symmetric or landau_y')

@dataclass(frozen=True)
class PeierlsResult:
    hamiltonian:np.ndarray
    magnetic_vertex:np.ndarray
    line_phase_per_B:np.ndarray
    acceptance:str='VERIFIED_FINITE_GRAPH_STRAIGHT_PEIERLS_MAGNETIC_VERTEX'

def peierls_field(hoppings,positions,field,charge_over_hbar,*,gauge='symmetric'):
    """Finite open graph, signed q/hbar, straight-link source j to row i."""
    H0=np.asarray(hoppings,complex);pos=np.asarray(positions,float)
    if H0.ndim!=2 or H0.shape[0]!=H0.shape[1] or H0.shape[0]<2:
        raise ValueError('finite open-boundary Hermitian graph required')
    if pos.shape!=(len(H0),2) or not np.isfinite(pos).all():
        raise ValueError('one physical (x,y) coordinate per graph node required')
    if not np.isfinite([field,charge_over_hbar]).all() or charge_over_hbar==0:
        raise ValueError('field and nonzero q/hbar must be finite')
    _check_hermitian(H0,'hoppings')
    phase_derivative=_line_phase_derivative(pos,gauge)*charge_over_hbar
    H=H0*np.exp(1j*field*phase_derivative)
    M=1j*phase_derivative*H
    _check_hermitian(H,'Peierls Hamiltonian');_check_hermitian(M,'Peierls magnetic vertex')
    return PeierlsResult(H,M,phase_derivative)

def field_gauge_unitary(positions,field,charge_over_hbar):
    """G maps symmetric A to Landau A for Lambda=B*x*y/2."""
    pos=np.asarray(positions,float)
    if pos.ndim!=2 or pos.shape[1]!=2 or not np.isfinite(pos).all():
        raise ValueError('invalid site coordinates')
    if not np.isfinite([field,charge_over_hbar]).all():raise ValueError('invalid magnetic parameters')
    theta=.5*charge_over_hbar*pos[:,0]*pos[:,1]
    G=np.diag(np.exp(1j*field*theta))
    dG=np.diag(1j*theta*np.exp(1j*field*theta))
    return G,dG

def magnetic_vertex_frame_law(H_s,M_s,G,dG):
    """If G depends on B: M'=G M G†+dG H G†+G H dG†; not simply G M G†."""
    H_s=np.asarray(H_s,complex);M_s=np.asarray(M_s,complex)
    G=np.asarray(G,complex);dG=np.asarray(dG,complex)
    if not (H_s.shape==M_s.shape==G.shape==dG.shape):
        raise ValueError('matrix shapes')
    if not np.allclose(_adj(G)@G,np.eye(len(G)),atol=2e-10):
        raise ValueError('G is not unitary')
    return G@M_s@_adj(G)+dG@H_s@_adj(G)+G@H_s@_adj(dG)

def plaquette_wilson(H_B,indices):
    """Link-phase product in indicated polygon order; includes any H0 flux."""
    H=np.asarray(H_B,complex);inds=list(indices)
    if len(inds)<3 or len(set(inds))!=len(inds):
        raise ValueError('simple closed polygon')
    result=1.+0j
    for a in range(len(inds)):
        source=inds[a];destination=inds[(a+1)%len(inds)]
        hop=H[destination,source]
        if abs(hop)<1e-14:raise ValueError('missing polygon hopping')
        result*=hop/abs(hop)
    return result
