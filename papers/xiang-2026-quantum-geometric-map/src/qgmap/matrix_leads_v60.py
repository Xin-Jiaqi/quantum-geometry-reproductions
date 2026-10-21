"""v60: independently validated semi-infinite multi-orbital principal-layer surface Green functions.

H_00 is Hermitian and H_01 is an arbitrary complex interlayer hopping.
Retarded g(E+iη)=[z-H00-H01*g*H01†]^-1 uses accelerated Lopez-Sancho
decimation at EXPLICIT eta>0. Physical device-to-surface coupling C may be
nonlocal in device orbital space; Sigma=CgC†, Gamma=i(Sigma-Sigma†).
Magnetic device-frame G(B) requires C'=GC, d_BC'=(d_BG)C.
All magnetic selfenergy-contact terms retained in full Caroli derivatives.
Finite-device algorithm only, not a bulk orbital mixed EB conductivity.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import quad_vec
from .embedding_magnetic_v55 import peierls_field,field_gauge_unitary
from .lead_transport_v57 import four_terminal

def _dag(a):return a.conj().T

def surface_g_matrix(energy,h00,h01,eta=1e-6,tol=1e-12,maxiter=160):
    h=np.asarray(h00,complex);v=np.asarray(h01,complex)
    n=h.shape[0] if h.ndim==2 else 0
    if n<1 or h.shape!=(n,n) or v.shape!=(n,n) or not np.isfinite(h).all() or not np.isfinite(v).all() or not np.allclose(h,_dag(h),atol=1e-12):
        raise ValueError('finite Hermitian h00 and same-size h01 required')
    if not np.isfinite([energy,eta,tol]).all() or eta<=0 or tol<=0 or maxiter<1:
        raise ValueError('positive retarded regulator and convergence settings required')
    z=(energy+1j*eta)*np.eye(n)
    bulk=h.copy();surface=h.copy();a=v.copy();b=_dag(v);ok=False
    for _ in range(maxiter):
        g=np.linalg.inv(z-bulk)
        ag=a@g@b;bg=b@g@a
        surface+=ag;bulk+=ag+bg
        a,b=a@g@a,b@g@b
        if max(np.linalg.norm(a),np.linalg.norm(b))<tol:
            ok=True;break
    if not ok:raise RuntimeError('matrix principal-layer decimation failed')
    gs=np.linalg.inv(z-surface)
    if np.min(np.linalg.eigvalsh(1j*(gs-_dag(gs))))<-5e-9:
        raise RuntimeError('retarded spectral density not positive semidefinite')
    residual=np.linalg.norm(np.linalg.inv(gs)-(z-h-v@gs@_dag(v)))
    if residual>max(1e-7,tol*1e3):raise RuntimeError('surface Dyson residual too large')
    return gs

@dataclass(frozen=True)
class MatrixLead:
    h00:np.ndarray
    h01:np.ndarray
    coupling:np.ndarray
    eta:float=1e-6
    def __post_init__(self):
        h=np.asarray(self.h00,complex);v=np.asarray(self.h01,complex);c=np.asarray(self.coupling,complex)
        if h.ndim!=2 or h.shape[0]<1 or h.shape[0]!=h.shape[1] or v.shape!=h.shape or c.ndim!=2 or c.shape[1]!=len(h):
            raise ValueError('principal-layer and device-to-lead block dimensions')
        if not np.isfinite(h).all() or not np.isfinite(v).all() or not np.isfinite(c).all() or not np.allclose(h,_dag(h),atol=1e-12) or np.linalg.norm(c)==0:
            raise ValueError('Hermitian onsite and finite nonzero coupling required')
        if not np.isfinite(self.eta) or self.eta<=0:raise ValueError('positive retarded eta required')

def selfenergy_matrix(energy,lead,coupling_derivative=None):
    c=np.asarray(lead.coupling,complex)
    dc=np.zeros_like(c) if coupling_derivative is None else np.asarray(coupling_derivative,complex)
    if dc.shape!=c.shape or not np.isfinite(dc).all():raise ValueError('coupling derivative dimensions')
    g=surface_g_matrix(energy,lead.h00,lead.h01,lead.eta)
    sig=c@g@_dag(c);dsig=dc@g@_dag(c)+c@g@_dag(dc)
    gamma=1j*(sig-_dag(sig));dgamma=1j*(dsig-_dag(dsig))
    if np.min(np.linalg.eigvalsh(gamma))<-5e-9:raise RuntimeError('negative lead Γ')
    return sig,dsig,gamma,dgamma

def matrix_caroli(energy,h,m,leads,dcouplings=None):
    """Exact arbitrary Hermitian Γ_a Caroli and dB with Γ'_a and Σ'_a."""
    h=np.asarray(h,complex);m=np.asarray(m,complex);n=h.shape[0] if h.ndim==2 else 0
    if n<2 or h.shape!=(n,n) or m.shape!=(n,n) or not np.allclose(h,_dag(h),atol=1e-11) or not np.allclose(m,_dag(m),atol=1e-11) or len(leads) not in (2,4):
        raise ValueError('finite Hermitian H,M and 2 or 4 leads required')
    if any(l.coupling.shape[0]!=n for l in leads):raise ValueError('lead coupling device dimension')
    if dcouplings is None:dcouplings=[None]*len(leads)
    if len(dcouplings)!=len(leads):raise ValueError('one derivative per lead')
    parts=[selfenergy_matrix(energy,l,d) for l,d in zip(leads,dcouplings)]
    sig=sum((a[0] for a in parts),np.zeros_like(h))
    dsig=sum((a[1] for a in parts),np.zeros_like(h))
    gammas=[a[2] for a in parts];dgammas=[a[3] for a in parts]
    gr=np.linalg.inv(energy*np.eye(n)-h-sig)
    dgr=gr@(m+dsig)@gr
    ga=_dag(gr);dga=_dag(dgr)
    k=len(leads);t=np.zeros((k,k));dt=t.copy()
    for a in range(k):
        for b in range(k):
            if a==b:continue
            x=gammas[a]@gr@gammas[b]@ga
            dx=(dgammas[a]@gr@gammas[b]@ga+
                gammas[a]@dgr@gammas[b]@ga+
                gammas[a]@gr@dgammas[b]@ga+
                gammas[a]@gr@gammas[b]@dga)
            t[a,b]=np.trace(x).real;dt[a,b]=np.trace(dx).real
    return t,dt

def embedded_four_site_leads(hoppings,positions,t=2.,orbital_matrix=False,eta=1e-8):
    n=len(hoppings)
    if n!=4 or np.asarray(positions).shape!=(4,2):raise ValueError('four-site reference graph only')
    leads=[]
    for site in range(4):
        if orbital_matrix:
            h00=np.array([[.15,.31+.06j],[.31-.06j,-.22]])
            h01=np.array([[t,.12+.15j],[.03-.06j,.77*t]],complex)
            c=np.zeros((4,2),complex)
            c[site,:]=[np.sqrt(.55*t/2)*.87,np.sqrt(.55*t/2)*.49j]
        else:
            h00=np.zeros((1,1),complex);h01=np.array([[t]],complex)
            c=np.zeros((4,1),complex);c[site,0]=np.sqrt(.55*t/2)
        leads.append(MatrixLead(h00,h01,c,eta))
    return tuple(leads)

def integrated_matrix_leads(hoppings,positions,leads,*,field=.31,q_over_hbar=-.87,
                            mu=.15,temp=.16,window=18.,epsrel=1e-7,gauge='symmetric'):
    """Real four-terminal Fermi integral; magnetic gauge transforms ALL lead C(B)."""
    pb=peierls_field(hoppings,positions,field,q_over_hbar,gauge=gauge)
    h=pb.hamiltonian;m=pb.magnetic_vertex
    if gauge=='landau_y':
        G,dG=field_gauge_unitary(positions,field,q_over_hbar)
        transformed=[MatrixLead(l.h00,l.h01,G@l.coupling,l.eta) for l in leads]
        dc=[dG@l.coupling for l in leads]
    else:
        transformed=list(leads);dc=None
    low=mu-window*temp;high=mu+window*temp
    if not np.isfinite([mu,temp,window,epsrel]).all() or temp<=0 or window<8 or epsrel<=0:
        raise ValueError('positive thermal/integration parameters')
    eig=list(np.linalg.eigvalsh(h))
    for l in leads:
        v=l.h01;hh=l.h00;tt=np.linalg.norm(v,2)
        if np.allclose(hh,hh[0,0]*np.eye(len(hh)),atol=1e-12) and np.allclose(v,v[0,0]*np.eye(len(v)),atol=1e-12):
            eig.extend([float(hh[0,0].real-2*tt),float(hh[0,0].real+2*tt)])
    points=sorted({float(e) for e in eig if low<e<high})
    count=[0]
    def integrand(E):
        count[0]+=1
        t,dt=matrix_caroli(E,h,m,transformed,dc)
        w=1/(4*temp*np.cosh((E-mu)/(2*temp))**2)
        return np.r_[t.ravel(),dt.ravel()]*w
    nl=len(leads)
    value,error,info=quad_vec(integrand,low,high,points=points,epsrel=epsrel,epsabs=1e-10,limit=600,full_output=True)
    if not info.success:raise RuntimeError('energy quadrature did not converge')
    T=value[:nl*nl].reshape(nl,nl);dT=value[nl*nl:].reshape(nl,nl)
    curr=float(np.max(abs(T.sum(axis=0)-T.sum(axis=1))))
    if curr>max(1e-7,50*error):raise RuntimeError('lead current conservation failure')
    return {'T':T,'dT_dB':dT,'observables':four_terminal(T,dT) if nl==4 else {'two_terminal':float(T[0,1])},
            'quadrature_error':float(error),'evaluations':count[0],'lead_current_conservation_error':curr,
            'status':'VERIFIED_SPECIFIED_MATRIX_LEADS_CONDITIONAL_FINITE_DEVICE_NOT_BULK'}
