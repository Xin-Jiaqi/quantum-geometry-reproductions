"""Independent extension of the Zou 2025 three-orbital TB model.

SM IV.5 / S2 / S10, fixed orthonormal orbital basis, one occupied band.
All results use source energy units and dimensionless lattice momentum.
The controlled paths below are NOT claimed to be a real sliding coordinate.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
from scipy.optimize import brentq, root

@dataclass(frozen=True)
class Params:
    mass: float = .0993215269145
    epsilon1: float = 1.03
    beta: float = 0.
    rho: float = 0.
    d1: float = .3
    d2: float = .3
    name: str = 'constant'
    def __post_init__(self):
        if not np.isfinite([self.mass,self.epsilon1,self.beta,self.rho,self.d1,self.d2]).all():
            raise ValueError('All model parameters must be finite')

FAMILIES = {
 'constant': Params(),
 'bond_y': Params(beta=.4,name='bond_y'),
 'mixed_orbital': Params(rho=4.,name='mixed_orbital'),
 'dark': Params(epsilon1=.98,d1=0.,name='dark'),
}

def matrices(x,y,p: Params):
    """H, analytic first and second derivatives, axes preceding band indices."""
    x,y=np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float))
    h=np.zeros(x.shape+(3,3),complex)
    v=np.zeros(x.shape+(2,3,3),complex)
    w=np.zeros(x.shape+(2,2,3,3),complex)
    sx,cx,sy,cy=np.sin(x),np.cos(x),np.sin(y),np.cos(y)
    a=-.75+.5*cy+cx
    h[...,0,0]=a; h[...,1,1]=-a
    h[...,2,2]=-6+p.epsilon1*(-.75+.5*cy)
    h[...,0,1]=p.mass*(1+p.beta*(1+cy))-1j*sx
    h[...,0,2]=-1j*p.d1*sy
    h[...,1,2]=p.rho*p.mass-1j*p.d2*sx
    v[...,0,0,0]=-sx; v[...,0,1,1]=sx
    v[...,0,0,1]=-1j*cx; v[...,0,1,2]=-1j*p.d2*cx
    v[...,1,0,0]=-.5*sy; v[...,1,1,1]=.5*sy
    v[...,1,2,2]=-.5*p.epsilon1*sy
    v[...,1,0,1]=-p.mass*p.beta*sy
    v[...,1,0,2]=-1j*p.d1*cy
    w[...,0,0,0,0]=-cx; w[...,0,0,1,1]=cx
    w[...,0,0,0,1]=1j*sx; w[...,0,0,1,2]=1j*p.d2*sx
    w[...,1,1,0,0]=-.5*cy; w[...,1,1,1,1]=.5*cy
    w[...,1,1,2,2]=-.5*p.epsilon1*cy
    w[...,1,1,0,1]=-p.mass*p.beta*cy
    w[...,1,1,0,2]=1j*p.d1*sy
    def herm(a): return a+np.swapaxes(np.triu(a,1).conj(),-2,-1)
    return herm(h),herm(v),herm(w)

def h_only(x,y,p):
    x,y=np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float))
    h=np.zeros(x.shape+(3,3),complex);a=-.75+.5*np.cos(y)+np.cos(x)
    h[...,0,0]=a;h[...,1,1]=-a;h[...,2,2]=-6+p.epsilon1*(-.75+.5*np.cos(y))
    h[...,0,1]=p.mass*(1+p.beta*(1+np.cos(y)))-1j*np.sin(x)
    h[...,0,2]=-1j*p.d1*np.sin(y);h[...,1,2]=p.rho*p.mass-1j*p.d2*np.sin(x)
    return h+np.swapaxes(np.triu(h,1).conj(),-2,-1)

def energy(x,y,p):
    e=np.linalg.eigvalsh(h_only(x,y,p));return e[...,1:]-e[...,:1]

def bands(x,y,p):
    H,V,W=matrices(x,y,p);E,U=np.linalg.eigh(H);Ud=np.swapaxes(U.conj(),-2,-1)
    v=Ud[...,None,:,:]@V@U[...,None,:,:]
    w=Ud[...,None,None,:,:]@W@U[...,None,None,:,:]
    return E,U,v,w

def derivative(E,v,w,n,m,a,b):
    om=E[...,n]-E[...,m]
    if np.any(np.abs(om)<1e-12): raise ValueError('Unresolved optical gap')
    q=(v[...,b,n,m]*(v[...,a,n,n]-v[...,a,m,m])+v[...,a,n,m]*(v[...,b,n,n]-v[...,b,m,m]))/om-w[...,a,b,n,m]
    dbc=np.zeros_like(q)
    for l in range(3):
        if l in (n,m): continue
        d1=E[...,l]-E[...,m];d2=E[...,n]-E[...,l]
        if np.any(abs(d1)<1e-12) or np.any(abs(d2)<1e-12):
            raise ValueError('Exact degeneracy: no cutoff is inserted')
        t1=v[...,b,n,l]*v[...,a,l,m]/d1
        t2=-v[...,a,n,l]*v[...,b,l,m]/d2
        q=q+t1+t2
        if n in (1,2) and m==0: dbc+=t2
    return 1j*q/om,1j*dbc/om

def kernels_from_bands(E,v,w,axes=(0,0,1)):
    a,b,c=axes; vals=[];dbcs=[];osc=[]
    for n in (1,2):
        en=E[...,n]-E[...,0];rb=1j*v[...,b,0,n]/en;rc=1j*v[...,c,0,n]/en
        D1,T1=derivative(E,v,w,n,0,a,c);D2,T2=derivative(E,v,w,n,0,a,b)
        vals.append(-np.imag(rb*D1+rc*D2));dbcs.append(-np.imag(rb*T1+rc*T2))
        osc.append((abs(v[...,0,0,n]/en)**2+abs(v[...,1,0,n]/en)**2).real)
    return E[...,1:]-E[...,:1],np.stack(vals,axis=-1),np.stack(dbcs,axis=-1),np.stack(osc,axis=-1),E

def kernels(x,y,p,axes=(0,0,1)):
    E,U,v,w=bands(x,y,p);return kernels_from_bands(E,v,w,axes)

def deriv_x(x,y,p):
    e,u=np.linalg.eigh(h_only(x,y,p));dx=np.zeros_like(u)
    dx[...,0,0]=-np.sin(x);dx[...,1,1]=np.sin(x)
    dx[...,0,1]=-1j*np.cos(x);dx[...,1,0]=1j*np.cos(x)
    dx[...,1,2]=-1j*p.d2*np.cos(x);dx[...,2,1]=1j*p.d2*np.cos(x)
    v=np.einsum('...in,...ij,...jn->...n',u.conj(),dx,u).real
    return v[...,1:]-v[...,:1]

def grad_hess(x,y,p):
    E,U,V,W=bands(x,y,p);H=[]
    for n in (0,1):
        value=W[:,:,n,n].real.copy()
        for j in range(3):
            if j!=n: value+=2*np.real(V[:,n,j,None]*V[None,:,j,n])/(E[n]-E[j])
        H.append(value)
    return (V[:,1,1]-V[:,0,0]).real,H[1]-H[0],float(E[1]-E[0])

def critical(p: Params,span=.05):
    def cur(m):return np.linalg.eigvalsh(grad_hess(0.,np.pi,replace(p,mass=m))[1])[-1]
    ms=np.linspace(.001,.65,261);cv=[cur(m) for m in ms]
    ix=np.flatnonzero(np.array(cv[:-1])*np.array(cv[1:])<0)
    if len(ix)!=1: raise ValueError(f'Expected a single critical root, got {len(ix)} for {p}')
    i=ix[0];mc=brentq(cur,ms[i],ms[i+1],xtol=3e-15);pc=replace(p,mass=mc)
    gr,H,ec=grad_hess(0.,np.pi,pc);ev,U=np.linalg.eigh(H);hard,soft=U[:,0],U[:,1]
    s=np.linspace(-span,span,81);f=energy(s*soft[0],np.pi+s*soft[1],pc)[...,0]-ec
    co=np.polynomial.polynomial.polyfit(s/span,f,8)
    A=-ev[0]/2;B=-co[4]/span**4
    ds=1e-5;sl=(cur(mc+ds)-cur(mc-ds))/(4*ds) # derivative of a=lambda_s/2
    K=float(kernels(0.,np.pi,pc)[1][0]);gap=np.linalg.eigvalsh(h_only(0.,np.pi,pc))
    if A<=0 or B<=0 or sl>=0:raise ValueError('Not the intended quadratic-quartic unfolding')
    return dict(name=p.name,mass=mc,energy=ec,A=float(A),B=float(B),a_prime=float(sl),kernel=K,
                hard=hard.tolist(),soft=soft.tolist(),hessian_eigenvalues=ev.tolist(),
                optical_gradient=gr.tolist(),cc_gap=float(gap[2]-gap[1]),span=span,
                parameters=vars(pc))

def stationary_separation(p:Params,c=None):
    c=critical(p) if c is None else c
    center=np.array([0.,np.pi]);soft=np.array(c['soft']);gr,H,et=grad_hess(*center,p)
    a=np.linalg.eigvalsh(H)[-1]/2
    if a<=0:return dict(separation=0.,center_energy=et,points=[center.tolist()])
    t=np.sqrt(a/(2*c['B']));pts=[];ens=[]
    for sign in [-1,1]:
        sol=root(lambda k:grad_hess(*k,p)[0],center+sign*t*soft,tol=1e-10)
        gr,h,e=grad_hess(*sol.x,p)
        if np.linalg.norm(gr)>1e-9 or np.max(np.linalg.eigvalsh(h))>=0:raise RuntimeError('Stationary maximum not found')
        pts.append(sol.x.tolist());ens.append(e)
    return dict(separation=float(np.mean(ens)-et),center_energy=et,points=pts,energies=ens)
