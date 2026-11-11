"""Local prediction including a competing quadratic optical threshold.

For the strong mixing test, the lower optical pair has the quartic unfolding
at T=(0,pi), while the upper optical pair has a nondegenerate minimum at Gamma.
Both coefficients are evaluated from the Hamiltonian; no full-BZ response fit.
Gaussian optical broadening only.
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from scipy.special import ndtr
from scipy.optimize import minimize
from .models import bands,kernels,energy,critical
from .lineshapes import phi


def upper_edge(p):
    E,U,V,W=bands(0.,0.,p);Hs=[]
    for n in (0,2):
        H=W[:,:,n,n].real.copy()
        for j in range(3):
            if j!=n:H+=2*np.real(V[:,n,j,None]*V[None,:,j,n])/(E[n]-E[j])
        Hs.append(H)
    H=Hs[1]-Hs[0];ev=np.linalg.eigvalsh(H)
    if min(ev)<=0:raise ValueError('Gamma is not a positive quadratic optical minimum')
    K=float(kernels(0.,0.,p)[1][1]);coef=K/(2*np.pi*np.sqrt(np.linalg.det(H)))
    return dict(energy=float(E[2]-E[0]),hessian=H.tolist(),hessian_eigenvalues=ev.tolist(),kernel=K,coefficient=float(coef),multiplicity=1)


def approximate(p,gamma,xi,c):
    et=float(energy(0.,np.pi,p)[0]);ph=et+xi*gamma
    zeta=c['a_prime']*(p.mass-c['mass'])/np.sqrt(c['B']*gamma)
    pref=c['kernel']*gamma**(-.25)/(np.sqrt(c['A'])*c['B']**.25*(2*np.pi)**2)
    lower=pref*phi(xi,zeta,'gaussian',240)
    e=upper_edge(p);upper=e['coefficient']*ndtr((ph-e['energy'])/gamma)
    return dict(photon=ph,lower=float(lower),upper=float(upper),total=float(lower+upper),xi=float(xi),zeta=float(zeta),edge=e)


def predict_two_feature(p,gamma):
    c=critical(p);scale=np.sqrt(c['B']*gamma)/abs(c['a_prime'])
    def obj(z):
        if not(1.<z[0]<6.) or not(-4.<z[1]<8.):return 100+float(np.sum(np.square(z)))
        pp=replace(p,mass=c['mass']-scale*z[0]);return -approximate(pp,gamma,z[1],c)['total']
    sols=[]
    for init in [[2.84,-.07],[3.7,0.],[4.4,0.]]:
        s=minimize(obj,init,method='Nelder-Mead',options={'xatol':1e-6,'fatol':1e-12,'maxiter':150})
        sols.append(s)
    s=min(sols,key=lambda s:s.fun);pp=replace(p,mass=c['mass']-scale*s.x[0]);res=approximate(pp,gamma,s.x[1],c)
    return {'parameters':vars(pp),'gamma':gamma,'critical':c,'prediction':res,'optimizer_success':bool(s.success),
            'multistart': [{'zeta':float(t.x[0]),'xi':float(t.x[1]),'amplitude':float(-t.fun)} for t in sols],
            'scope':'Sum of a local quadratic-quartic feature and a Gaussian-broadened quadratic threshold; no full-BZ response fitting.'}
