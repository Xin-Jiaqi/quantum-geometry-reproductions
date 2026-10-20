"""Main Fig. 1(a), *printed reduced* intrinsic spin Hall Eq. (13).

Code migrated from the separately audited V27 paper-model implementation.
This is a two-band reduced formula with the paper's four-band terms dropped.
It is not an independent full finite-eta/tau^0 driven-DC solver.
Energy is meV, momentum inverse nm, natural e=hbar=1.
"""
from __future__ import annotations
import numpy as np
from numpy.polynomial.legendre import leggauss
from paper_models import SX,SY,SZ,I,fig1a_hamiltonian


def _gap(beta,alpha,vy):
    t=alpha/vy
    if not abs(t)<1:raise ValueError('Requires tilted type-I Dirac cone')
    return abs(beta)*np.sqrt(1.-t*t)


def _contour(mu,s,vx,vy,alpha,beta,nodes):
    """Exact T=0 elliptical Fermi line (used only for Eq13 FS term)."""
    t=alpha/vy; l=1.-t*t
    if abs(mu)<=_gap(beta,alpha,vy):return None
    r=np.sqrt(mu*mu/l-beta*beta)
    theta=2*np.pi*(np.arange(nodes)+0.5)/nodes
    qx=r*np.cos(theta)
    qy=r*np.sin(theta)/np.sqrt(l)-s*t*mu/l
    h=np.sqrt(qx*qx+qy*qy+beta*beta)
    weight=h/(vx*vy*np.sqrt(l))/(2*np.pi*nodes)
    return qx,qy,h,weight


def intrinsic_reduced_eq13(mu,vx,vy,alpha,beta,nr=72,ntheta=144):
    """The *printed reduced* Eq.(13): separate Fermi sea and FS terms.

    The occupied-valence reference integral is zero by ky parity; the integral
    is evaluated over the conduction-electron or valence-hole elliptical pocket.
    No empirical factor of two is applied to paper's Fig.1(a) ordinate.
    """
    if vx<=0 or vy<=0 or nr<8 or ntheta<16:
        raise ValueError('Invalid velocity or integration mesh')
    if abs(mu)<=_gap(beta,alpha,vy):
        return {'sea':0.,'surface':0.,'total':0.}
    b=1 if mu>0 else -1
    t=alpha/vy;l=1.-t*t
    radius=np.sqrt(mu*mu/l-beta*beta)
    z,w=leggauss(nr)
    theta=2*np.pi*(np.arange(ntheta)+.5)/ntheta
    rr=(z[:,None]+1)*radius/2
    qx=rr*np.cos(theta)[None,:]
    sea=0.;surface=0.
    for s in (-1,1):
        qy=rr*np.sin(theta)[None,:]/np.sqrt(l)-s*t*mu/l
        h=np.sqrt(qx*qx+qy*qy+beta*beta)
        # Q = 2*Delta_mn^y*g_spin_nm^{yx;y} + barDelta_mn*g_charge_nm^{xy}
        g_spin=s*vx*alpha*(qy*qy+beta*beta)/(8*h**4)
        delta_mn=-2*b*vy*qy/h
        bar_delta=-s*alpha*b*qx/h
        g_charge=-vx*vy*qx*qy/(4*h**4)
        Q=2*delta_mn*g_spin+bar_delta*g_charge
        integrand=b*Q/(4*h*h)
        jacobian=rr/(vx*vy*np.sqrt(l))
        sea+=float(np.sum(np.mean(integrand*jacobian,axis=1)*w)*radius/2/(2*np.pi))
        qx_s,qy_s,hs,weights=_contour(mu,s,vx,vy,alpha,beta,ntheta)
        g_yy=-s*alpha*vy*qx_s*qy_s/(8*hs**4)
        eps_nm=2*b*hs
        vnx=b*vx*qx_s/hs
        surface+=float(np.sum(weights*(2*g_yy/eps_nm)*(-vnx)))
    return {'sea':sea,'surface':surface,'total':sea+surface}


def kernel_matrix_at_k(kx,ky,s,b,vx,vy,alpha,beta):
    """Independent Bloch-eigenstate matrix-element check of the Eq13 Q term."""
    H=fig1a_hamiltonian(kx,ky,s,vx,vy,alpha,beta)
    e,U=np.linalg.eigh(H)
    n=1 if b>0 else 0;m=1-n
    vxop=vx*SY;vyop=-s*vy*SX+s*alpha*I
    jyop=(vyop@SY+SY@vyop)/4
    Vx=U.conj().T@vxop@U
    Vy=U.conj().T@vyop@U
    Jy=U.conj().T@jyop@U
    eps=e[n]-e[m]
    g_spin=np.real(Jy[n,m]*Vx[m,n])/eps**2
    g_charge=np.real(Vx[n,m]*Vy[m,n])/eps**2
    g_yy=np.real(Jy[n,m]*Vy[m,n])/eps**2
    delta_mn=Vy[m,m].real-Vy[n,n].real
    bar_delta=Jy[m,m].real-Jy[n,n].real
    Q=2*delta_mn*g_spin+bar_delta*g_charge
    return {'Q_over_eps2':float(Q/eps**2),'g_yy':float(g_yy)}


def kernel_analytic_at_k(kx,ky,s,b,vx,vy,alpha,beta):
    """Independent Pauli-basis algebra for the matrix oracle above."""
    qx=vx*kx;qy=vy*ky
    h=np.sqrt(qx*qx+qy*qy+beta*beta)
    gs=s*vx*alpha*(qy*qy+beta*beta)/(8*h**4)
    gc=-vx*vy*qx*qy/(4*h**4)
    delta_mn=-2*b*vy*qy/h
    bar_delta=-s*alpha*b*qx/h
    Q=2*delta_mn*gs+bar_delta*gc
    return {'Q_over_eps2':float(Q/(4*h*h)),
            'g_yy':float(-s*alpha*vy*qx*qy/(8*h**4))}
