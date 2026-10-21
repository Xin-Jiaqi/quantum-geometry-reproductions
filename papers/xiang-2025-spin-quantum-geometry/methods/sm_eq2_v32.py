"""Nondegenerate multiband SM Eq.(2) four-source audit; independent V14 migration.

Eq.(2) literal with omega_1=omega_2=0, first switch eta, second 2*eta.
A fixed-orbital Liouville recursion is an independent check.
No statement about physical scattering-driven DC or PT degeneracy is made.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm
from scipy.special import expit

@dataclass
class Model:
    N: int = 4
    seed: int = 20251020
    phi_slope: tuple[float,...] | None = None
    phi_offset: tuple[float,...] | None = None

    def __post_init__(self):
        if self.N < 2 or self.N > 7: raise ValueError('2 <= N <= 7')
        rng=np.random.default_rng(self.seed+self.N)
        a=rng.normal(size=(self.N,self.N))+1j*rng.normal(size=(self.N,self.N))
        self.G=(a+a.conj().T)/2*.17
        self.G+=np.diag(np.linspace(-.4,.3,self.N))
        self.centers=np.linspace(-2.4,2.4,self.N)
        self.slopes=np.linspace(-.24,.29,self.N)
        self.curvature=np.linspace(.04,.13,self.N)
        self.slopes_phi=np.array(self.phi_slope if self.phi_slope is not None else [0.]*self.N,dtype=float)
        self.offset_phi=np.array(self.phi_offset if self.phi_offset is not None else [0.]*self.N,dtype=float)
        if len(self.slopes_phi)!=self.N or len(self.offset_phi)!=self.N:raise ValueError('phase length')

    def data(self,k:float,mu:float=.25,T:float=.45):
        if T<=0:raise ValueError('finite T required')
        en=self.centers+self.slopes*k+self.curvature*np.cos(k)
        ep=self.slopes-self.curvature*np.sin(k)
        epp=-self.curvature*np.cos(k)
        f=expit((mu-en)/T)
        fp=-f*(1-f)*ep/T
        fpp=f*(1-f)*((1-2*f)*ep**2/T**2-epp/T)
        phase=self.slopes_phi*k+self.offset_phi
        V=np.diag(np.exp(1j*phase))
        U=expm(-1j*k*self.G)@V
        A=V.conj().T@self.G@V-np.diag(self.slopes_phi)
        r=A.copy();np.fill_diagonal(r,0.)
        return dict(k=k,en=en,ep=ep,epp=epp,f=f,fp=fp,fpp=fpp,U=U,A=A,r=r)

def rho1_band(d,eta):
    if eta<=0:raise ValueError('eta>0 required')
    N=len(d['en']);f=d['f'];fp=d['fp'];en=d['en'];r=d['r']
    out=np.zeros((N,N),complex)
    for m in range(N):
        out[m,m]=fp[m]/eta
        for n in range(N):
            if m==n:continue
            out[m,n]=(f[n]-f[m])*r[m,n]/(-(en[m]-en[n])+1j*eta)
    return out

def orbital_liouville_second(model:Model,k,eta,mu=.25,T=.45,h=2e-5):
    """Different algorithm: fixed orbital rho1 derivative; no bandwise Eq2 terms."""
    c=model.data(k,mu,T)
    def at(q):
        d=model.data(q,mu,T)
        return d['U']@rho1_band(d,eta)@d['U'].conj().T
    d1=(at(k+h)-at(k-h))/(2*h)
    source=c['U'].conj().T@(1j*d1)@c['U']
    en=c['en'];N=len(en);out=np.empty((N,N),complex)
    for m in range(N):
        for n in range(N):
            out[m,n]=source[m,n]/(-(en[m]-en[n])+2j*eta)
    return out

def eq2_four_sources(model:Model,k,eta,mu=.25,T=.45,h=2e-5):
    """Four printed Eq.(2) groups, with iη,2iη retained exactly."""
    c=model.data(k,mu,T)
    f=c['f'];fp=c['fp'];fpp=c['fpp'];r=c['r'];A=c['A'];en=c['en'];N=len(en)
    p=model.data(k+h,mu,T);q=model.data(k-h,mu,T)
    def qband(d,m,n):
        x=d['en'][m]-d['en'][n]
        return 1j*(d['f'][n]-d['f'][m])*d['r'][m,n]/(-x+1j*eta)
    ts=np.zeros((4,N,N),complex)
    for m in range(N):
        for n in range(N):
            x=en[m]-en[n];outden=-x+2j*eta
            if m==n:
                ts[0,m,m]=-fpp[m]/((1j*eta)*(2j*eta))
            else:
                dq=(qband(p,m,n)-qband(q,m,n))/(2*h)
                dq-=1j*(A[m,m]-A[n,n])*qband(c,m,n)
                ts[1,m,n]=dq/outden
                ts[2,m,n]=r[m,n]*(fp[n]-fp[m])/eta/outden
            summand=0j
            for l in range(N):
                f_nl=f[n]-f[l];f_lm=f[l]-f[m]
                summand+=r[m,l]*r[l,n]*(f_nl/(-(en[l]-en[n])+1j*eta)
                    -f_lm/(-(en[m]-en[l])+1j*eta))
            ts[3,m,n]=summand/outden
    return ts

def eta0_eq2_and_printed_eq6(model:Model,k,mu=.25,T=.45,h=2e-5):
    """Two different printed algebra paths: Eq2 eta^0 vs Eq6 as typeset."""
    c=model.data(k,mu,T);p=model.data(k+h,mu,T);q=model.data(k-h,mu,T)
    en=c['en'];f=c['f'];fp=c['fp'];r=c['r'];A=c['A'];N=len(en)
    d2=np.zeros((N,N),complex);f2=d2.copy();tr=d2.copy();missing=d2.copy();drude=d2.copy()
    for m in range(N):
        for n in range(N):
            if m==n:continue
            x=en[m]-en[n];Fp=fp[n]-fp[m]
            def foo(z):
                return (z['f'][n]-z['f'][m])*z['r'][m,n]/(z['en'][m]-z['en'][n])
            df=(foo(p)-foo(q))/(2*h)-1j*(A[m,m]-A[n,n])*foo(c)
            c2=1j*df/x
            c3=-2j*r[m,n]*Fp/x**2
            printed=-1j*r[m,n]*Fp/x**2
            k0=0j
            for l in range(N):
                if l in (m,n):continue
                k0+=r[m,l]*r[l,n]*((f[n]-f[l])/(-(en[l]-en[n]))
                    -(f[l]-f[m])/(-(en[m]-en[l])))/(-x)
            tr[m,n]=k0
            d2[m,n]=c2+c3+k0
            f2[m,n]=c2+printed+k0
            missing[m,n]=c3-printed
            drude[m,n]=-r[m,n]*Fp/x
    return dict(eq2_C=d2,eq6_print_C=f2,triple=tr,extra=missing,
                eq4_leading=drude,bands=en,U=c['U'])

def printed_eq6_triple_independent(model:Model,k,mu=.25,T=.45):
    """Recompute Eq.(6)'s l!=m,n sum without calling the Eq2 triplet code."""
    c=model.data(k,mu,T);en=c['en'];f=c['f'];r=c['r'];N=len(en)
    out=np.zeros((N,N),complex)
    for m in range(N):
        for n in range(N):
            if m==n:continue
            x=en[m]-en[n];s=0j
            for l in range(N):
                if l in (m,n):continue
                s+=r[m,l]*r[l,n]*((f[l]-f[n])/(en[l]-en[n])
                    -(f[m]-f[l])/(en[m]-en[l]))
            out[m,n]=-s/x
    return out

def probe(N=4,seed=20251020,k=.23,eta=.02,mu=.25,T=.45,h=2e-5,
          phi_slope=None,phi_offset=None):
    m=Model(N,seed,phi_slope,phi_offset)
    d=m.data(k,mu,T)
    four=eq2_four_sources(m,k,eta,mu,T,h)
    orbital=orbital_liouville_second(m,k,eta,mu,T,h)
    lc=eta0_eq2_and_printed_eq6(m,k,mu,T,h)
    off=~np.eye(N,dtype=bool)
    z=np.diag(np.array(([.5,-.5]*((N+1)//2))[:N]))
    S=d['U'].conj().T@z@d['U']
    trace=lambda a:float(np.trace(S@a).real)
    return dict(n=N,eta=eta,k=k,
        eq2_vs_orbital_max=float(np.max(abs(four.sum(axis=0)-orbital))),
        eq2_minus_eq6_max=float(np.max(abs((lc['eq2_C']-lc['eq6_print_C']-lc['extra'])[off]))),
        triple_max=float(np.max(abs(lc['triple'][off]))),
        extra_max=float(np.max(abs(lc['extra'][off]))),
        delta_spin_trace=trace(lc['extra']),
        eq2_offdiag_spin_trace=trace(lc['eq2_C']),
        eq6_offdiag_spin_trace=trace(lc['eq6_print_C']),
        full_orbital_spin_trace=trace(orbital),
        protocol='common_eta_double_eta_dc_non_degenerate_multiband',
        gauge=dict(phi_slope=tuple(m.slopes_phi),phi_offset=tuple(m.offset_phi)))
