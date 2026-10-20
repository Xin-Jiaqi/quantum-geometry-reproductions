"""v23 independent homogeneous length-gauge two-color EB Liouville kernel.

For periodic orbital Bloch gauge (r_y=i partial_ky+A_y):
 rho_E = R_E(i D_y rho0), rho_B=R_B([Z,rho0]),
 rho_EB=R_sum(i D_y rho_B + [Z,rho_E]),
 R(w,eta)(X)=(w+i eta-[H,.])^-1 X, D_y X=partial_ky X-i[A_y,X].

Electric field and bbar=Zeeman energy are complex Fourier AMPLITUDES, with no 1/2.
No collision integral, length/time SI conversion, physical experimental DC, or exact
interacting Ward statement is implied. All quantities dimensionless; this is an
independent periodic *lattice* regularization rather than the paper's continuum Dirac.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np

I=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)

@dataclass(frozen=True)
class UniformSpec:
    n:int=48
    beta:float=.24
    alpha:float=.13
    wilson:float=.35
    mass:float=.32
    mu:float=.35
    temp:float=.12
    omega_e:float=0.
    omega_b:float=0.
    eta_e:float=.06
    eta_b:float=.03
    gauge_amplitude:float=0.
    def __post_init__(self):
        if self.n<8 or self.n%2 or self.temp<=0 or self.eta_e<=0 or self.eta_b<=0:
            raise ValueError('even n>=8, positive T and regulators required')


def comm(a,b):return a@b-b@a

def kmesh(n):
    k=2*np.pi*np.arange(n)/n
    return np.stack(np.meshgrid(k,k,indexing='ij'),axis=-1)

def h0_and_vel(spec:UniformSpec):
    k=kmesh(spec.n);x=k[...,0];y=k[...,1]
    s=2-np.cos(x)-np.cos(y)
    h=(spec.beta*np.sin(x)+spec.alpha*s)[...,None,None]*I
    h=h+np.sin(x)[...,None,None]*SY-np.sin(y)[...,None,None]*SX
    h=h+(spec.mass+spec.wilson*s)[...,None,None]*SZ
    vx=(spec.beta*np.cos(x)+spec.alpha*np.sin(x))[...,None,None]*I+np.cos(x)[...,None,None]*SY+(spec.wilson*np.sin(x))[...,None,None]*SZ
    vy=(spec.alpha*np.sin(y))[...,None,None]*I-np.cos(y)[...,None,None]*SX+(spec.wilson*np.sin(y))[...,None,None]*SZ
    return h,vx,vy

def fft_dy(a):
    n=a.shape[1]
    modes=np.fft.fftfreq(n,1/n)
    shape=(1,n)+(1,)*(a.ndim-2)
    return np.fft.ifft(1j*modes.reshape(shape)*np.fft.fft(a,axis=1),axis=1)

def fermi(e,mu,temp):
    return 1/(1+np.exp(np.clip((e-mu)/temp,-500,500)))

def diag_to_orb(u,diag):
    return (u*diag[...,None,:])@np.swapaxes(u.conj(),-1,-2)

def orb_to_band(u,mat):
    ud=np.swapaxes(u.conj(),-1,-2)
    return ud@mat@u

def resolvent(e,u,source,omega,eta):
    x=orb_to_band(u,source)
    den=omega+1j*eta-(e[..., :,None]-e[...,None,:])
    return u@(x/den)@np.swapaxes(u.conj(),-1,-2)

def gauge_transform(spec,h,vx,vy,z):
    """U=e^{i a sin ky sigma_z}; transformed position connection A'= -a cos ky sigma_z."""
    a=spec.gauge_amplitude
    if not a:return h,vx,vy,z,np.zeros_like(h)
    y=kmesh(spec.n)[...,1]
    phi=a*np.sin(y)
    u=np.zeros_like(h);u[...,0,0]=np.exp(1j*phi);u[...,1,1]=np.exp(-1j*phi)
    ud=np.swapaxes(u.conj(),-1,-2)
    tr=lambda m:ud@m@u
    A=-(a*np.cos(y))[...,None,None]*SZ
    return tr(h),tr(vx),tr(vy),tr(z),A

def inputs(spec):
    h,vx,vy=h0_and_vel(spec)
    z=np.broadcast_to(SX,h.shape).copy()
    h,vx,vy,z,conn=gauge_transform(spec,h,vx,vy,z)
    return h,vx,vy,z,conn

def cov_dy(mat,conn):
    return fft_dy(mat)-1j*comm(conn,mat)

def rho_equil(h,spec,zero_occupations=None):
    e,u=np.linalg.eigh(h)
    fn=fermi(e,spec.mu,spec.temp) if zero_occupations is None else zero_occupations
    return diag_to_orb(u,fn),e,u,fn

def source_and_channels(spec:UniformSpec):
    h,vx,vy,z,conn=inputs(spec)
    rho,e,u,f=rho_equil(h,spec)
    de=cov_dy(rho,conn)
    RE=resolvent(e,u,1j*de,spec.omega_e,spec.eta_e)
    RB=resolvent(e,u,comm(z,rho),spec.omega_b,spec.eta_b)
    o_eb=resolvent(e,u,1j*cov_dy(RB,conn),spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    o_be=resolvent(e,u,comm(z,RE),spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    val=lambda x:complex(np.real_if_close(np.mean(np.trace(vx@x,axis1=-2,axis2=-1))))
    return dict(b_then_e=val(o_eb),e_then_b=val(o_be),total=val(o_eb+o_be),
                rho_e=RE,rho_b=RB,rho_eb=o_eb+o_be,
                input_h=h,charge_current=vx,occ=f,gap_min=float(np.min(e[...,1]-e[...,0])))

def leading_fermi_pole(spec:UniformSpec):
    """Analytic coefficient lim_{eta_E=eta_B->0} eta_E * j_EB.

    Only the intraband E-drive contributes to the pole. Magnetic amplitude
    here multiplies Z=Pauli sigma_x, NOT B in Tesla.
    """
    h,vx,vy,z,conn=inputs(spec)
    e,u=np.linalg.eigh(h)
    f=fermi(e,spec.mu,spec.temp)
    vd=np.real(np.diagonal(orb_to_band(u,vy),axis1=-2,axis2=-1))
    dfn=-f*(1-f)*vd/spec.temp
    zb=orb_to_band(u,z)
    delta=e[..., :,None]-e[...,None,:]
    dfn_nm=dfn[...,None,:]-dfn[...,:,None]
    mask=np.eye(2,dtype=bool)
    # Exactly remove diagonal before dividing; otherwise 1/eta amplifies roundoff.
    source=np.where(mask,0.,zb*dfn_nm)
    rho=np.where(mask,0.,source/np.where(mask,1.,-delta))
    vx_band=orb_to_band(u,vx)
    return complex(np.mean(np.trace(vx_band@rho,axis1=-2,axis2=-1)))

def prepared_electric(spec:UniformSpec,b=0.,frozen=False):
    """Static B-prepared electric susceptibility, linear in E with eta_E.
    If frozen=True, follow zero-field band occupations instead of re-thermalizing.
    """
    if spec.omega_b!=0:raise ValueError('prepared B comparison requires omega_b=0')
    h,vx,vy,z,conn=inputs(spec)
    f0=np.linalg.eigh(h)[0]
    frozen_f=fermi(f0,spec.mu,spec.temp) if frozen else None
    hb=h+b*z
    rho,e,u,_=rho_equil(hb,spec,zero_occupations=frozen_f)
    er=resolvent(e,u,1j*cov_dy(rho,conn),spec.omega_e,spec.eta_e)
    return complex(np.mean(np.trace(vx@er,axis1=-2,axis2=-1)))

def prepared_slope(spec:UniformSpec,b=.0005,frozen=False):
    return (prepared_electric(spec,b,frozen)-prepared_electric(spec,-b,frozen))/(2*b)

def population_preparation_contribution(spec:UniformSpec):
    """Difference of prepared equilibrium and prepared frozen, analytic at B=0.
    The derivative f'(e_n)<n|Z|n> enters through i D_y (drho_pop),
    resolvent R_E on both spin bands.
    """
    h,vx,vy,z,conn=inputs(spec)
    e,u=np.linalg.eigh(h)
    f=fermi(e,spec.mu,spec.temp)
    zd=np.real(np.diagonal(orb_to_band(u,z),axis1=-2,axis2=-1))
    pop=diag_to_orb(u,-f*(1-f)*zd/spec.temp)
    rho=resolvent(e,u,1j*cov_dy(pop,conn),spec.omega_e,spec.eta_e)
    return complex(np.mean(np.trace(vx@rho,axis1=-2,axis2=-1)))

def divided_difference_Drho(e,u,vy,spec):
    """Analytic Fréchet derivative D_y f(H) in occupied band basis, no FFT."""
    f=fermi(e,spec.mu,spec.temp)
    df=-f*(1-f)/spec.temp
    dife=e[..., :,None]-e[...,None,:]
    diff=f[..., :,None]-f[...,None,:]
    deriv=np.empty_like(dife)
    mask=np.abs(dife)>1e-10
    deriv[mask]=diff[mask]/dife[mask]
    deriv[~mask]=np.broadcast_to(df[..., :,None],dife.shape)[~mask]
    return u@(deriv*orb_to_band(u,vy))@np.swapaxes(u.conj(),-1,-2)

def algebraic_channels(spec:UniformSpec):
    """No spatial numerical derivatives; differentiates magnetic resolvent algebraically.

    D R_B[Z,rho0] = R_B([Z,D rho0]+[D H, R_B[Z,rho0]])
    for k-independent physical Zeeman Z and D Z=0.
    """
    h,vx,vy,z,conn=inputs(spec)
    rho,e,u,_=rho_equil(h,spec)
    dr=divided_difference_Drho(e,u,vy,spec)
    re=resolvent(e,u,1j*dr,spec.omega_e,spec.eta_e)
    rb=resolvent(e,u,comm(z,rho),spec.omega_b,spec.eta_b)
    drb=resolvent(e,u,comm(z,dr)+comm(vy,rb),spec.omega_b,spec.eta_b)
    be=resolvent(e,u,1j*drb,spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    eb=resolvent(e,u,comm(z,re),spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    val=lambda x:complex(np.mean(np.trace(vx@x,axis1=-2,axis2=-1)))
    return dict(b_then_e=val(be),e_then_b=val(eb),total=val(be+eb))


def integrated_e_charge(spec:UniformSpec):
    """Single E response current (no Zeeman drive)."""
    h,vx,vy,z,conn=inputs(spec)
    rho,e,u,_=rho_equil(h,spec)
    re=resolvent(e,u,1j*cov_dy(rho,conn),spec.omega_e,spec.eta_e)
    return complex(np.mean(np.trace(vx@re,axis1=-2,axis2=-1)))

def eq20_direct_band_current(spec:UniformSpec):
    """Direct SM Eq.(20) rho_mn=-tau Z_mn partial_y f_nm/eps_mn.

    Return current / tau, with epsilon_mn=epsilon_m-epsilon_n and
    rho_mn V_nm ordering. Bbar=one energy unit. Do not silently assign SI units.
    """
    h,vx,vy,z,conn=inputs(spec)
    e,u=np.linalg.eigh(h)
    f=fermi(e,spec.mu,spec.temp)
    vyb=orb_to_band(u,vy)
    zb=orb_to_band(u,z)
    vxb=orb_to_band(u,vx)
    df=-f*(1-f)*np.real(np.diagonal(vyb,axis1=-2,axis2=-1))/spec.temp
    # row index m, column index n. f_nm = f_n-f_m.
    f_nm=df[...,None,:]-df[...,:,None]
    eps_mn=e[..., :,None]-e[...,None,:]
    eye=np.eye(2,dtype=bool)
    rho=np.where(eye,0.,-zb*f_nm/np.where(eye,1.,eps_mn))
    return complex(np.mean(np.einsum('...mn,...nm->...',rho,vxb)))
