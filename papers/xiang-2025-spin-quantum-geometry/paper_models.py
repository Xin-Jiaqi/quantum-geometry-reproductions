"""Independent Xiang/Jin/Wang PRL 135, 146303 (2025) paper-model equations.

Curated offline V27 subset. Printed/source-integral/static branches remain
separate. Natural e=hbar=1 unless parameters explicitly labelled.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)
I=np.eye(2,dtype=complex)

def tilted_dirac(kx:float,ky:float,v:float=1.,tilt:float=.2,delta:float=.04):
    return tilt*kx*I+v*(kx*SX+ky*SY)+delta*SZ

def tilted_bands(kx:float,ky:float=0.,v:float=1.,tilt:float=.2,delta:float=.04):
    return np.linalg.eigvalsh(tilted_dirac(kx,ky,v,tilt,delta))

def spin_qgt(kx:float,ky:float,v:float=1.,tilt:float=.2,delta:float=.04,
             spin_axis:int=1):
    """Interband spin quantum geometry (real(T), -2*imag(T)), spin j={vx,s}/4."""
    energy,U=np.linalg.eigh(tilted_dirac(kx,ky,v,tilt,delta))
    gap=energy[1]-energy[0]
    if gap<1e-10:raise ValueError('Degenerate: use subspace projectors')
    vx=tilt*I+v*SX;vy=v*SY
    s=(SX,SY,SZ)[spin_axis]
    j=(vx@s+s@vx)/4
    T=(U.conj().T@j@U)[0,1]*(U.conj().T@vy@U)[1,0]/gap**2
    return float(T.real),float(-2*T.imag)

def sm_metric(kx,ky,v=1.,tilt=.2,delta=.04):
    r2=v*v*(kx*kx+ky*ky)+delta*delta
    return tilt*v*(v*v*kx*kx+delta*delta)/(8*r2*r2)

def sm_curvature(kx,ky,v=1.,tilt=.2,delta=.04):
    r2=v*v*(kx*kx+ky*ky)+delta*delta
    return tilt*v*delta/(4*r2**1.5)

def _freq_check(omega,delta):
    if not 0<omega<2*abs(delta):raise ValueError('subgap: 0 < omega < 2|Delta|')

def sm_metric_printed(omega,delta=40.,ratio=.2):
    """Published SM Eq.(56), sigma/(e/8pi), NOT independently corrected."""
    _freq_check(omega,delta)
    D=abs(delta);L=math.log((2*D+omega)/(2*D-omega));a=(1-math.pi)/math.pi
    return ratio*((.25+a*(D/omega)**2)*L-a*D/omega)

def sm_metric_from_eq55(omega,delta=40.,ratio=.2):
    """Exact 2D angular/radial integral of the integrand SM Eq.(55)."""
    _freq_check(omega,delta)
    D=abs(delta);L=math.log((2*D+omega)/(2*D-omega))
    return ratio*((.25+(D/omega)**2)*L-D/omega)

def sm_berry_printed(omega,delta=40.,ratio=.2):
    _freq_check(omega,delta)
    D=abs(delta)
    return ratio*((D/omega)*math.log((2*D+omega)/(2*D-omega))-1)

def sm_metric_independent_2d(omega,delta=40.,ratio=.2,nr=64,ntheta=48):
    """Full polar Gauss integration BEFORE simplifying the Eq.(55) angle."""
    _freq_check(omega,delta)
    D=abs(delta)
    x,w=leggauss(nr);t=(x+1)/2
    q=D*t/(1-t);jac=D/(2*(1-t)**2)
    theta=(np.arange(ntheta)+.5)*2*np.pi/ntheta
    E=np.hypot(q,D)
    kernel=(q[:,None]**2*np.cos(theta[None,:])**2+D**2)/(
        E[:,None]**3*(4*E[:,None]**2-omega**2))
    return float(ratio*omega/math.pi*np.dot(w*jac*q,kernel.mean(axis=1)*2*math.pi))

def main_fig1b(mu,delta=40.,tau_ps=.5):
    """Six distinct source protocols. mu, delta in meV; common mu_B e^2 omitted."""
    D=abs(delta);u=abs(mu)
    if D==0:raise ValueError('gapless')
    tau=tau_ps/.6582119569
    drude=0. if u<=D else -tau*tau*delta*(u*u-D*D)/(16*math.pi*u**3)
    if u<=D:
        main=(3-2*math.pi)/(960*math.pi**2*D**2)
        sm70=3*(1-2*math.pi)/(960*math.pi**2*D**2)
        literal=delta/(240*math.pi*D**3)
    else:
        sm70=((3*(1+math.pi)*delta**3-5*math.pi*delta*u*u)/(
            960*math.pi**2*u**5)+delta*(u*u-D*D)/(32*math.pi*u**5))
        main=sm70
        literal=delta*(25*u*u-21*D*D)/(960*math.pi*u**5)
    U=max(u,D)
    source=delta*(7*U**2-3*D**2)/(192*math.pi*U**5)
    return {'mu_meV':mu,'drude_SM67_natural':drude,
            'intrinsic_maintext_natural':main,'intrinsic_SM70_natural':sm70,
            'intrinsic_literal_SM_integrand_natural':literal,
            'independent_density_matrix_static_natural':source,
            'protocol_note':'static density matrix != uniquely specified driven DC'}

def fig1a_hamiltonian(kx,ky,valley,vx,vy,alpha,beta):
    return vx*kx*SY-valley*vy*ky*SX+beta*SZ+valley*alpha*ky*I

def fig1a_spin_metric(kx,ky,valley,vx,vy,alpha,beta):
    h2=(vx*kx)**2+(vy*ky)**2+beta**2
    if h2==0:raise ValueError('degenerate')
    return valley*vx*alpha*((vy*ky)**2+beta**2)/(8*h2*h2)

def fig1a_matrix_metric(kx,ky,valley,vx,vy,alpha,beta):
    energy,U=np.linalg.eigh(fig1a_hamiltonian(kx,ky,valley,vx,vy,alpha,beta))
    vxop=vx*SY;vyop=-valley*vy*SX+valley*alpha*I
    j=(vyop@SY+SY@vyop)/4
    jb=U.conj().T@j@U;vb=U.conj().T@vxop@U
    return float(np.real(jb[0,1]*vb[1,0])/(energy[1]-energy[0])**2)

def fig1a_ndsc(mu,vx,vy,alpha,beta,tau=1.,nodes=128):
    """Eq.(9) two-valley NDSC; exact zero-temperature elliptical Fermi shell."""
    t=alpha/vy;l=1-t*t
    if l<=0:raise ValueError('type-II cone unsupported')
    if abs(mu)<=abs(beta)*math.sqrt(l):return 0.
    band=1 if mu>0 else -1
    radius=math.sqrt(mu*mu/l-beta*beta)
    theta=2*math.pi*(np.arange(nodes)+.5)/nodes
    qx=radius*np.cos(theta);total=0.
    for valley in (-1,1):
        qy=radius*np.sin(theta)/math.sqrt(l)-valley*t*mu/l
        h=np.sqrt(qx*qx+qy*qy+beta*beta)
        weights=h/(vx*vy*math.sqrt(l))/(2*math.pi*nodes)
        velocity=valley*alpha+band*vy*qy/h
        g=valley*vx*alpha*(qy*qy+beta*beta)/(8*h**4)
        total+=float(np.sum(weights*(2*band*h)*g*velocity))
    return tau*tau*total

def reproduce():
    """Machine-readable model subset, NOT publisher vector-curve extraction."""
    vx=.6582119569*1e5*1e-3;alpha=.1*vx;tau=.5/.6582119569
    return {'paper':'PRL 135 146303 (2025)','status':'PARTIAL',
            'fig1a':{'v_hbar_meV_nm':vx,'ndsc_mu20':fig1a_ndsc(20,vx,vx,alpha,10.,tau)},
            'fig1b':[main_fig1b(mu) for mu in (0.,60.,80.)],
            'sm_fig1b':[{'omega_meV':w,'printed_eq56':sm_metric_printed(w),
                        'integrated_eq55':sm_metric_independent_2d(w),
                        'closed_eq55':sm_metric_from_eq55(w),
                        'printed_eq58':sm_berry_printed(w)} for w in (20.,40.,60.)],
            'scientific_gates':{'fig1a_absolute':'OPEN_FACTOR_TWO',
                'fig1b_published_formula':'PASS',
                'SM_eq55_vs_eq56':'OPEN_SOURCE_INCONSISTENCY',
                'metal_DC_protocol':'NOT_VALIDATED','DFT':'OUT_OF_SCOPE'}}
