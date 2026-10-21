"""v45: original-paper tilted Dirac PHE protocol stress test (not a correction).

Source: Xiang et al. PRB 113 L201406 (2026), main Eqs 7,8,10,11.
The paper's magnetic linear PHE is SYMMETRIC xy=yx, order tau^1.
The preceding v41-v44 QWZ charge Hall is ANTISYMMETRIC and a separate
method benchmark; it does not establish this nonlinear E_y*B_x coefficient.

For the *fixed-tilt*, unbounded continuum Dirac model, a uniform spin
Zeeman Bx enters entirely as ky -> ky + delta. This proves that the FULL
RTA equilibrium-population Drude conductivity is Bx-invariant when integrated
over the entire momentum plane. This protocol is not automatically identical
with the original paper's partial Liouville density-matrix response Eq.(8).

The frozen Eq.(8) spectral-angular and main Eq.(11) branch are separate,
including their sign convention; do not silently fix it or infer figure units.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import quad

I = np.eye(2, dtype=complex)
SX = np.array([[0.,1.],[1.,0.]],complex)
SY = np.array([[0.,-1j],[1j,0.]],complex)


def hamiltonian(kx:float,ky:float,beta=.125,delta=0.,helicity='main'):
    """Fixed beta=tx/vF and vF=1; delta=g mu_B Bx/vF for H_B=-g mu_B Bx SX."""
    if not np.isfinite([kx,ky,beta,delta]).all() or not abs(beta)<1:
        raise ValueError('finite k, |beta|<1 required')
    if helicity not in ('main','supplementary'):
        raise ValueError('select paper main Eq10 or SI Eq62 helicity explicitly')
    sign=1 if helicity=='main' else -1
    return beta*kx*I+sign*(kx*SY-ky*SX)-delta*SX


def shift_identity(kx,ky,beta=.125,delta=.07,helicity='main'):
    """Algebraic identity H(k;delta)=H(kx,ky+sign*delta;0)."""
    sign=1 if helicity=='main' else -1
    lhs=hamiltonian(kx,ky,beta,delta,helicity)
    rhs=hamiltonian(kx,ky+sign*delta,beta,0,helicity)
    return float(np.max(abs(lhs-rhs)))


def zeeman_curvature(kx,ky,band=1,helicity='main'):
    """Direct eigenbasis Z^{xx},Z^{yx} via interband velocity/spin matrices.

    Here band 0=valence and 1=conduction, Zeeman curvature is defined by
    main Eq.(8): Z^{a x}_{nm} = -2 Im[r^a_nm sigma^x_mn],
    r^a_nm=i v^a_nm/(eps_m-eps_n). No finite difference of phases needed.
    """
    if band not in (0,1):raise ValueError('select nondegenerate band')
    if np.hypot(kx,ky)<1e-12: raise ValueError('Dirac degeneracy')
    H=hamiltonian(kx,ky,helicity=helicity)
    sign=1 if helicity=='main' else -1
    e,U=np.linalg.eigh(H)
    Ud=U.conj().T
    sx=Ud@SX@U
    vx=Ud@(sign*SY)@U
    vy=Ud@(-sign*SX)@U
    m=1-band
    zx=-2*np.imag(1j*vx[band,m]*sx[m,band]/(e[m]-e[band]))
    zy=-2*np.imag(1j*vy[band,m]*sx[m,band]/(e[m]-e[band]))
    return float(zx),float(zy)


def printed_eq11_geometry(beta):
    """Main Eq.11 dimensionless geometric factor, no tau, e, hbar, or g muB."""
    if not 0<beta<1:raise ValueError('printed Eq11 assumes 0<beta<1')
    return -beta/(4*np.pi*(1+np.sqrt(1-beta*beta))**2)


def angular_eq8_geometry(beta,mu_sign=1):
    """Independently integrate main Eq.8 using its printed Z and band velocity.

    Normalize by (tau*g*mu_B) with vF=1; B is NOT inserted here.
    Main Eq10: E_s(k,theta)=k(s+beta cos(theta)),
    Zxx_s=s cos(theta)sin(theta)/k, Zyx_s=-s cos^2(theta)/k.
    The negative Fermi derivative fixes an overall minus sign. For |beta|<1,
    mu>0 samples conduction only and mu<0 samples valence only.
    """
    if not 0<beta<1 or mu_sign not in (-1,1):raise ValueError('beta and mu_sign')
    s=int(mu_sign)
    def f(theta):
        c=np.cos(theta); t=np.sin(theta)
        va=beta+s*c; vb=s*t
        Zxx=s*c*t; Zyx=-s*c*c  # radial k cancels k dk measure
        return -(va*Zyx+vb*Zxx)/abs(s+beta*c)
    val,err=quad(f,0,2*np.pi,epsabs=1e-12,epsrel=1e-12,limit=250)
    return float(val/(2*(2*np.pi)**2)), float(err/(2*(2*np.pi)**2))


def equilibrium_drude_xy(beta=.125,mu=.30,thermal=.055,delta=.0,
                         cutoff=2.5,grid=301,helicity='main'):
    """Full Fermi-population constant-tau Drude cross channel on a square grid.

    sigma_xy/tau = int d²k/(2pi)² sum_s v_sx v_sy (-f'(E_s)).
    For fixed-beta unlimited continuum and Bx shift, exactly zero by
    ky+sign*delta -> qy and odd parity in qy. Finite-cutoff/grid only.
    """
    if not np.isfinite([beta,mu,thermal,delta,cutoff]).all() or abs(beta)>=1 or thermal<=0 or cutoff<=0 or not isinstance(grid,int) or grid<21:
        raise ValueError('invalid finite-T Drude integration')
    if helicity not in ('main','supplementary'):raise ValueError('helicity')
    sign=1 if helicity=='main' else -1
    k=np.linspace(-cutoff,cutoff,grid)
    kx,ky=np.meshgrid(k,k,indexing='ij')
    qy=ky+sign*delta
    R=np.hypot(kx,qy)
    Rs=np.where(R<1e-13,1.,R)
    total=np.zeros_like(R)
    for band in (-1,1):
        energy=beta*kx+band*R
        vx=beta+band*kx/Rs
        vy=band*qy/Rs
        x=np.clip((energy-mu)/(2*thermal),-200,200)
        dfermi=1/(4*thermal*np.cosh(x)**2)
        total+=vx*vy*dfermi
    return float(total.sum()*(k[1]-k[0])**2/(2*np.pi)**2)



def drude_delta_derivative_terms(beta=.125,mu=.30,thermal=.055,cutoff=2.5,grid=600):
    """Finite-T derivative at delta=0 separated into velocity and population.

    An *independent* calculus proof of cancellation, from differentiating
    sigma_xy=integral v_x v_y (-f'). The derivative with respect to the
    magnetic spin momentum shift delta equals integral d_{ky}(integrand),
    which vanishes on the full momentum plane. This differs from taking
    only the paper's selected density-matrix geometric contribution.
    """
    if not 0<beta<1 or thermal<=0 or cutoff<=0 or not isinstance(grid,int) or grid<100 or grid%2:
        raise ValueError('Require 0<beta<1, T>0, even grid>=100')
    axis=np.linspace(-cutoff,cutoff,grid)
    kx,ky=np.meshgrid(axis,axis,indexing='ij')
    r=np.hypot(kx,ky)  # even mesh avoids nodal k=0
    part_velocity=0.
    part_occupation=0.
    for s in (-1,1):
        e=beta*kx+s*r
        vx=beta+s*kx/r
        vy=s*ky/r
        dvx=-s*kx*ky/r**3
        dvy=s*kx*kx/r**3
        u=np.clip((e-mu)/(2*thermal),-200,200)
        negfp=1/(4*thermal*np.cosh(u)**2)
        dnegfp=-negfp*np.tanh(u)/thermal
        part_velocity+=float(np.sum((dvx*vy+vx*dvy)*negfp))
        part_occupation+=float(np.sum(vx*vy*dnegfp*vy))
    fac=(axis[1]-axis[0])**2/(2*np.pi)**2
    return {'velocity_change':part_velocity*fac,
            'equilibrium_occupation_change':part_occupation*fac,
            'full_derivative':(part_velocity+part_occupation)*fac}

def audit(beta=.125,mu=.3,temp=.055,delta=.10,grid=301):
    eq8,p_err=angular_eq8_geometry(beta,1 if mu>0 else -1)
    eq11=printed_eq11_geometry(beta)*np.sign(mu)
    p=equilibrium_drude_xy(beta,mu,temp,delta,grid=grid)
    m=equilibrium_drude_xy(beta,mu,temp,-delta,grid=grid)
    zero=equilibrium_drude_xy(beta,mu,temp,0,grid=grid)
    split=drude_delta_derivative_terms(beta,mu,temp,grid=400)
    return {
        'protocol':'fixed_tilt_uniform_Bx_RTA_phe_shift_vs_printed_Eq8_Eq11_v45',
        'scientific_status':'PARTIAL',
        'beta':beta,'mu':mu,'temp_model':temp,'delta_model':delta,
        'main_eq8_angular_coefficient_per_tau_gmuB':eq8,
        'paper_main_eq11_geometry':float(eq11),
        'eq8_plus_eq11_sign_residual':float(eq8+eq11),
        'drude_xy_plus_B_per_tau':p,'drude_xy_minus_B_per_tau':m,
        'drude_xy_zero_B_per_tau':zero,
        'finite_B_central_dsigma_d_delta_per_tau':(p-m)/(2*delta),
        'drude_derivative_split':split,
        'shift_matrix_residual_main':shift_identity(.21,-.31,beta,delta,'main'),
        'shift_matrix_residual_sm':shift_identity(.21,-.31,beta,delta,'supplementary'),
        'distinct_protocols':[
            'Eq8 is the paper density-matrix geom term; a finite-T Drude derivative includes equilibrium field preparation',
            'The fixed-tilt infinite Dirac cone has a pure ky shift under uniform Bx, which enforces zero equilibrium Drude derivative',
            'Finite sample/UV, B-dependent tilt, or baths/collisions can change this, and require explicit declaration',
            'Eq8-versus-Eq11 sign depends on field/helicity conventions; not silently repaired',
            'The previously frozen Fig2(b) ~2pi vertical discrepancy is not resolved'],
        'scientific_open':'which field-preparation and relaxation protocol reproduces the published finite PHE despite the momentum-shift Ward identity'
    }
