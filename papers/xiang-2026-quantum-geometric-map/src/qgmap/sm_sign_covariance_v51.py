"""v51 scientific sign and spin-field covariance audit for Xiang et al. PRB SM.

Literal primary SM p.5 Eq(20), p.8 Eq(31)-(33), p.13 Eq(62),
main PRB Eq(10)-(11). This is algebra / basis-covariance and finite-T
Dirac-model method audit, NOT a manuscript erratum claim or physical transport.
"""
from __future__ import annotations
import numpy as np
from numpy.polynomial.legendre import leggauss

I=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)


def source_sign_identity(e, v, spin, grad_f):
    """Generic pointwise two-band Eq20 current vs Eq31 curvature.

    e sorted arbitrary nondegenerate energies; v and spin Hermitian *band*
    matrices; grad_f real occupied band derivatives. Per E*Bbar and tau=1.
    Not derived from Dirac or using pre-programmed response sign.
    """
    e=np.asarray(e,float); v=np.asarray(v,complex); spin=np.asarray(spin,complex); grad_f=np.asarray(grad_f,float)
    if e.shape!=(2,) or v.shape!=(2,2) or spin.shape!=(2,2) or grad_f.shape!=(2,):
        raise ValueError('two-band input required')
    if abs(e[0]-e[1])<1e-12 or not np.isfinite(e).all() or not np.isfinite(grad_f).all():
        raise ValueError('non-degenerate finite bands required')
    if not np.allclose(v,v.conj().T) or not np.allclose(spin,spin.conj().T):
        raise ValueError('Hermitian operators required')
    j20=0j; zpart=0j; j31=0j
    for m,n in ((0,1),(1,0)):
        epsmn=e[m]-e[n]
        rho20=-spin[m,n]*(grad_f[n]-grad_f[m])/epsmn
        j20+=rho20*v[n,m]
        # the actual author's velocity relation v_nm=i*epsilon_nm*r_nm
        rnm=v[n,m]/(1j*(e[n]-e[m]))
        rmn=v[m,n]/(1j*(e[m]-e[n]))
        znm=1j*(rnm*spin[m,n]-rmn*spin[n,m])
        zpart+=znm*grad_f[n]
        j31-=znm*grad_f[n]
    return dict(j20_per_Bbar=float(j20.real),Z_df=float(zpart.real),
                printed_eq31_per_Bbar=float(j31.real),
                direct_minus_Z=float((j20-zpart).real),
                eq20_plus_printed_eq31=float((j20+j31).real),
                imaginary_spurious=float(max(abs(j20.imag),abs(zpart.imag))))


def hamiltonian(kx,ky,beta=.125,helicity='main'):
    if helicity not in ('main','sm'): raise ValueError('main or sm')
    sign=1 if helicity=='main' else -1
    return beta*kx*I+sign*(kx*SY-ky*SX)


def spin_operator(spin_sign=1):
    if spin_sign not in (-1,1): raise ValueError('spin sign +/-1 required')
    return spin_sign*SX


def single_point(kx=.29,ky=.41,beta=.125,helicity='main',spin_sign=1,mu=.3,T=.04):
    """Source literal two-current coefficients at a Bloch momentum and finite T."""
    if not np.isfinite([kx,ky,beta,mu,T]).all() or T<=0 or not abs(beta)<1 or np.hypot(kx,ky)<1e-10:
        raise ValueError('valid nondegenerate point required')
    h=hamiltonian(kx,ky,beta,helicity)
    eps,U=np.linalg.eigh(h)
    v_sign=1 if helicity=='main' else -1
    vs={'x':beta*I+v_sign*SY,'y':-v_sign*SX}
    S=U.conj().T@spin_operator(spin_sign)@U
    f=1/(1+np.exp(np.clip((eps-mu)/T,-500,500)))
    fprime=-f*(1-f)/T
    out={}
    for a,b in [('x','y'),('y','x')]:
        va=U.conj().T@vs[a]@U
        vb=U.conj().T@vs[b]@U
        grad=fprime*np.real(np.diag(vb))
        out[a+b]=source_sign_identity(eps,va,S,grad)
    return out


def finite_difference_spin_ward(kx=.29,ky=.41,beta=.125,helicity='main',spin_sign=1,step=5e-6):
    """Source SM Eq33 pointwise completeness: sum_{m!=n} Z_nm^ax = d_a <n|Sx|n>."""
    if step<=0 or not np.isfinite(step):raise ValueError('step > 0')
    hs=1 if helicity=='main' else -1
    H=hamiltonian(kx,ky,beta,helicity)
    e,U=np.linalg.eigh(H)
    S=spin_operator(spin_sign)
    Ud=U.conj().T
    Ss=Ud@S@U
    out={}
    for alpha in ('x','y'):
        vx=beta*I+hs*SY if alpha=='x' else -hs*SX
        vv=Ud@vx@U
        z=[]
        for n,m in ((0,1),(1,0)):
            rnm=vv[n,m]/(1j*(e[n]-e[m]))
            rmn=vv[m,n]/(1j*(e[m]-e[n]))
            z.append(float(np.real(1j*(rnm*Ss[m,n]-rmn*Ss[n,m]))))
        xp=step if alpha=='x' else 0.
        yp=step if alpha=='y' else 0.
        def spinexp(xx,yy):
            _,u=np.linalg.eigh(hamiltonian(xx,yy,beta,helicity))
            return np.real(np.diag(u.conj().T@S@u))
        fd=(spinexp(kx+xp,ky+yp)-spinexp(kx-xp,ky-yp))/(2*step)
        out[alpha]={'Z_band':z,'d_spin_fd':fd.tolist(),
                    'maximum_residual':float(np.max(abs(np.asarray(z)-fd)))}
    return out


def integrated_table(beta=.125,mu=.3,T=.03,nr=80,nphi=112,cutoff=2.5):
    """Separate main/SM helicity and magnetic source sign on same numerical grid.

    Report per Bbar; physical B converts both by the SAME minus sign.
    """
    if not (0<beta<1 and T>0 and cutoff>0 and nr>=20 and nphi>=32):raise ValueError('parameters')
    xx,ww=leggauss(nr)
    r=cutoff*(xx+1)/2; dr=cutoff*ww/2
    th=2*np.pi*(np.arange(nphi)+.317)/nphi
    kx=(r[:,None]*np.cos(th)).ravel()
    ky=(r[:,None]*np.sin(th)).ravel()
    weights=np.broadcast_to((r*dr)[:,None]/(2*np.pi*nphi),(nr,nphi)).ravel()
    table={}
    for hname in ('main','sm'):
        hsign=1 if hname=='main' else -1
        hh=beta*kx[:,None,None]*I+hsign*(kx[:,None,None]*SY-ky[:,None,None]*SX)
        e,u=np.linalg.eigh(hh)
        ud=u.conj().swapaxes(-1,-2)
        vx=beta*I+hsign*SY;vy=-hsign*SX
        vv_x=ud@vx@u; vv_y=ud@vy@u
        f=1/(1+np.exp(np.clip((e-mu)/T,-500,500)))
        fp=-f*(1-f)/T
        for ssign in (1,-1):
            sp=ud@(ssign*SX)@u
            bxy=fp*np.real(np.diagonal(vv_y,axis1=-2,axis2=-1))
            byx=fp*np.real(np.diagonal(vv_x,axis1=-2,axis2=-1))
            J={}
            for key,va,df in (('xy',vv_x,bxy),('yx',vv_y,byx)):
                j20=np.zeros(len(kx),complex)
                j31=np.zeros(len(kx),complex)
                for m,n in ((0,1),(1,0)):
                    de=e[:,m]-e[:,n]
                    rho=-sp[:,m,n]*(df[:,n]-df[:,m])/de
                    j20+=rho*va[:,n,m]
                    rn=va[:,n,m]/(1j*(e[:,n]-e[:,m]))
                    rm=va[:,m,n]/(1j*(e[:,m]-e[:,n]))
                    zn=1j*(rn*sp[:,m,n]-rm*sp[:,n,m])
                    j31-=zn*df[:,n]
                J[key]={'eq20_per_Bbar':float(np.real(np.dot(weights,j20))),
                        'printed_eq31_per_Bbar':float(np.real(np.dot(weights,j31)))}
            J['sym_eq20']=.5*(J['xy']['eq20_per_Bbar']+J['yx']['eq20_per_Bbar'])
            J['sym_printed_eq31']=.5*(J['xy']['printed_eq31_per_Bbar']+J['yx']['printed_eq31_per_Bbar'])
            table[hname+('_spin_x' if ssign==1 else '_negative_spin_x')]=J
    covariant_pairs=[('main_spin_x','sm_negative_spin_x'),('main_negative_spin_x','sm_spin_x')]
    return {'protocol':'v51_exact_spin_basis_covariance_and_SM20_31_relative_sign',
            'scientific_status':'PARTIAL','beta':beta,'mu':mu,'temperature':T,
            'table':table,
            'basis_covariance_residuals':{a+'_vs_'+b:{kind:table[a][kind]-table[b][kind]
               for kind in ('sym_eq20','sym_printed_eq31')} for a,b in covariant_pairs},
            'limitation':'Changing helicity without transforming sigma_x changes the physical model; neither operation repairs the Eq20/31 relative sign.'}


def analyze():
    t=integrated_table()
    t['local_point']=single_point()
    t['eq33_local_completeness']=finite_difference_spin_ward()
    return t
