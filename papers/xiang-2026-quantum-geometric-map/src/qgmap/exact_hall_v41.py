"""v41: exact finite-torus, zero-T *linear* Hall Kubo benchmark.

This is a full single-particle Lehmann Hall sum for each frozen scalar-disorder
sample, not an isolated RA bubble, a SCBA ladder, a numerical Streda I/II
split, or the original nonlinear E*B PHE.  e=hbar=a=1, sigma in e^2/h.
QWZ H(k)=sin(kx)sx+sin(ky)sy+(m+cos(kx)+cos(ky))sz.
Twist theta is a boundary flux: uniform vector potential A=theta/L.
"""
from __future__ import annotations
import numpy as np

_I=np.eye(2,dtype=complex)
_SX=np.array([[0,1],[1,0]],complex)
_SY=np.array([[0,-1j],[1j,0]],complex)
_SZ=np.diag([1.,-1.]).astype(complex)


def torus(L=4, mass=-1., W=.0, seed=0, theta=(0.,0.)):
    """(H,Jx,Jy) in a FIXED real-site basis, shared disorder for every twist.

    Disorder potential V_r has variance W per site and is the same for both
    orbital components. It is independent of the vector potential.
    """
    if not isinstance(L,int) or not 3<=L<=9 or W<0 or not np.isfinite([mass,W,*theta]).all():
        raise ValueError('invalid torus parameters')
    if min(abs(mass-x) for x in (-2.,0.,2.)) < 1e-10:
        raise ValueError('QWZ gap-closing mass requires explicit treatment')
    n=L*L
    k=-np.pi+2*np.pi*np.arange(L)/L
    px,py=np.meshgrid(k,k,indexing='ij')
    xx=px+theta[0]/L; yy=py+theta[1]/L
    d=(np.sin(xx),np.sin(yy),mass+np.cos(xx)+np.cos(yy))
    dx=(np.cos(xx),np.zeros_like(xx),-np.sin(xx))
    dy=(np.zeros_like(xx),np.cos(yy),-np.sin(yy))
    def blocks(v):
        return v[0][...,None,None]*_SX + v[1][...,None,None]*_SY + v[2][...,None,None]*_SZ
    sites=np.stack(np.meshgrid(np.arange(L),np.arange(L),indexing='ij'),-1).reshape(n,2)
    momenta=np.stack((px,py),-1).reshape(n,2)
    U=np.kron(np.exp(1j*sites@momenta.T)/np.sqrt(n),_I)
    def real_space(mat):
        b=np.zeros((2*n,2*n),complex)
        for i,m in enumerate(mat.reshape(n,2,2)):
            b[2*i:2*i+2,2*i:2*i+2]=m
        return U@b@U.conj().T
    H=real_space(blocks(d))
    if W>0:
        rng=np.random.default_rng(seed)
        H+=np.diag(np.repeat(rng.normal(0,np.sqrt(W),size=n),2))
    return H,real_space(blocks(dx)),real_space(blocks(dy))


def hall_lehmann(H,Jx,Jy,mu=0.,omega=0.,eta=0.):
    """T=0, full occupied-to-empty spectral sum of antisymmetric sigma_xy.

    Hall orientation matches v38: Omega_occ = 2 Im[Jx_nm Jy_mn]/Delta^2.
    Complex frequency omega+i*eta regularizes RESPONSE, not one-particle G.
    For an energy gap at mu this is an exact finite-volume single-particle
    transverse Kubo response; sigma(0,0) is a Berry-curvature sum.
    """
    H=np.asarray(H); Jx=np.asarray(Jx); Jy=np.asarray(Jy)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or Jx.shape!=H.shape or Jy.shape!=H.shape or H.shape[0]%2:
        raise ValueError('matching even-dimensional matrices required')
    if eta<0 or not np.isfinite([mu,omega,eta]).all():raise ValueError('invalid response arguments')
    for x in (H,Jx,Jy):
        if not np.allclose(x,x.conj().T,atol=1e-10):raise ValueError('Hermitian matrices required')
    eig,U=np.linalg.eigh(H)
    occ=eig<mu
    if not occ.any() or occ.all():raise ValueError('occupied and empty states required')
    gap=float(eig[~occ].min()-eig[occ].max())
    Vx=U.conj().T@Jx@U; Vy=U.conj().T@Jy@U
    prod=Vx[np.ix_(occ,~occ)]*Vy.T[np.ix_(occ,~occ)]
    delta=eig[~occ][None,:]-eig[occ][:,None]
    z=complex(omega,eta)
    if omega and eta==0 and np.min(abs(delta-abs(omega)))<1e-10:
        raise ValueError('undamped interband pole: set eta>0')
    # sigma/(e^2/h) = 2pi/L^2 * sum_occ,empty 2 Im(P)/(Delta^2-z^2)
    L=int(round(np.sqrt(H.shape[0]/2)))
    if 2*L*L!=H.shape[0]:raise ValueError('requires LxL two-orbital torus')
    val=2*np.pi/(L*L)*np.sum(2*np.imag(prod)/(delta**2-z**2))
    return {'hall_real':float(val.real),'hall_imag':float(val.imag),
            'occupied':int(occ.sum()),'minimum_fermi_gap':gap,
            'minimum_transition':float(delta.min())}


def hall_retarded_crosscheck(H,Jx,Jy,mu=0.,omega=.1,eta=.03):
    """Independent Lehmann chi_xy and chi_yx, then conductivity -chi_H/(i*z).

    This deliberately sums ALL occupied/empty ordered state pairs rather than
    directly evaluating the simplified Hall denominator. The mixed A_x A_y
    diamagnetic contact vanishes for QWZ because d_x d_y H/dA_x dA_y=0.
    """
    z=complex(omega,eta)
    if abs(z)==0:raise ValueError('choose finite omega or eta for Kubo quotient')
    H=np.asarray(H)
    eig,U=np.linalg.eigh(H)
    f=(eig<mu).astype(float)
    if f.sum() in (0,len(f)):raise ValueError('needs occupied/empty states')
    xx=U.conj().T@Jx@U
    yy=U.conj().T@Jy@U
    ff=f[:,None]-f[None,:]
    den=z+eig[:,None]-eig[None,:]
    chi_xy=np.sum(ff*xx*yy.T/den)
    chi_yx=np.sum(ff*yy*xx.T/den)
    L=int(round(np.sqrt(len(eig)/2)))
    if 2*L*L!=len(eig):raise ValueError('requires two-orbital square torus')
    sigma=-2*np.pi/(L*L)*(chi_xy-chi_yx)/(2j*z)
    return complex(sigma)


def _seam_gauge(L,which):
    coords=np.stack(np.meshgrid(np.arange(L),np.arange(L),indexing='ij'),axis=-1).reshape(-1,2)
    return np.repeat(np.exp(-2j*np.pi*coords[:,which]/L),2)


def twist_chern(L=4,mass=-1.,W=0.,seed=0,ntheta=5,mu=0.):
    """Independent occupied-subspace FHS determinant, including flux seam gauges.

    The same disorder realization and a fixed site basis are used at every
    twist. Occupied rank must not change. Requires a gap on the twist grid.
    """
    if ntheta<3:raise ValueError('ntheta>=3')
    eigens={}; nocc=None; gap=np.inf; idxs=np.arange(ntheta)
    for ix in idxs:
        for iy in idxs:
            theta=(2*np.pi*ix/ntheta,2*np.pi*iy/ntheta)
            H,_,_=torus(L,mass,W,seed,theta)
            e,u=np.linalg.eigh(H)
            occ=e<mu
            if nocc is None:nocc=int(occ.sum())
            if nocc!=int(occ.sum()) or nocc in (0,len(e)):
                raise ValueError('occupied rank changes or is trivial: gap not certified')
            gap=min(gap,float(e[~occ].min()-e[occ].max()))
            eigens[ix,iy]=u[:,occ]
    links=np.empty((ntheta,ntheta,2),complex)
    for ix in idxs:
        for iy in idxs:
            state=eigens[ix,iy]
            for axis in (0,1):
                jx=(ix+int(axis==0))%ntheta
                jy=(iy+int(axis==1))%ntheta
                neighbor=eigens[jx,jy]
                if (ix if axis==0 else iy)==ntheta-1:
                    neighbor=_seam_gauge(L,axis)[:,None]*neighbor
                overlap=state.conj().T@neighbor
                det=np.linalg.slogdet(overlap)
                if det[0]==0 or det[1]<-40:
                    raise ValueError('singular occupied overlap: refine twist grid')
                links[ix,iy,axis]=det[0]
    ux=links[:,:,0];uy=links[:,:,1]
    plaquette=ux*np.roll(uy,-1,axis=0)/(
        np.roll(ux,-1,axis=1)*uy)
    C=float(np.sum(np.angle(plaquette))/(2*np.pi))
    return {'chern':C,'nearest_integer':int(round(C)),
            'minimum_fermi_gap_twist_grid':float(gap),'occupied':nocc}


def audit(L=4,mass=-1.,W=0.0,seed=0,ntheta=5,mu=0.,eta=0.):
    """Fixed-flux exact Hall, twist-average exact Hall and independent FHS C."""
    H,Jx,Jy=torus(L,mass,W,seed)
    zero=hall_lehmann(H,Jx,Jy,mu,eta=eta)
    vals=[]
    for ix in range(ntheta):
        for iy in range(ntheta):
            th=(2*np.pi*ix/ntheta,2*np.pi*iy/ntheta)
            a,b,c=torus(L,mass,W,seed,th)
            vals.append(hall_lehmann(a,b,c,mu,eta=eta)['hall_real'])
    topo=twist_chern(L,mass,W,seed,ntheta,mu)
    return {'protocol':'exact_finite_torus_linear_Hall_Lehmann_twist_FHS',
            'scientific_status':'PARTIAL','L':L,'mass':mass,'W':W,
            'seed':seed,'ntheta':ntheta,'mu':mu,'eta_response':eta,
            'sigma_units':'e^2/h','hall_zero_twist':zero['hall_real'],
            'hall_twist_average':float(np.mean(vals)),
            'chern_fhs':topo['chern'],
            'twist_ave_minus_chern':float(np.mean(vals)-topo['chern']),
            'minimum_fermi_gap_twist_grid':topo['minimum_fermi_gap_twist_grid'],
            'limitations':['exact finite L and finite sample, not thermodynamic disorder average',
                           'this is a complete exact linear Lehmann reference, not Streda I/II decomposition',
                           'SCBA self-energy/current vertices/contact protocol not yet matched',
                           'no nonlinear EB PHE evaluated']}
