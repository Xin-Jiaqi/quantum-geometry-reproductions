"""Independent N-band SM(16-19) spin density-matrix reproduction.

This is deliberately distinct from the v4 two-band implementation. The fixed
four-band reference is orbital(2) x physical Pauli spin(2), generically
nondegenerate and periodic. Units: e=hbar=E=1, current divided by E_b Bbar_c;
Bbar=-g mu_B B. Full BZ integration uses periodic quadrature.

This implementation does NOT address degeneracies, collisions, non-Abelian
projectors or real-material Wannier operator definitions.
"""
from __future__ import annotations
import numpy as np
from ._spin_grid import SX, SY, SZ, I, fermi, mesh


def model_four_band(N=61,seed=713,scale=1.,offset=(.31,.23)):
    """Deterministic Hermitian 4x4 Bloch model + exact velocities and Pauli S_z."""
    if N<4 or scale<=0:raise ValueError('positive model scale and N>=4 required')
    x,y=mesh(N,offset)
    nk=x.size
    rng=np.random.default_rng(seed)
    def random_hermitian():
        z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        return (z+z.conj().T)/2
    onsite=np.diag([-2.7,-1.05,.6,2.5]).astype(complex)
    C=[random_hermitian() for _ in range(6)]
    cs=np.stack([np.sin(x),np.cos(x),np.sin(y),np.cos(y),np.sin(x+y),np.cos(x-y)])
    csx=np.stack([np.cos(x),-np.sin(x),np.zeros_like(x),np.zeros_like(x),np.cos(x+y),-np.sin(x-y)])
    csy=np.stack([np.zeros_like(y),np.zeros_like(y),np.cos(y),-np.sin(y),np.cos(x+y),np.sin(x-y)])
    C=np.stack(C)*.095*scale
    H=onsite+np.einsum('pij,pkl->klij',C,cs,optimize=True)
    vx=np.einsum('pij,pkl->klij',C,csx,optimize=True)
    vy=np.einsum('pij,pkl->klij',C,csy,optimize=True)
    S=np.kron(np.eye(2),SZ)
    return H,vx,vy,S


def density_matrix_audit(H,vx,vy,S,mu=.15,T=.29,a=0,b=1,phase_seed=None):
    """Band-wise rho^(EB;1-4) and a separately evaluated SM Eq.(29)."""
    H=np.asarray(H,complex)
    nb=H.shape[-1]
    if H.ndim<3 or H.shape[-2:]!=(nb,nb) or H.shape!=vx.shape or H.shape!=vy.shape:
        raise ValueError('H,vx,vy must have identical batch and square matrix shape')
    if S.shape!=(nb,nb) or T<=0:raise ValueError('constant Hermitian S and T>0 required')
    if not np.allclose(H,H.swapaxes(-1,-2).conj(),atol=1e-10):raise ValueError('H non-Hermitian')
    if not np.allclose(S,S.conj().T,atol=1e-10):raise ValueError('spin non-Hermitian')
    if (a,b) not in ((0,1),(1,0)):raise ValueError('transverse a,b only')
    e,U=np.linalg.eigh(H)
    min_gap=float(np.min(np.diff(e,axis=-1)))
    if min_gap<=1e-6:raise ValueError('degenerate / too close bands unsupported')
    if phase_seed is not None:
        rng=np.random.default_rng(phase_seed)
        U=U*np.exp(1j*rng.uniform(-np.pi,np.pi,e.shape))[...,None,:]
    ud=U.conj().swapaxes(-1,-2)
    v=[ud@v0@U for v0 in (vx,vy)]; spin=ud@S@U
    vdiag=[np.real(np.diagonal(v0,axis1=-2,axis2=-1)) for v0 in v]
    r=[]
    # r_nm=i<v_n m>/(E_m-E_n), n != m
    for va in v:
        rr=np.zeros_like(va)
        for n in range(nb):
            for m in range(nb):
                if n!=m:rr[...,n,m]=1j*va[...,n,m]/(e[...,m]-e[...,n])
        r.append(rr)
    f,fp=fermi(e,mu,T)
    rho=[np.zeros_like(spin) for _ in range(4)]
    # SM Eq 16 diagonal density; local Zeeman metric Q^bc_nm.
    for n in range(nb):
        for ell in range(nb):
            if n==ell:continue
            eps=e[...,n]-e[...,ell]
            q=np.real(r[b][...,n,ell]*spin[...,ell,n])
            rho[0][...,n,n] += -2*(f[...,n]-f[...,ell])*q/eps**2
    # SM Eqs 17-19, with the off-diagonal covariant derivative from Eq 25.
    Ds=1j*(r[b]@spin-spin@r[b])
    for m in range(nb):
        for n in range(nb):
            if m==n:continue
            eps=e[...,m]-e[...,n]
            fnm=f[...,n]-f[...,m]
            f_grad=fp[...,n]*vdiag[b][...,n]-fp[...,m]*vdiag[b][...,m]
            d_gap=vdiag[b][...,m]-vdiag[b][...,n]
            scalar=fnm/eps
            d_scalar=f_grad/eps-fnm*d_gap/eps**2
            dprod=d_scalar*spin[...,m,n]+scalar*Ds[...,m,n]
            term17=np.zeros_like(e[...,0],complex)
            term18=np.zeros_like(term17)
            for ell in range(nb):
                if ell!=n:
                    dln=e[...,ell]-e[...,n]
                    fnl=f[...,n]-f[...,ell]
                    term17+=r[b][...,m,ell]*spin[...,ell,n]*fnl/dln
                    term18+=spin[...,m,ell]*r[b][...,ell,n]*fnl/dln
                if ell!=m:
                    dml=e[...,m]-e[...,ell]
                    flm=f[...,ell]-f[...,m]
                    term17-=flm*spin[...,m,ell]*r[b][...,ell,n]/dml
                    term18-=flm*r[b][...,m,ell]*spin[...,ell,n]/dml
            rho[1][...,m,n]=1j/eps*(dprod-1j*term17)
            rho[2][...,m,n]=term18/eps
            rho[3][...,m,n]=-1j*spin[...,m,n]*f_grad/eps**2
    vals=[np.mean(np.sum(rr*v[a].swapaxes(-1,-2),axis=(-2,-1))) for rr in rho]
    imags=[float(z.imag) for z in vals]
    vals=[float(z.real) for z in vals]
    # main Eq7, independent of the intermediate rho, via two-band local metrics
    target=np.zeros(e.shape[:-1])
    for n in range(nb):
        for m in range(nb):
            if n==m:continue
            eps=e[...,n]-e[...,m]
            qa=np.real(r[a][...,n,m]*spin[...,m,n])
            qb=np.real(r[b][...,n,m]*spin[...,m,n])
            target+=2*(qa*vdiag[b][...,n]-qb*vdiag[a][...,n])*fp[...,n]/eps
    expected=float(np.mean(target))
    return dict(protocol='SM16_19_4band_spin_density_vs_SM29',
        n_band=nb,mu=mu,T=T,min_gap=min_gap,mesh_count=int(e.shape[0]*e.shape[1]),
        j_rho_terms=vals,imag_rho_currents=imags,total_rho=sum(vals),
        direct_SM29=expected,closing_residual=sum(vals)-expected,
        tau_order='tau^0',magnetic_convention='Bbar=-g*mu_B*B (Pauli spin)',
        units='dimensionless e=hbar=1, momentum integral BZ mean')


def four_band_audit(N=61,mu=.15,T=.29,phase_seed=None,seed=713,scale=1.,a=0,b=1):
    H,vx,vy,S=model_four_band(N,seed=seed,scale=scale)
    return density_matrix_audit(H,vx,vy,S,mu,T,a,b,phase_seed)
