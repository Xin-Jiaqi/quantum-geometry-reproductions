"""Clean-limit intrinsic E*B spin response for finite-dimensional Bloch Hamiltonians.

Conventions: H in eV, k in 1/nm, v_beta = dH/dk_beta in eV nm,
spin matrices equal Pauli sigma (not hbar/2),
Zeeman perturbation H_B = -g mu_B B_gamma sigma_gamma,
response (spin magnetic moment / mu_B)/area per (E[V/nm] B[T])
reported numerically in (nm V T)^-1 for a 2D model.

The two-band + three-band formula represents SM Eqs (13)+(14),
while FULL_STATIC adds the disputed diagonal contribution of SM Eq (9).
Not a finite-scattering or degenerate-subspace theory.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Mapping
import numpy as np

MU_B = 5.7883818060e-5  # eV/T
KB = 8.617333262145e-5  # eV/K


@dataclass(frozen=True)
class BlochModel:
    hamiltonian: Callable[[np.ndarray], np.ndarray]  # (Nk,d) -> (Nk,N,N)
    velocity: Callable[[np.ndarray, int], np.ndarray] # dH/dk_beta
    spin: Mapping[str, Callable[[np.ndarray], np.ndarray]] # spin component, k basis
    reciprocal: np.ndarray # rows reciprocal vectors (1/nm)
    g_factor: float = 10.
    name: str = 'custom'
    velocity_spin: Mapping[str, Callable[[np.ndarray, int], np.ndarray]] | None = None


def mesh(reciprocal: np.ndarray, n: int, offset: tuple[float,...] | None=None):
    g=np.asarray(reciprocal, float)
    d=g.shape[0]
    if g.shape!=(d,d) or not 1 <= d <= 3: raise ValueError('reciprocal must be square, 1D-3D')
    if n<2: raise ValueError('n must be at least 2')
    if offset is None: offset=(0.37, 0.23, 0.41)[:d]
    if len(offset)!=d: raise ValueError('offset wrong dimension')
    axes=[(np.arange(n)+o)/n for o in offset]
    frac=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,d)
    return frac@g


def fermi(energy, chemical_potential, temp_K=0.):
    energy=np.asarray(energy,float)
    if temp_K<0: raise ValueError('negative temperature')
    if temp_K==0:return (energy<chemical_potential).astype(float)
    x=np.clip((energy-chemical_potential)/(KB*temp_K),-500,500)
    return 1/(1+np.exp(x))


def _validate_matrices(H, V, S):
    if H.ndim!=3 or H.shape[1]!=H.shape[2]: raise ValueError('H shape invalid')
    if V.shape!=H.shape or any(s.shape!=H.shape for s in S):raise ValueError('H,V,spin shape mismatch')
    for name, a in [('H',H),('velocity',V), *[(f'spin[{i}]',s) for i,s in enumerate(S)]]:
        if not np.allclose(a,a.conj().swapaxes(1,2),atol=1e-10,rtol=1e-10):
            raise ValueError(f'{name} not Hermitian')


def _bands(H,V,Sa,Sc, degeneracy_tol=1e-10):
    _validate_matrices(H,V,[Sa,Sc]); e,U=np.linalg.eigh(H)
    # The printed band-sum formulas require nondegenerate bands. Do not mask physical divergences.
    smallest=float(np.min(np.diff(e,axis=1)))
    if smallest <degeneracy_tol:
        raise ValueError(f'Nondegenerate-band expansion unsafe: minimum adjacent splitting {smallest:.3e} eV. '
                         'Use a symmetry-adapted degenerate-subspace/projector formalism.')
    ud=U.conj().swapaxes(1,2)
    ve=ud@V@U; a=ud@Sa@U; c=ud@Sc@U
    return e,ve,a,c,smallest


def _kernels_from_eigen(e,ve,a,c):
    """Dimensionless band kernels [n,k] for printed 2-band, 3-band, diagonal.

    Matches previous independent SM(13),(14) reproductions; no physical prefactor yet.
    """
    nk, nb=e.shape
    de=e[:,:,None]-e[:,None,:]; eye=np.eye(nb,dtype=bool)[None,:,:]
    den=np.where(eye,1.,de)
    r=np.where(eye,0.,-1j*ve/den)
    vdiag=np.real(np.diagonal(ve,axis1=1,axis2=2))
    vbar=vdiag[:,:,None]-vdiag[:,None,:]
    sdiag_c=np.real(np.diagonal(c,axis1=1,axis2=2)); sdiag_a=np.real(np.diagonal(a,axis1=1,axis2=2))
    Q_a=np.real(r*np.swapaxes(a,1,2));Q_c=np.real(r*np.swapaxes(c,1,2))
    F=-2*np.imag(a*np.swapaxes(c,1,2))
    two=np.where(eye,0.,4*(sdiag_c[:,:,None]-sdiag_c[:,None,:])*Q_a/(den**2)
                 +vbar*F/den**3).sum(axis=2)
    diagonal=np.where(eye,0.,2*(sdiag_a[:,:,None]-sdiag_a[:,None,:])*Q_c/(den**2)).sum(axis=2)
    three=np.zeros((nk,nb),float)
    for n in range(nb):
        for m in range(nb):
            if n==m:continue
            for ell in range(nb):
                if ell==n or ell==m:continue
                enm=de[:,n,m];enl=de[:,n,ell];eml=de[:,m,ell]
                aa=a[:,n,m]*(c[:,m,ell]*r[:,ell,n]-r[:,m,ell]*c[:,ell,n])
                bb=a[:,n,m]*(c[:,m,ell]*r[:,ell,n]+r[:,m,ell]*c[:,ell,n])
                cc=a[:,m,ell]*(c[:,n,m]*r[:,ell,n]+r[:,n,m]*c[:,ell,n])
                three[:,n]+=np.real(-2*aa/enm**2-2*bb/(enm*enl)+2*cc/(enl*eml))
    return {'two_band':two,'three_band':three,'diagonal_disputed':diagonal,
            'printed':two+three,'full_static':two+three+diagonal}


def intrinsic_kernel(model: BlochModel, points: np.ndarray, alpha='x', gamma='z', beta=0,
                     degeneracy_tol=1e-10):
    points=np.asarray(points,float)
    H=np.asarray(model.hamiltonian(points),complex)
    V=np.asarray(model.velocity(points,beta),complex)
    Sa=np.asarray(model.spin[alpha](points),complex)
    Sc=np.asarray(model.spin[gamma](points),complex)
    for label, op in [('alpha',Sa),('gamma',Sc)]:
        if not np.allclose(op,op[:1],atol=1e-11,rtol=1e-11):
            raise ValueError(f'k-dependent {label} spin operator not supported by SM(13)/(14) band kernels; include operator covariant derivatives in an extended theory')
    e,ve,a,c,minsep=_bands(H,V,Sa,Sc,degeneracy_tol)
    out=_kernels_from_eigen(e,ve,a,c)
    return e,out,minsep


def integrate(model: BlochModel, points:np.ndarray, *, alpha='x',gamma='z',beta=0,
              mu=0.,temperature_K=0., degeneracy_tol=1e-10):
    if len(points)==0:raise ValueError('no k points')
    e,out,minsep=intrinsic_kernel(model,points,alpha,gamma,beta,degeneracy_tol)
    f=fermi(e,mu,temperature_K)
    g=np.asarray(model.reciprocal)
    if g.shape[0]!=points.shape[1]:raise ValueError('k dimension mismatch')
    # d^d k/(2pi)^d; for 2D output (nm T V)^(-1) with Pauli magnetic moment unit.
    pref=model.g_factor*MU_B * abs(float(np.linalg.det(g)))/(2*np.pi)**g.shape[0]/len(points)
    coeff={k:float(pref*np.sum(f*v)) for k,v in out.items()}
    coeff['prefactor']=float(pref)
    coeff['min_adjacent_gap_eV']=minsep
    coeff['occupation_min']=float(np.min(f));coeff['occupation_max']=float(np.max(f))
    coeff['charge_gap_eV']=float(np.min(e[:,e.shape[1]//2])-np.max(e[:,e.shape[1]//2-1])) if e.shape[1]%2==0 else float('nan')
    return coeff


def chi_k(energy, U, velocity, obs, frozen_f, *, degeneracy_tol=1e-10):
    """Static spin-electric interband susceptibility at individual k, frozen occupations.
    Without the spin magnetic moment factor mu_B and electric charge (eV/V convention).
    """
    ud=U.conj().swapaxes(1,2)
    se=ud@obs@U; v=ud@velocity@U
    d=energy[:,:,None]-energy[:,None,:]; eye=np.eye(energy.shape[1],dtype=bool)[None]
    den=np.where(eye,1.,d)
    # For coincident eigenvalues with identical f: zero exactly by projectors.
    near=(abs(d)<degeneracy_tol)&(~eye)
    occdiff=frozen_f[:,:,None]-frozen_f[:,None,:]
    if np.any(near & (abs(occdiff)>1e-9)):
        raise ValueError('Degenerate occupied/unoccupied pair: static susceptibility undefined without regularization')
    den=np.where(near,1.,den)
    r=np.where(eye|near,0.,-1j*v/den)
    return np.sum(occdiff*np.real(se*np.swapaxes(r,1,2))/den,axis=(1,2))


def finite_field_reference(model:BlochModel, points:np.ndarray, *,alpha='x',gamma='z',beta=0,
                           mu=0.,temperature_K=0.,B_step_T=0.02, occupancy='frozen'):
    """Independent finite-field derivative of static linear E spin susceptibility.

    frozen: fixed H0 occupations, matches the clean intrinsic interband response;
    equilibrium: re-evaluate Fermi occupations IN THE INTERBAND FORMULA ONLY.
    This is a diagnostic, NOT a complete thermodynamic equilibrium response
    (missing explicit intraband f' susceptibility).
    Requires constant-in-k gamma Zeeman spin matrix unless velocity_spin[gamma] supplied.
    """
    if occupancy not in ('frozen','equilibrium'):raise ValueError('unknown occupancy')
    if B_step_T<=0:raise ValueError('B_step_T must be positive')
    k=np.asarray(points,float);H0=np.asarray(model.hamiltonian(k),complex)
    v0=np.asarray(model.velocity(k,beta),complex)
    Oa=np.asarray(model.spin[alpha](k),complex);Og=np.asarray(model.spin[gamma](k),complex)
    if model.velocity_spin is not None and gamma in model.velocity_spin:
        vog=np.asarray(model.velocity_spin[gamma](k,beta),complex)
    else:
        # An orbital spin operator depends on k only if user explicitly supplies derivative.
        if not np.allclose(Og,Og[:1],atol=1e-12,rtol=1e-12):
            raise ValueError('k-dependent Zeeman operator requires velocity_spin[gamma]')
        vog=0.
    e0=np.linalg.eigvalsh(H0);f0=fermi(e0,mu,temperature_K)
    measure=abs(float(np.linalg.det(np.asarray(model.reciprocal))))/(2*np.pi)**k.shape[1]/len(k)
    vals=[]
    for sign in (+1,-1):
        B=sign*B_step_T
        h=H0-model.g_factor*MU_B*B*Og
        v=v0-model.g_factor*MU_B*B*vog
        e,U=np.linalg.eigh(h)
        f=f0 if occupancy=='frozen' else fermi(e,mu,temperature_K)
        val=chi_k(e,U,v,Oa,f)
        vals.append(float(measure*val.sum()))
    return {'response':(vals[0]-vals[1])/(2*B_step_T),
            'chi_plus':vals[0],'chi_minus':vals[1], 'occupation':occupancy}


def layer_spin(spin_op:np.ndarray, projector:np.ndarray):
    """Hermitian layer-local spin: 1/2(P_l S + S P_l).
    Both must be in exactly the same Bloch/Wannier basis.
    """
    return .5*(projector@spin_op+spin_op@projector)


def integrate_chunked(model:BlochModel,points:np.ndarray,*,alpha='x',gamma='z',beta=0,
                      mu=0.,temperature_K=0.,chunk_size=4096,degeneracy_tol=1e-10):
    """Same physical integral as integrate(), streaming batches to bound RAM.
    Caller must provide a uniform k grid or a mathematically valid equal-weight quadrature.
    """
    points=np.asarray(points,float)
    if chunk_size<=0 or len(points)==0:raise ValueError('bad chunk size or empty mesh')
    totals={name:0. for name in ('two_band','three_band','diagonal_disputed','printed','full_static')}
    minimum=np.inf
    vtop=-np.inf;cbottom=np.inf
    for i in range(0,len(points),chunk_size):
        sub=points[i:i+chunk_size]
        e,terms,minsep=intrinsic_kernel(model,sub,alpha,gamma,beta,degeneracy_tol)
        minimum=min(minimum,minsep)
        if e.shape[1]%2==0:
            mid=e.shape[1]//2
            vtop=max(vtop,float(e[:,mid-1].max()))
            cbottom=min(cbottom,float(e[:,mid].min()))
        f=fermi(e,mu,temperature_K)
        for name,t in terms.items():totals[name]+=float(np.sum(f*t))
    pref=model.g_factor*MU_B*abs(float(np.linalg.det(model.reciprocal)))/(2*np.pi)**points.shape[1]/len(points)
    vals={name:total*pref for name,total in totals.items()}
    vals['min_adjacent_gap_eV']=minimum
    vals['charge_gap_eV']=cbottom-vtop if cbottom!=np.inf else float('nan')
    return vals
