"""v65: static EQUILIBRIUM orbital electric/magnetic Maxwell derivative, not Hall.

Explicit finite open Peierls graph, a SUPPLIED full finite projected polarization
P=q diag(position)+Hermitian off-diagonal dipoles, T>0. The complete mixed
grand-potential derivative includes d_B P in field-dependent representations.
Fixed-N canonical response includes mu_B so that P+cI origin shifts cancel.
Not a driven EB conductivity or a full realistic Wannier embedding.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from .embedding_magnetic_v55 import peierls_field, field_gauge_unitary


def _check_herm(x,label):
    x=np.asarray(x,complex)
    if x.ndim!=2 or x.shape[0]!=x.shape[1] or not np.isfinite(x).all() or not np.allclose(x,x.conj().T,atol=1e-11,rtol=0):
        raise ValueError(label+' must be finite Hermitian matrix')
    return x


def fermi(e,mu,temp):
    return np.exp(-np.logaddexp(0,(np.asarray(e,dtype=float)-mu)/temp))


def finite_fermi_frechet(H,M,mu,temp):
    """Spectral Fermi divided differences, including exactly degenerate bands."""
    H=_check_herm(H,'H');M=_check_herm(M,'M')
    if H.shape!=M.shape or not np.isfinite([mu,temp]).all() or temp<=0:
        raise ValueError('matched Hermitian operators, T>0 and finite mu')
    e,u=np.linalg.eigh(H);f=fermi(e,mu,temp)
    de=e[:,None]-e[None,:];df=f[:,None]-f[None,:]
    fm=fermi((e[:,None]+e[None,:])/2,mu,temp)
    deriv=-fm*(1-fm)/temp
    tol=1e-10*max(1.,np.max(abs(e)))
    dd=np.divide(df,de,out=np.zeros_like(de),where=abs(de)>tol)
    dd=np.where(abs(de)>tol,dd,deriv)
    return u@(dd*(u.conj().T@M@u))@u.conj().T


@dataclass(frozen=True)
class Spec:
    mu:float=.13
    temperature:float=.37
    field:float=.23
    charge_over_hbar:float=-.84
    def __post_init__(self):
        if not np.isfinite([self.mu,self.temperature,self.field,self.charge_over_hbar]).all() or self.temperature<=0 or self.charge_over_hbar==0:
            raise ValueError('finite mu,B,q and positive T')


def full_projected_polarization(positions,charge,axis='y',dipoles=None):
    xy=np.asarray(positions,float)
    if xy.ndim!=2 or xy.shape[1]!=2 or len(xy)<2 or not np.isfinite(xy).all() or not np.isfinite(charge) or charge==0 or axis not in ('x','y'):
        raise ValueError('finite (N,2) positions and signed charge')
    p=np.diag(charge*xy[:,0 if axis=='x' else 1]).astype(complex)
    if dipoles is not None:
        d=_check_herm(dipoles,'projected position dipoles')
        if d.shape!=p.shape:raise ValueError('dipole dimensions')
        p+=d
    return p


def operators(hoppings,positions,P,spec=Spec(),frame='symmetric',K=None):
    """H,M=d_BH,P,dP=d_B P, including B-dependent magnetic/orbital frame terms."""
    H0=_check_herm(hoppings,'hoppings');P=_check_herm(P,'projected polarization')
    xy=np.asarray(positions,float)
    if P.shape!=H0.shape or xy.shape!=(len(H0),2):raise ValueError('dimensions')
    s=peierls_field(H0,xy,spec.field,spec.charge_over_hbar,gauge='symmetric')
    h,m=s.hamiltonian,s.magnetic_vertex
    if frame=='symmetric':
        U=np.eye(len(h),dtype=complex);dU=np.zeros_like(U)
    elif frame=='landau_y':
        U,dU=field_gauge_unitary(xy,spec.field,spec.charge_over_hbar)
    elif frame=='orbital':
        K=_check_herm(K,'generator K')
        if K.shape!=h.shape:raise ValueError('generator dimension')
        U=expm(1j*spec.field*K);dU=1j*K@U
    else:raise ValueError('frame symmetric, landau_y or orbital')
    Ud=U.conj().T;dUd=dU.conj().T
    H=U@h@Ud
    M=U@m@Ud+dU@h@Ud+U@h@dUd
    Pb=U@P@Ud
    dPb=dU@P@Ud+U@P@dUd
    for mat,name in ((H,'H'),(M,'M'),(Pb,'P'),(dPb,'dP')):_check_herm(mat,name)
    return {'H':H,'M':M,'P':Pb,'dP':dPb,'U':U,'dU':dU,'frame':frame}


def grand_potential(h,mu,temp):
    e=np.linalg.eigvalsh(_check_herm(h,'h'))
    if not np.isfinite([mu,temp]).all() or temp<=0:raise ValueError('T>0')
    return float(-temp*np.sum(np.logaddexp(0,-(e-mu)/temp)))


def equilibrium_derivatives(hoppings,positions,P,spec=Spec(),frame='symmetric',K=None,E=0.):
    if not np.isfinite(E):raise ValueError('finite E')
    op=operators(hoppings,positions,P,spec,frame,K)
    H=op['H']-E*op['P'];M=op['M']-E*op['dP']
    e,u=np.linalg.eigh(H);f=fermi(e,spec.mu,spec.temperature)
    rho=(u*f[None,:])@u.conj().T
    return {'F':grand_potential(H,spec.mu,spec.temperature),
            'dF_dE':float(np.trace(rho@(-op['P'])).real),
            'dF_dB':float(np.trace(rho@M).real),
            'rho':rho,**op}


def mixed_static_EB(hoppings,positions,P,spec=Spec(),frame='symmetric',K=None):
    """-Tr[P L_f(H;M)]-Tr[f(H) d_B P], thermodynamic STATIC, fixed mu."""
    state=equilibrium_derivatives(hoppings,positions,P,spec,frame,K)
    drho=finite_fermi_frechet(state['H'],state['M'],spec.mu,spec.temperature)
    bulk=-float(np.trace(state['P']@drho).real)
    contact=-float(np.trace(state['rho']@state['dP']).real)
    return {'mixed_d2F_dEdB':bulk+contact,'fermi_density_derivative':bulk,
            'position_contact_dP_dB':contact,'dF_dE':state['dF_dE'],
            'dF_dB':state['dF_dB'],'frame':frame,
            'status':'VERIFIED_STATIC_FINITE_EQUILIBRIUM_MAXWELL_NOT_DRIVEN_ORBITAL_EB'}


def finite_difference_mixed(hoppings,positions,P,spec=Spec(),frame='symmetric',K=None,step=2e-4,route='magnetic_of_electric'):
    """Independent re-diagonalization in E,B; no Frechet operator is used."""
    if not np.isfinite(step) or step<=0 or route not in ('magnetic_of_electric','electric_of_magnetic','free_energy'):
        raise ValueError('step and route')
    b=spec.field;h=step
    def obs(B,E,which):
        return equilibrium_derivatives(hoppings,positions,P,replace(spec,field=B),frame,K,E)[which]
    if route=='magnetic_of_electric':return (obs(b+h,0,'dF_dE')-obs(b-h,0,'dF_dE'))/(2*h)
    if route=='electric_of_magnetic':return (obs(b,h,'dF_dB')-obs(b,-h,'dF_dB'))/(2*h)
    return (obs(b+h,h,'F')-obs(b+h,-h,'F')-obs(b-h,h,'F')+obs(b-h,-h,'F'))/(4*h*h)


def canonical_mixed_static_EB(hoppings,positions,P,spec=Spec(),frame='symmetric',K=None):
    """Fixed N susceptibility and analytic magnetic mu_B Ward correction."""
    op=operators(hoppings,positions,P,spec,frame,K)
    e,u=np.linalg.eigh(op['H']);f=fermi(e,spec.mu,spec.temperature)
    rho=(u*f[None,:])@u.conj().T;fp=-f*(1-f)/spec.temperature
    nmu=float(-fp.sum())
    if nmu<1e-14:raise ValueError('fixed N compressibility too small')
    drho=finite_fermi_frechet(op['H'],op['M'],spec.mu,spec.temperature)
    nB=float(np.trace(drho).real);muB=-nB/nmu
    drho_N=drho+u@np.diag(-fp*muB)@u.conj().T
    residual=float(abs(np.trace(drho_N)))
    if residual>2e-10:raise RuntimeError('particle-number Ward failed')
    bulk=-float(np.trace(op['P']@drho_N).real)
    contact=-float(np.trace(rho@op['dP']).real)
    return {'mixed_d2F_dEdB_fixed_N':bulk+contact,'fixed_N_fermi_term':bulk,
            'position_contact_dP_dB':contact,'reference_N':float(f.sum()),
            'mu_B_fixed_N':float(muB),'dN_dB_at_mu':nB,'dN_dmu':nmu,
            'particle_number_derivative_residual':residual,
            'status':'VERIFIED_STATIC_CANONICAL_MAXWELL_NOT_TRANSPORT'}


def finite_difference_canonical(hoppings,positions,P,spec=Spec(),frame='symmetric',K=None,
                                step=3e-4,route='magnetic_of_electric'):
    """Numerical refit mu(E,B) to fixed N, independent of analytic mu_B."""
    if not np.isfinite(step) or step<=0 or route not in ('magnetic_of_electric','electric_of_magnetic','free_energy'):
        raise ValueError('step and route')
    ini=operators(hoppings,positions,P,spec,frame,K)
    target=float(fermi(np.linalg.eigvalsh(ini['H']),spec.mu,spec.temperature).sum())
    if not 1e-12<target<len(ini['H'])-1e-12:raise ValueError('nontrivial canonical N')
    def obs(B,E,which):
        news=replace(spec,field=B)
        op=operators(hoppings,positions,P,news,frame,K)
        ev=np.linalg.eigvalsh(op['H']-E*op['P']);t=spec.temperature
        mu=brentq(lambda m:float(fermi(ev,m,t).sum()-target),
                  float(min(ev)-50*t-50),float(max(ev)+50*t+50),xtol=1e-13)
        val=equilibrium_derivatives(hoppings,positions,P,replace(news,mu=mu),frame,K,E)
        return val['F']+mu*target if which=='F' else val[which]
    b=spec.field;h=step
    if route=='magnetic_of_electric':return (obs(b+h,0,'dF_dE')-obs(b-h,0,'dF_dE'))/(2*h)
    if route=='electric_of_magnetic':return (obs(b,h,'dF_dB')-obs(b,-h,'dF_dB'))/(2*h)
    return (obs(b+h,h,'F')-obs(b+h,-h,'F')-obs(b-h,h,'F')+obs(b-h,-h,'F'))/(4*h*h)


def benchmark_graph(tr_breaking=True):
    """Finite model only: explicitly supplied Peierls graph and dipoles."""
    xy=np.array([[0.,0.],[1.15,.08],[1.07,1.12],[-.12,1.03]])
    h=np.diag([-.8,-.27,.34,.93]).astype(complex)
    for i,j,z in ((0,1,.72+.09j),(1,2,.55-.13j),(2,3,.65+.12j),(3,0,.81+.14j),(0,2,.19-.07j)):
        if not tr_breaking:z=z.real
        h[i,j]=z;h[j,i]=np.conj(z)
    d=np.zeros((4,4),complex)
    for i,j,z in ((0,1,.13+.05j),(2,3,-.17+.07j),(0,2,.09-.04j)):
        if not tr_breaking:z=z.real
        d[i,j]=z;d[j,i]=np.conj(z)
    return h,xy,d
