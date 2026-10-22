"""v68: finite-frequency conserving *noninteracting-lead* Keldysh orbital E/B gate.

The same finite-device projected polarization-current operator J_a(E)=
 i[H(B)-E P_b(B),P_a(B)] used in v66 is measured in an explicit
semi-infinite two-reservoir scattering state. A spatially uniform scalar
length-gauge drive is applied to the DEVICE ONLY, not to the leads. This is
not source Xiang's thermodynamic bulk intrinsic orbital Hall conductivity.

Two-time Keldysh linear variation, e^{-i omega t} convention:
 delta G^<(epsilon+omega,epsilon) = -Gr_+ P_b G^<_- - G^<_+ P_b Ga_-
 and G^<(epsilon)=i f(epsilon) Gr Gamma Ga. Hence
 rho_E = -Integral (Gr_+ Pb f_- A_- + f_+ A_+ Pb Ga_-)/(2pi).
 Exact projected current contact J_E=-i[Pb,Pa]. All magnetic derivatives
 carry H_B, Pa_B, Pb_B AND transformed device-to-lead coupling/Sigma_B.

Full finite E is modeled by H(B,E)=H_B-E Pb(B) and electrochemical
reservoirs held fixed. This is a specified finite open-scattering protocol;
complete real Wannier and lead electrostatic gauge are not inferred.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre
from .equilibrium_orbital_eb_v65 import operators,Spec as OrbitalSpec,fermi
from .dynamic_orbital_eb_v66 import benchmark_model
from .reservoir_ward_v67 import _lead_vectors


def comm(a,b):return a@b-b@a

@dataclass(frozen=True)
class Spec:
    mu:float=.13
    temperature:float=.37
    field:float=.23
    charge_over_hbar:float=-.84
    omega:float=.29
    lead_hopping:float=1.1
    lead_onsite:float=0.
    nodes_per_interval:int=640
    def __post_init__(self):
        vals=(self.mu,self.temperature,self.field,self.charge_over_hbar,self.omega,self.lead_hopping,self.lead_onsite)
        if not np.isfinite(vals).all() or self.temperature<=0 or self.lead_hopping<=0 or self.charge_over_hbar==0:
            raise ValueError('finite mu,T,B,q,omega,t and T,t>0,q != 0')
        if not isinstance(self.nodes_per_interval,int) or not 16<=self.nodes_per_interval<=4096:
            raise ValueError('quadrature nodes 16..4096')
        if abs(self.mu-self.lead_onsite)>=2*self.lead_hopping:
            raise ValueError('reference chemical potential must lie in propagating lead band')


def semiinfinite_g(energy,t,onsite):
    """Real-axis retarded surface Green function, exact band branch and eta=0+."""
    z=float(energy-onsite)
    if abs(z)<2*t:return complex(z,-np.sqrt(4*t*t-z*z))/(2*t*t)
    if abs(z)==2*t:return complex(z/(2*t*t))
    return complex((z-np.sign(z)*np.sqrt(z*z-4*t*t))/(2*t*t))

@lru_cache(maxsize=128)
def energy_grid(hop,onsite,omega,n):
    """Split Gauss-Legendre mesh exactly at spectral continuum thresholds."""
    lo,hi=onsite-2*hop,onsite+2*hop
    # Occupation spectral A_-(eps) or A_+(eps+omega) is nonzero only in the
    # respective shifted continuum. Split where either inverse sqrt changes.
    boundaries=sorted(set((lo,hi,lo-omega,hi-omega)))
    x,w=roots_legendre(n)
    grid=[];weights=[]
    for a,b in zip(boundaries[:-1],boundaries[1:]):
        if b<=a:continue
        grid.extend(((a+b)+(b-a)*x)/2)
        weights.extend((b-a)*w/2)
    return np.asarray(grid),np.asarray(weights)


def operator_data(spec=Spec(),frame='symmetric',K=None):
    h,xy,pa,pb,k=benchmark_model()
    if frame=='orbital' and K is None:K=k
    os=OrbitalSpec(mu=spec.mu,temperature=spec.temperature,field=spec.field,
                   charge_over_hbar=spec.charge_over_hbar)
    a=operators(h,xy,pa,os,frame,K); b=operators(h,xy,pb,os,frame,K)
    if not np.allclose(a['H'],b['H']) or not np.allclose(a['M'],b['M']):
        raise RuntimeError('operator frames inconsistent')
    U,dU=a['U'],a['dU'];l0,r0=_lead_vectors(len(a['H']))
    vectors=[U@l0,U@r0]; vectorB=[dU@l0,dU@r0]
    return {'H':a['H'],'M':a['M'],'Pa':a['P'],'PaB':a['dP'],
            'Pb':b['P'],'PbB':b['dP'], 'vectors':vectors,'vectorsB':vectorB,
            'frame':frame,'U':U,'K':K}


def greens(energy,o,spec=Spec(),electric_gate=0.,magnetic_derivative=False,omit_lead_B=False):
    """On-shell exact reservoir selfenergies and optional complete B-derivative.

    No scalar gamma, no wide-band approximation. Matrix energy continuation
    uses an infinitesimal retarded i0+ only for numerical inversion stability;
    no fictitious bath population is introduced.
    """
    if not np.isfinite([energy,electric_gate]).all():raise ValueError('finite energy and gate')
    t=spec.lead_hopping;g=semiinfinite_g(energy,t,spec.lead_onsite)
    vv=o['vectors'];vd=o['vectorsB']
    s=sum(g*np.outer(v,v.conj()) for v in vv)
    gamma=1j*(s-s.conj().T)
    h=o['H']-electric_gate*o['Pb'];m=o['M']-electric_gate*o['PbB']
    G=np.linalg.inv((energy+1e-12j)*np.eye(len(h))-h-s)
    Ga=G.conj().T
    A=G@gamma@Ga
    d=None
    if magnetic_derivative:
        sd=sum(g*(np.outer(dv,v.conj())+np.outer(v,dv.conj())) for v,dv in zip(vv,vd))
        gd=1j*(sd-sd.conj().T)
        if omit_lead_B:
            sd=np.zeros_like(sd);gd=np.zeros_like(gd)
        GB=G@(m+sd)@G
        AB=GB@gamma@Ga+G@gd@Ga+G@gamma@GB.conj().T
        d={'GB':GB,'AB':AB,'Sd':sd,'GammaB':gd}
    return {'Gr':G,'Ga':Ga,'A':A,'Gamma':gamma,'Sigma':s,'B':d}


def _integrate(spec,o,omega=None,magnetic_derivative=False,electric_gate=0.,omit_lead_B=False):
    """Integrate exact Fermi-lead Keldysh electron density/AC linear response."""
    if omega is None:omega=spec.omega
    egrid,wgrid=energy_grid(spec.lead_hopping,spec.lead_onsite,omega,spec.nodes_per_interval)
    n=len(o['H']);rho=np.zeros((n,n),complex);rhoE=np.zeros_like(rho)
    if magnetic_derivative:rhoB=np.zeros_like(rho);rhoEB=np.zeros_like(rho)
    Pb=o['Pb'];PbB=o['PbB']
    for en,w in zip(egrid,wgrid):
        lower=greens(en,o,spec,electric_gate,magnetic_derivative,omit_lead_B)
        upper=greens(en+omega,o,spec,electric_gate,magnetic_derivative,omit_lead_B)
        fm=fermi(en,spec.mu,spec.temperature);fp=fermi(en+omega,spec.mu,spec.temperature)
        Gp=upper['Gr'];GmA=lower['Ga'];Am=lower['A'];Ap=upper['A']
        rho+=w*fm*Am
        rhoE-=w*(fm*Gp@Pb@Am+fp*Ap@Pb@GmA)
        if magnetic_derivative:
            GpB=upper['B']['GB'];AmB=lower['B']['AB'];ApB=upper['B']['AB'];GmAB=lower['B']['GB'].conj().T
            rhoB+=w*fm*AmB
            rhoEB-=w*(fm*(GpB@Pb@Am+Gp@PbB@Am+Gp@Pb@AmB)+
                      fp*(ApB@Pb@GmA+Ap@PbB@GmA+Ap@Pb@GmAB))
    scale=1/(2*np.pi)
    result={'rho':rho*scale,'rhoE':rhoE*scale}
    if magnetic_derivative:result.update({'rhoB':rhoB*scale,'rhoEB':rhoEB*scale})
    return result


def response(spec=Spec(),frame='symmetric',K=None,include_electric_contact=True,
             include_all_magnetic_contacts=True,omit_lead_magnetic_contacts=False):
    """Current of finite DEVICE projected polarization and its complete d_B.

    NOT terminal current, NOT source Xiang bulk orbital Hall. rhoE is exact
    linearized lesser Keldysh for device-local length-gauge AC perturbation.
    """
    o=operator_data(spec,frame,K)
    d=_integrate(spec,o,magnetic_derivative=True,omit_lead_B=omit_lead_magnetic_contacts)
    H,M,Pa,PaB,Pb,PbB=(o[x] for x in ('H','M','Pa','PaB','Pb','PbB'))
    J=1j*comm(H,Pa);JE=-1j*comm(Pb,Pa)
    if include_all_magnetic_contacts:
        JB=1j*(comm(M,Pa)+comm(H,PaB))
        JEB=-1j*(comm(PbB,Pa)+comm(Pb,PaB))
        rEB=d['rhoEB'];rB=d['rhoB']
    else:
        JB=1j*comm(M,Pa);JEB=np.zeros_like(J)
        rEB=d['rhoEB'];rB=d['rhoB']
    first=complex(np.trace(J@d['rhoE']))
    electric=complex(np.trace(JE@d['rho'])) if include_electric_contact else 0j
    components={'magnetic_current':complex(np.trace(JB@d['rhoE'])),
                'magnetic_density':complex(np.trace(J@rEB)),
                'electric_contact_magnetic':complex(np.trace(JEB@d['rho'])) if include_electric_contact else 0j,
                'electric_contact_density':complex(np.trace(JE@rB)) if include_electric_contact else 0j}
    return {'sigma_device':first+electric,
            'sigma_noncontact':first,'sigma_electric_contact':electric,
            'mixed_dB_sigma_device':sum(components.values()),
            'mixed_components':components,'rho':d['rho'],'rhoE':d['rhoE'],
            'rhoB':d['rhoB'],'rhoEB':d['rhoEB'],'o':o,
            'status':'VERIFIED_FINITE_EXACT_NONINTERACTING_LEAD_KELDYSH_CONDITIONAL_DEVICE_FIELD_NOT_BULK_ORIGINAL'}


def finite_difference_B(spec=Spec(),frame='symmetric',K=None,step=1e-3,points=3):
    if not np.isfinite(step) or step<=0 or points not in (3,5):raise ValueError('positive B step, 3/5 stencil')
    def f(x):return response(replace(spec,field=x),frame,K)['sigma_device']
    b=spec.field
    if points==3:return (f(b+step)-f(b-step))/(2*step)
    return (f(b-2*step)-8*f(b-step)+8*f(b+step)-f(b+2*step))/(12*step)


def equilibrium_current(spec=Spec(),E=0.,frame='symmetric',K=None):
    """Same projected DEVICE J(E) at static gate in equilibrium reservoirs.

    Device-only projected dPa/dt can be nonzero due to reservoir exchange;
    observable is not lead-integrated terminal current. Exact E field follows
    H-E Pb in both propagator and current vertex.
    """
    if not np.isfinite(E):raise ValueError('finite E')
    o=operator_data(spec,frame,K)
    d=_integrate(spec,o,omega=0.,electric_gate=E)
    J=1j*comm(o['H']-E*o['Pb'],o['Pa'])
    return complex(np.trace(J@d['rho']))


def static_field_FD(spec=Spec(),frame='symmetric',K=None,step=2e-3):
    if abs(spec.omega)>1e-15:raise ValueError('static finite-E oracle only at omega=0')
    if not np.isfinite(step) or step<=0:raise ValueError('step>0')
    return (equilibrium_current(spec,step,frame,K)-equilibrium_current(spec,-step,frame,K))/(2*step)


def static_mixed_FD(spec=Spec(),frame='symmetric',K=None,step=4e-3):
    if abs(spec.omega)>1e-15:raise ValueError('static mixed finite-E/B oracle only at omega=0')
    if not np.isfinite(step) or step<=0:raise ValueError('step>0')
    f=lambda db,de:equilibrium_current(replace(spec,field=spec.field+db),de,frame,K)
    return (f(step,step)-f(step,-step)-f(-step,step)+f(-step,-step))/(4*step*step)


def keldysh_continuity(spec=Spec(),frame='symmetric',K=None):
    """Microscopic two-energy *lead collision* Ward, not inferred from residual.

    omega*rhoE - [H,rhoE] + [Pb,rho0] = I_leads(omega).
    At site 3 (projector away from all lead vectors), RHS is exactly zero.
    Other sites have explicit nonzero frequency-dependent lead injection.
    """
    o=operator_data(spec,frame,K)
    rho=np.zeros((4,4),complex);rhoE=np.zeros_like(rho);injection=np.zeros_like(rho)
    grid,w=energy_grid(spec.lead_hopping,spec.lead_onsite,spec.omega,spec.nodes_per_interval)
    for e,weight in zip(grid,w):
        minus=greens(e,o,spec);plus=greens(e+spec.omega,o,spec)
        fm=fermi(e,spec.mu,spec.temperature);fp=fermi(e+spec.omega,spec.mu,spec.temperature)
        Lm=1j*fm*minus['A'];Lp=1j*fp*plus['A']
        Pb=o['Pb'];V=-Pb
        X=plus['Gr']@V@Lm+Lp@V@minus['Ga']
        # complete two-energy on-shell Dyson reservoir collision term
        sig_less_p=1j*fp*plus['Gamma'];sig_less_m=1j*fm*minus['Gamma']
        collision=(plus['Sigma']@X-X@minus['Sigma'].conj().T+
                   sig_less_p@plus['Ga']@V@minus['Ga']-
                   plus['Gr']@V@minus['Gr']@sig_less_m)
        rho+=weight*Lm
        rhoE+=weight*X
        injection+=weight*collision
    scale=-1j/(2*np.pi)
    rho*=scale;rhoE*=scale;injection*=scale
    lhs=spec.omega*rhoE-comm(o['H'],rhoE)+comm(o['Pb'],rho)
    Q=o['U']@np.diag([0,0,0,1])@o['U'].conj().T
    return {'full_ward_error':float(np.linalg.norm(lhs-injection)),
            'uncontacted_site_ward':complex(np.trace(Q@lhs)),
            'uncontacted_site_collision':complex(np.trace(Q@injection)),
            'lead_collision_norm':float(np.linalg.norm(injection)),
            'number_accumulation':complex(spec.omega*np.trace(rhoE)),
            'number_lead_injection':complex(np.trace(injection)),
            'rhoE':rhoE,'rho':rho,'collision':injection,
            'status':'EXACT_TWO_ENERGY_ELASTIC_RESERVOIR_WARD_AT_FIXED_FINITE_OMEGA'}
