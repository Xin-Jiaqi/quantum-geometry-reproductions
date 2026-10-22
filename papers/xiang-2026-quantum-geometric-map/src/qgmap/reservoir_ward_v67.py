"""v67 independent conserving two-reservoir Peierls-scattering benchmark.

Two *physical semi-infinite, single-channel, noninteracting tight-binding
reservoirs*, not phenomenological scalar-gamma relaxation. The full finite
Peierls device H(B), supplied projected P_x(B),P_y(B), and B-frame connections
come from v65/v66. A static E gate H(B,E)=H(B)-E P_y(B) is DISTINCT from
reservoir voltage bias, which drives the two-terminal Landauer current.

Observable: transmission T(E,B;mu) and mixed d_E d_B T, not original orbital
sigma_xy, nor a same-observable match to v66 finite-frequency electric drive.
Transmission is a conserving coherent elastic noninteracting benchmark; no
microscopic many-body collisions or thermodynamic bulk limit are claimed.
"""
from dataclasses import dataclass
import numpy as np
from .equilibrium_orbital_eb_v65 import operators,Spec as OrbitalSpec
from .dynamic_orbital_eb_v66 import benchmark_model


def _comm(a,b): return a@b-b@a


@dataclass(frozen=True)
class ReservoirSpec:
    mu: float=.13
    field: float=.23
    charge_over_hbar:float=-.84
    electric_gate:float=.04
    lead_hopping:float=1.1
    lead_onsite:float=0.
    lead_regulator:float=1e-10
    def __post_init__(self):
        vals=(self.mu,self.field,self.charge_over_hbar,self.electric_gate,
              self.lead_hopping,self.lead_onsite,self.lead_regulator)
        if (not np.isfinite(vals).all() or self.lead_hopping<=0 or self.lead_regulator<0
                or self.charge_over_hbar==0 or abs(self.mu-self.lead_onsite)>=2*self.lead_hopping):
            raise ValueError('finite parameters, physical hopping>0, regulator>=0, charge!=0')


def chain_surface(spec):
    """Retarded g00 for an exactly semi-infinite single-orbital nearest-neighbor chain.

    g=1/(z-eps0-t²g), choose Im g < 0 for propagating channel. No
    wide-band self-energy approximation or fitted device relaxation rate.
    """
    z=complex(spec.mu-spec.lead_onsite,spec.lead_regulator)
    root=np.sqrt(z*z-4*spec.lead_hopping**2+0j)
    if root.imag<0:root=-root
    g=(z-root)/(2*spec.lead_hopping**2)
    if g.imag>1e-12:raise RuntimeError('nonretarded surface solution')
    return g


def _lead_vectors(n=4):
    if n!=4: raise ValueError('four-site benchmark; physical lead couplings must be supplied anew for other graphs')
    l=np.array([.48+.07j,.39-.06j,0,0],complex)
    r=np.array([0,0,.52-.09j,0],complex)
    return l,r


def setup(hoppings,positions,Pa,Pb,spec=ReservoirSpec(),frame='symmetric',K=None,
          transform_contacts=True):
    """Exact B-dependent contact selfenergy and full magnetic derivative.

    transform_contacts=False is a scientific NEGATIVE control: changes physics
    in B-dependent frames and is not a valid gauge transformation.
    """
    o=OrbitalSpec(mu=spec.mu,temperature=.37,field=spec.field,
                  charge_over_hbar=spec.charge_over_hbar)
    a=operators(hoppings,positions,Pa,o,frame,K)
    b=operators(hoppings,positions,Pb,o,frame,K)
    H0=a['H'];M0=a['M']; A=a['P'];dA=a['dP']; B=b['P'];dB=b['dP']
    H=H0-spec.electric_gate*B
    M=M0-spec.electric_gate*dB
    l0,r0=_lead_vectors(len(H))
    U,dU=a['U'],a['dU']
    if transform_contacts:
        vs=[U@x for x in (l0,r0)]; dvs=[dU@x for x in (l0,r0)]
    else:
        vs=[l0,r0];dvs=[np.zeros_like(l0),np.zeros_like(r0)]
    g=chain_surface(spec)
    S=[g*np.outer(v,v.conj()) for v in vs]
    Sd=[g*(np.outer(dv,v.conj())+np.outer(v,dv.conj())) for v,dv in zip(vs,dvs)]
    G=[1j*(s-s.conj().T) for s in S]
    Gd=[1j*(s-s.conj().T) for s in Sd]
    Gr=np.linalg.inv(spec.mu*np.eye(len(H))-H-sum(S))
    GrB=Gr@(M+sum(Sd))@Gr
    Ag=Gr.conj().T;AgB=GrB.conj().T
    spectral_err=np.linalg.norm(1j*(Gr-Ag)-Gr@sum(G)@Ag)
    return {'H':H,'M':M,'Pa':A,'dPa':dA,'Pb':B,'dPb':dB,'U':U,'dU':dU,
            'S':S,'Sd':Sd,'Gamma':G,'GammaB':Gd,'Gr':Gr,'GrB':GrB,
            'Ga':Ag,'GaB':AgB,'vectors':vs,'vectorB':dvs,'g':g,
            'spectral_identity_error':float(spectral_err),'spec':spec}


def conductance(s):
    GL,GR=s['Gamma']; G=s['Gr']; A=s['Ga']
    return float(np.trace(GL@G@GR@A).real)


def conductance_B(s,omit_contact_B=False):
    """Analytic derivative of exact two-terminal Landauer transmission."""
    GL,GR=s['Gamma'];G=s['Gr'];A=s['Ga']
    if omit_contact_B:
        dGL=np.zeros_like(GL);dGR=np.zeros_like(GR)
        GB=G@s['M']@G
    else:
        dGL,dGR=s['GammaB'];GB=s['GrB']
    dA=GB.conj().T
    terms=(dGL@G@GR@A+GL@GB@GR@A+GL@G@dGR@A+GL@G@GR@dA)
    return float(np.trace(terms).real)


def conductance_E(s):
    G=s['Gr'];A=s['Ga'];GL,GR=s['Gamma'];B=s['Pb']
    GE=-G@B@G
    return float(np.trace(GL@GE@GR@A+GL@G@GR@GE.conj().T).real)


def conductance_EB(s,omit_magnetic_contacts=False,omit_position_contact=False):
    """d_B d_E T with all eight differentiated Caroli terms.

    Contact B derivatives include *both* selfenergy and linewidth. Pb_B is
    essential under orbital magnetic basis rotations. The omitted branches
    exist only as negative controls, never as physical gauge equivalents.
    """
    G=s['Gr']; A=s['Ga']; GL,GR=s['Gamma'];B=s['Pb']
    if omit_magnetic_contacts:
        GB=G@s['M']@G; Lb=np.zeros_like(GL);Rb=np.zeros_like(GR)
    else:
        GB=s['GrB'];Lb,Rb=s['GammaB']
    dB=np.zeros_like(B) if omit_position_contact else s['dPb']
    GE=-G@B@G
    GEB=-(GB@B@G+G@dB@G+G@B@GB)
    t=(Lb@GE@GR@A + GL@GEB@GR@A+GL@GE@Rb@A+GL@GE@GR@GB.conj().T+
       Lb@G@GR@GE.conj().T+GL@GB@GR@GE.conj().T+
       GL@G@Rb@GE.conj().T+GL@G@GR@GEB.conj().T)
    return float(np.trace(t).real)


def local_polarization_bias(s):
    """Gate-independent small *voltage bias* rho, NOT homogeneous E-field response.

    rho_delta=(G ΓL G†-G ΓR G†)/2 for V_L=+1/2,V_R=-1/2.
    P_a-current Tr(i[H,P_a] rho_delta) is an internal diagnostic, NOT
    automatically the physical terminal current or original orbital Hall.
    """
    G=s['Gr'];A=s['Ga'];GL,GR=s['Gamma'];GB=s['GrB'];AB=s['GaB']
    net=GL-GR;dnet=s['GammaB'][0]-s['GammaB'][1]
    rho=.5*G@net@A
    rhoB=.5*(GB@net@A+G@dnet@A+G@net@AB)
    H,Pa,M,dPa=s['H'],s['Pa'],s['M'],s['dPa']
    J=1j*_comm(H,Pa);JB=1j*(_comm(M,Pa)+_comm(H,dPa))
    value=float(np.trace(J@rho).real)
    derivative=float(np.trace(JB@rho+J@rhoB).real)
    # Local physical site3 number projector, transformed with device B frame.
    # Physical lead coupling to site3 is zero, and [H,rho]_site3=0.
    U=s['U'];dU=s['dU'];Q0=np.diag([0.,0.,0.,1.]);Q=U@Q0@U.conj().T
    QB=dU@Q0@U.conj().T+U@Q0@dU.conj().T
    continuity=np.trace(Q@_comm(H,rho))
    continuity_B=np.trace(QB@_comm(H,rho)+Q@(_comm(M,rho)+_comm(H,rhoB)))
    return {'internal_pa_bias_current':value,'internal_pa_bias_current_B':derivative,
            'uncontacted_site3_continuity':float(abs(continuity)),
            'uncontacted_site3_continuity_B':float(abs(continuity_B)),
            'rho_bias':rho,'rho_bias_B':rhoB,'J':J,'JB':JB}


def scattering(s):
    """Exact two-channel S, independent unitarity oracle for lead current Ward."""
    g=s['g'];gamma=-2*g.imag
    if gamma<=0:raise ValueError('no propagating lead channel')
    W=np.stack(s['vectors'],axis=1)*np.sqrt(gamma)
    S=np.eye(2,dtype=complex)-1j*W.conj().T@s['Gr']@W
    return {'S':S,'unitarity_error':float(np.linalg.norm(S.conj().T@S-np.eye(2))),
            'cross_probability':float(abs(S[1,0])**2),
            'left_probability_sum':float(np.sum(abs(S[:,0])**2)),
            'right_probability_sum':float(np.sum(abs(S[:,1])**2))}


def reference(frame='symmetric',K=None,spec=ReservoirSpec(),transform_contacts=True):
    h,xy,pa,pb,k=benchmark_model()
    if frame=='orbital' and K is None:K=k
    return setup(h,xy,pa,pb,spec,frame,K,transform_contacts)


def finite_difference(spec=ReservoirSpec(),frame='symmetric',K=None,
                      parameter='B',step=1e-4,points=5,transform_contacts=True):
    """Independent finite-field fully recomputed scattering (not analytic derivatives)."""
    from dataclasses import replace
    if not np.isfinite(step) or step<=0 or points not in (3,5) or parameter not in ('B','E'):
        raise ValueError('3/5, B/E and positive finite step')
    def val(t):
        p=replace(spec,field=t) if parameter=='B' else replace(spec,electric_gate=t)
        return conductance(reference(frame,K,p,transform_contacts))
    x=spec.field if parameter=='B' else spec.electric_gate
    if points==3:return (val(x+step)-val(x-step))/(2*step)
    return (val(x-2*step)-8*val(x-step)+8*val(x+step)-val(x+2*step))/(12*step)


def finite_difference_mixed(spec=ReservoirSpec(),frame='symmetric',K=None,
                            step=2e-3,points=5,transform_contacts=True):
    """Independent mixed conductance FD; 4 corner (3pt), 16 corner (5pt)."""
    from dataclasses import replace
    if not np.isfinite(step) or step<=0 or points not in (3,5):raise ValueError('positive step and 3/5 finite difference')
    b=spec.field;e=spec.electric_gate
    def f(db,de):return conductance(reference(frame,K,replace(spec,field=b+db,electric_gate=e+de),transform_contacts))
    if points==3:return (f(step,step)-f(step,-step)-f(-step,step)+f(-step,-step))/(4*step*step)
    offsets=(-2,-1,1,2);weights=(1,-8,8,-1)
    return sum(wb*we*f(db*step,de*step) for db,wb in zip(offsets,weights)
               for de,we in zip(offsets,weights))/(144*step*step)
