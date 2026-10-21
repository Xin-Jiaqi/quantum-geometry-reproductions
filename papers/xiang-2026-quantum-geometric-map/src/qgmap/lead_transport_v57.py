"""v57 finite coherent multiterminal Landauer-Buettiker transport.

G^r=[E-H(B)+i Gamma_total/2]^-1, wide-band local site leads.
T_ab=Tr[Gamma_a G^r Gamma_b G^a] (a != b).
B derivative G'_B=G^r (d_B H) G^r (field-INDEPENDENT Gamma).
4-lead voltage probes satisfy I_top=I_bottom=0, bias V_left=+1/2,
V_right=-1/2. Units: energy and signed qB/hbar are model units;
physical conductance multiplies dimensionless lead currents by q^2/h.
This is a coherent finite-device transport model, NOT a bulk MNHE tensor.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
from numpy.polynomial.legendre import leggauss
from .embedding_magnetic_v55 import peierls_field

@dataclass(frozen=True)
class LeadSpec:
    field:float=0.31
    charge_over_hbar:float=-0.87
    mu:float=0.15
    temperature:float=0.16
    broadening:float=0.55
    energy_points:int=160
    window_in_kT:float=20.0
    def __post_init__(self):
        vals=[self.field,self.charge_over_hbar,self.mu,self.temperature,self.broadening,self.window_in_kT]
        if not np.isfinite(vals).all() or self.charge_over_hbar==0 or self.temperature<=0 or self.broadening<=0 or self.window_in_kT<10 or not isinstance(self.energy_points,int) or self.energy_points<24:
            raise ValueError('finite physical parameters; T, lead width positive; grid >=24; window >=10 kT')


def diamond_graph(n=1,t=1.0):
    """n concentric four-site square diamond rings, increasing radius by layer.

    Each layer has 4 vertices L,T,R,B; adjacent rings coupled radially.
    Four local leads attach only to outermost layer, not a bulk sample.
    """
    if not isinstance(n,int) or not 1<=n<=10 or not np.isfinite(t) or t<=0:
        raise ValueError('1<=rings<=10, hopping positive')
    pos=[]
    for ring in range(1,n+1):
        a=float(ring)
        pos.extend([(-a,0),(0,a),(a,0),(0,-a)])
    pos=np.asarray(pos,float)
    h=np.zeros((4*n,4*n),complex)
    for j in range(n):
        k=4*j
        for i in range(4):
            a,b=k+i,k+(i+1)%4
            h[a,b]=t;h[b,a]=t
        if j:
            for i in range(4):
                a,b=4*(j-1)+i,k+i
                h[a,b]=.63*t;h[b,a]=.63*t
    lead_sites=tuple(4*(n-1)+i for i in range(4))
    return pos,h,lead_sites


def _validate_graph(h0,positions,lead_sites,gammas):
    h0=np.asarray(h0,complex);positions=np.asarray(positions,float)
    n=h0.shape[0] if h0.ndim==2 else 0
    if h0.shape!=(n,n) or n<4 or not np.isfinite(h0).all() or not np.allclose(h0,h0.conj().T,atol=1e-12):
        raise ValueError('finite Hermitian hopping matrix of >=4 sites')
    if positions.shape!=(n,2) or not np.isfinite(positions).all():raise ValueError('site positions shape (N,2)')
    sites=tuple(lead_sites)
    if len(sites) not in (2,4) or len(set(sites))!=len(sites) or any(not isinstance(s,(int,np.integer)) or s<0 or s>=n for s in sites):
        raise ValueError('2 or 4 distinct valid lead site indices')
    gamma=np.array(gammas,float)
    if gamma.shape!=(len(sites),) or not np.isfinite(gamma).all() or np.any(gamma<=0):
        raise ValueError('independent positive wide-band lead broadenings required')
    return h0,positions,sites,gamma


def _weights(spec):
    nodes,gw=leggauss(spec.energy_points)
    x=spec.window_in_kT*nodes
    E=spec.mu+spec.temperature*x
    thermal=1/(4*np.cosh(x/2)**2)
    weights=gw*spec.window_in_kT*thermal
    return E,weights


def _transmission_at_energy(E,h,m,sites,gammas):
    """Analytic d_B T_ab for B-independent local lead self-energies."""
    n=len(h);nl=len(sites)
    width=np.zeros(n,float)
    for site,g in zip(sites,gammas):width[site]+=g
    denom=(E*np.eye(n,dtype=complex)-h+0.5j*np.diag(width))
    gr=np.linalg.inv(denom)
    dgr=gr@m@gr
    ts=np.zeros((nl,nl),float); dt=np.zeros_like(ts)
    for a in range(nl):
        for b in range(nl):
            if a==b:continue
            ga=gammas[a]*gammas[b]
            gab=gr[sites[a],sites[b]]
            dgb=dgr[sites[a],sites[b]]
            ts[a,b]=float(ga*abs(gab)**2)
            dt[a,b]=float(2*ga*np.real(gab.conjugate()*dgb))
    return ts,dt


def transmissions(h0,positions,lead_sites,gammas,spec=LeadSpec(),*,gauge='symmetric'):
    h0,positions,sites,gammas=_validate_graph(h0,positions,lead_sites,gammas)
    pb=peierls_field(h0,positions,spec.field,spec.charge_over_hbar,gauge=gauge)
    energies,weights=_weights(spec)
    t=np.zeros((len(sites),len(sites)),float);dt=t.copy()
    for E,w in zip(energies,weights):
        x,dx=_transmission_at_energy(E,pb.hamiltonian,pb.magnetic_vertex,sites,gammas)
        t+=w*x;dt+=w*dx
    return {'T':t,'dT_dB':dt,'source':'wide_band_local_leads_finite_open_graph',
            'status':'VERIFIED_COHERENT_LEAD_MODEL_NOT_BULK',
            'fermi_window_coverage':float(np.tanh(spec.window_in_kT/2)),
            'limits':'electrically biased four/two lead mesoscopic coherent conductor; energy-independent diagonal Sigma=-i Gamma/2'}


def laplacian(T):
    T=np.asarray(T,float)
    if T.ndim!=2 or T.shape[0]!=T.shape[1] or not np.isfinite(T).all():
        raise ValueError('finite square transmission matrix')
    return -T+np.diag(np.sum(T,axis=1))


def four_terminal(T,dT):
    """Lead order L,T,R,B; V_L=+.5, V_R=-.5, T/B floating."""
    if np.asarray(T).shape!=(4,4) or np.asarray(dT).shape!=(4,4):raise ValueError('4 terminals')
    L=laplacian(T);dL=laplacian(dT)
    fixed=np.array([0,2]);floating=np.array([1,3])
    V=np.zeros(4);V[fixed]=[.5,-.5]
    pp=L[np.ix_(floating,floating)]
    V[floating]=-np.linalg.solve(pp,L[np.ix_(floating,fixed)]@V[fixed])
    dV=np.zeros(4)
    dV[floating]=-np.linalg.solve(pp,dL[np.ix_(floating,np.arange(4))]@V)
    I=L@V
    dI=dL@V+L@dV
    return {'voltages':V,'dV_dB':dV,'lead_currents_h_over_q2':I,
            'dI_dB':dI,'hall_ratio':float(V[1]-V[3]),
            'dHall_dB':float(dV[1]-dV[3]),
            'longitudinal_conductance_q2_over_h':float(I[0]),
            'dLongitudinal_dB':float(dI[0]),
            'max_current_conservation_residual':float(abs(np.sum(I))),
            'max_floating_probe_current':float(max(abs(I[1]),abs(I[3]))),
            'max_derivative_current_conservation_residual':float(abs(np.sum(dI)))}


def device_response(h0,positions,lead_sites,gammas,spec=LeadSpec(),gauge='symmetric'):
    d=transmissions(h0,positions,lead_sites,gammas,spec,gauge=gauge)
    if len(lead_sites)==4:
        d['observables']=four_terminal(d['T'],d['dT_dB'])
    elif len(lead_sites)==2:
        d['observables']={'two_terminal_transmission_q2_over_h':float(d['T'][0,1]),
              'dTransmission_dB':float(d['dT_dB'][0,1])}
    return d


def full_finite_field_fd(h0,positions,lead_sites,gammas,spec=LeadSpec(),gauge='symmetric',step=1e-4):
    if not np.isfinite(step) or step<=0:raise ValueError('positive finite-difference step')
    def t(delta):
        return device_response(h0,positions,lead_sites,gammas,replace(spec,field=spec.field+delta),gauge)
    p2,p1,n1,n2=t(2*step),t(step),t(-step),t(-2*step)
    data={'dT_dB':(-p2['T']+8*p1['T']-8*n1['T']+n2['T'])/(12*step)}
    if len(lead_sites)==4:
        for key in ('hall_ratio','longitudinal_conductance_q2_over_h'):
            data['d_'+key]=(-p2['observables'][key]+8*p1['observables'][key]-8*n1['observables'][key]+n2['observables'][key])/(12*step)
    else:
        key='two_terminal_transmission_q2_over_h'
        data['d_'+key]=float((-p2['observables'][key]+8*p1['observables'][key]-8*n1['observables'][key]+n2['observables'][key])/(12*step))
    return data
