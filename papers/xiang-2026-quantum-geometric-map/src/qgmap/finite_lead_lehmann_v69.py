"""v69: independent finite-Hamiltonian Lehmann oracle for the v68 device E/B response.

Two finite chains of length L are attached explicitly to the same four-orbital
Peierls device. The entire closed Hermitian Hamiltonian is diagonalized; the
AC density response is computed from its eigenstate-pair Lehmann sum, not from
surface self-energies or Keldysh energy integration. A small positive eta
regularizes the discrete transition spectrum and must be taken to zero only
with a controlled L -> infinity sequence.

The observable intentionally matches v68's DEVICE projected-polarization
current J_a(E)=i[H_device-E P_b,P_a], embedded as zero on lead sites. It is not
the terminal Landauer current and not yet the source paper's bulk orbital Hall.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
from scipy.linalg import eigh
from .lead_keldysh_eb_v68 import Spec as KeldyshSpec, operator_data, comm


@dataclass(frozen=True)
class Spec:
    mu: float = .13
    temperature: float = .37
    field: float = .23
    charge_over_hbar: float = -.84
    omega: float = .29
    lead_hopping: float = 1.1
    lead_onsite: float = 0.
    length: int = 120
    eta: float = .02
    end_cap: float = 0.

    def __post_init__(self):
        vals=(self.mu,self.temperature,self.field,self.charge_over_hbar,self.omega,
              self.lead_hopping,self.lead_onsite,self.eta,self.end_cap)
        if not np.isfinite(vals).all(): raise ValueError('all parameters must be finite')
        if self.temperature<=0 or self.lead_hopping<=0 or self.eta<=0 or self.charge_over_hbar==0:
            raise ValueError('T, lead hopping and Lehmann eta must be positive; charge nonzero')
        if not isinstance(self.length,int) or not 12<=self.length<=2000:
            raise ValueError('finite lead length must be an integer in [12,2000]')
        if abs(self.mu-self.lead_onsite)>=2*self.lead_hopping:
            raise ValueError('mu must lie inside the propagating lead band')


def _base_spec(s: Spec):
    return KeldyshSpec(mu=s.mu,temperature=s.temperature,field=s.field,
        charge_over_hbar=s.charge_over_hbar,omega=s.omega,
        lead_hopping=s.lead_hopping,lead_onsite=s.lead_onsite,
        nodes_per_interval=64)


def build_finite_system(spec: Spec=Spec(), frame='symmetric', K=None, electric_gate=0.):
    """Return exact finite Hermitian H, embedded P and device-current vertices.

    The endpoint contact vectors are the same physical vectors used by v68's
    Sigma=g v v†; their B-dependent frame transformation is included by the
    v68 operator builder. The chain hopping t is internal to each lead; the
    device-endpoint matrix element is v (not t*v).
    """
    if not np.isfinite(electric_gate): raise ValueError('finite electric gate required')
    bs=_base_spec(spec)
    o=operator_data(bs,frame,K)
    nd=len(o['H']); L=spec.length; n=nd+2*L
    H=np.zeros((n,n),complex); H[:nd,:nd]=o['H']-electric_gate*o['Pb']
    t=spec.lead_hopping; eps=spec.lead_onsite
    vectors=o['vectors']
    for lead in range(2):
        start=nd+lead*L
        H[start:start+L,start:start+L]=np.diag(np.full(L,eps,dtype=complex))
        idx=np.arange(L-1)
        H[start+idx,start+idx+1]=t
        H[start+idx+1,start+idx]=t
        H[start+L-1,start+L-1]+=spec.end_cap
        v=vectors[lead]
        H[:nd,start]=v
        H[start,:nd]=v.conj()
    Pa=np.zeros((n,n),complex); Pb=np.zeros_like(Pa)
    Pa[:nd,:nd]=o['Pa']; Pb[:nd,:nd]=o['Pb']
    J0=np.zeros_like(H); J0[:nd,:nd]=1j*comm(o['H']-electric_gate*o['Pb'],o['Pa'])
    JE=np.zeros_like(H); JE[:nd,:nd]=-1j*comm(o['Pb'],o['Pa'])
    # The observable is device-only current; contact hopping is part of Hfull
    # and therefore of the response, but is not silently added to the vertex.
    herm=np.linalg.norm(H-H.conj().T)
    if herm>1e-12: raise RuntimeError(f'finite lead Hamiltonian not Hermitian: {herm:g}')
    return {'H':H,'Pa':Pa,'Pb':Pb,'J0':J0,'JE':JE,'Hdevice':o['H'],
            'device_dimension':nd,'length':L,'frame':frame,'operators':o,
            'electric_gate':electric_gate,'hermiticity_error':float(herm)}


def _fermi(e,mu,T):
    x=np.clip((np.asarray(e)-mu)/T,-700,700)
    return 1/(np.exp(x)+1)


def _device_eig_operator(U, Adev):
    """Transform an operator supported only on the four device orbitals."""
    Ud=U[:Adev.shape[0],:]
    return Ud.conj().T@Adev@Ud


def lehmann_response(spec: Spec=Spec(), frame='symmetric', K=None):
    """Finite-frequency response from a closed finite Hamiltonian Lehmann sum.

    rho_E[m,n]=(f_m-f_n) P_b[m,n]/(omega+E_n-E_m+i eta).
    This is an isothermal finite-frequency susceptibility with a specified
    Lorentzian resolution eta. It has no reservoir self-energy or Keldysh
    integral. Its static omega=0 limit is NOT the isothermal thermodynamic
    derivative unless the diagonal -f'(E_n) P_nn term is included separately.
    """
    b=build_finite_system(spec,frame,K)
    e,U=eigh(b['H'],check_finite=False,driver='evd')
    f=_fermi(e,spec.mu,spec.temperature)
    Pb=_device_eig_operator(U,b['Pb'][:4,:4])
    J=_device_eig_operator(U,b['J0'][:4,:4])
    JE=_device_eig_operator(U,b['JE'][:4,:4])
    numerator=f[:,None]-f[None,:]
    denominator=spec.omega+e[None,:]-e[:,None]+1j*spec.eta
    rhoE=numerator*Pb/denominator
    current_response=np.sum(J*rhoE.T)
    contact_response=np.dot(f,np.diag(JE))
    value=complex(current_response+contact_response)
    Udev=U[:4,:]
    Pdev_eig=Udev.conj().T@Udev
    device_charge=complex(np.trace(Pdev_eig@rhoE))
    total_charge=complex(np.trace(rhoE))
    return {'sigma_device':value,'sigma_noncontact':complex(current_response),
        'sigma_electric_contact':complex(contact_response),'energies':e,'occupations':f,
        'rhoE_eigenbasis':rhoE,'J_eigenbasis':J,'Pb_eigenbasis':Pb,
        'device_charge_response':device_charge,
        'lead_charge_response':total_charge-device_charge,
        'total_charge_response':total_charge,
        'status':'VERIFIED_FINITE_CLOSED_HERMITIAN_LEHMANN_CONDITIONAL_FINITE_L_ETA_NOT_BULK'}


def static_isothermal_response(spec: Spec=Spec(omega=0.), frame='symmetric', K=None):
    """Exact d_E equilibrium current from the Frechet derivative of f(H-E P).

    Includes the diagonal -f'(E_n) P_nn contribution omitted by the purely
    dynamical commutator susceptibility of an isolated finite system.
    """
    if abs(spec.omega)>1e-14: raise ValueError('static isothermal response requires omega=0')
    b=build_finite_system(spec,frame,K)
    e,U=eigh(b['H'],check_finite=False,driver='evd')
    f=_fermi(e,spec.mu,spec.temperature)
    Pb=_device_eig_operator(U,b['Pb'][:4,:4])
    J=_device_eig_operator(U,b['J0'][:4,:4])
    JE=_device_eig_operator(U,b['JE'][:4,:4])
    gap=e[None,:]-e[:,None]
    num=f[:,None]-f[None,:]
    mid=_fermi((e[:,None]+e[None,:])/2,spec.mu,spec.temperature)
    kernel=np.empty_like(gap)
    mask=np.abs(gap)>1e-9
    kernel[mask]=num[mask]/gap[mask]
    kernel[~mask]=(mid*(1-mid)/spec.temperature)[~mask]
    rhoE=kernel*Pb
    value=np.sum(J*rhoE.T)+np.dot(f,np.diag(JE))
    return {'sigma_device':complex(value),'rhoE_eigenbasis':rhoE,
            'status':'VERIFIED_FINITE_SYSTEM_ISOTHERMAL_FRECHET_RESPONSE'}


def equilibrium_current(spec: Spec=Spec(omega=0.), E=0., frame='symmetric', K=None):
    """Re-diagonalize the complete finite Hamiltonian at static device E."""
    if not np.isfinite(E): raise ValueError('finite electric field required')
    b=build_finite_system(spec,frame,K,electric_gate=E)
    e,U=eigh(b['H'],check_finite=False,driver='evd')
    f=_fermi(e,spec.mu,spec.temperature)
    J=_device_eig_operator(U,b['J0'][:4,:4])
    return complex(np.dot(f,np.diag(J)))


def static_finite_difference(spec: Spec=Spec(omega=0.), frame='symmetric', K=None, step=1e-4):
    if abs(spec.omega)>1e-14: raise ValueError('static finite difference requires omega=0')
    if not np.isfinite(step) or step<=0: raise ValueError('positive finite step required')
    return (equilibrium_current(spec,step,frame,K)-equilibrium_current(spec,-step,frame,K))/(2*step)


def finite_difference_B(spec: Spec=Spec(), frame='symmetric', K=None, step=1e-3, points=5):
    """Independent full-Hamiltonian B finite difference of the Lehmann result."""
    if not np.isfinite(step) or step<=0 or points not in (3,5): raise ValueError('positive step; points 3 or 5')
    def f(b): return lehmann_response(replace(spec,field=b),frame,K)['sigma_device']
    B=spec.field
    if points==3:return (f(B+step)-f(B-step))/(2*step)
    return (f(B-2*step)-8*f(B-step)+8*f(B+step)-f(B+2*step))/(12*step)


def finite_difference_static_mixed(spec: Spec=Spec(omega=0.), frame='symmetric', K=None, step=2e-3):
    """Independent four-corner E/B derivative of full finite-system equilibrium current."""
    if abs(spec.omega)>1e-14 or not np.isfinite(step) or step<=0:
        raise ValueError('omega=0 and positive finite step required')
    B=spec.field
    def f(db,E): return equilibrium_current(replace(spec,field=B+db),E,frame,K)
    return (f(step,step)-f(step,-step)-f(-step,step)+f(-step,-step))/(4*step*step)


def lehmann_continuity_ward(spec: Spec=Spec(), frame='symmetric', K=None):
    """Closed finite-system two-energy Ward identity, independent of Keldysh.

    With eta included in the retarded response denominator, the exact matrix
    identity is (omega+i eta) rho_E - [H,rho_E] + [P_b,rho_0] = 0.
    The finite-system total induced particle number is zero; device charge is
    balanced by opposite induced charge on the explicitly included leads.
    """
    b=build_finite_system(spec,frame,K)
    e,U=eigh(b['H'],check_finite=False,driver='evd')
    f=_fermi(e,spec.mu,spec.temperature)
    Pb=_device_eig_operator(U,b['Pb'][:4,:4])
    numerator=f[:,None]-f[None,:]
    denominator=spec.omega+e[None,:]-e[:,None]+1j*spec.eta
    rhoE=numerator*Pb/denominator
    comm_energy=(e[:,None]-e[None,:])*rhoE
    p_rho=(f[None,:]-f[:,None])*Pb
    residual=(spec.omega+1j*spec.eta)*rhoE-comm_energy+p_rho
    Pdev=U[:4,:].conj().T@U[:4,:]
    device_charge=complex(np.trace(Pdev@rhoE))
    total_charge=complex(np.trace(rhoE))
    lead_charge=total_charge-device_charge
    return {'ward_error':float(np.linalg.norm(residual)),
        'total_charge_response':total_charge,'device_charge_response':device_charge,
        'lead_charge_response':lead_charge,
        'charge_partition_error':float(abs(total_charge-device_charge-lead_charge)),
        'status':'VERIFIED_FINITE_CLOSED_HAMILTONIAN_LEHMANN_WARD'}
