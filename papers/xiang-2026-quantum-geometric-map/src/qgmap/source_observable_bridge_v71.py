"""v71: exact finite-lead polarization contact and literal source spin-sign gates.

VERIFIED only for specified finite Hermitian models and printed-index algebra.
Do not conflate device polarization current with bulk magnetotransport tensor.
"""
from dataclasses import dataclass
from dataclasses import replace
from typing import Mapping
import numpy as np
from scipy.linalg import eigh
from .finite_lead_lehmann_v69 import Spec,build_finite_system,_fermi


def _comm(a,b):
    return a@b-b@a


def _eigen(U,op):
    return U.conj().T@op@U


def finite_contact_identity(spec=Spec(length=24),frame='symmetric',K=None):
    """An exact bridge between device-only and FULL finite-lead Heisenberg current.

    Polarization is supported only on device (zero lead block): this is not
    a physical all-space coordinate or a terminal scattering current.
    The isothermal static derivative uses the full f(H-E Pb), including
    diagonal Fermi-population derivative, unlike the isolated dynamic limit.
    """
    b=build_finite_system(spec,frame,K)
    H,Pa,Pb=b['H'],b['Pa'],b['Pb']
    Jdev=b['J0']
    Jfull=1j*_comm(H,Pa)
    Jcontact=Jfull-Jdev
    JE=b['JE']
    e,U=eigh(H,check_finite=False,driver='evd')
    f=_fermi(e,spec.mu,spec.temperature)
    pb=_eigen(U,Pb)
    ja=_eigen(U,Pa)
    jd=_eigen(U,Jdev)
    jc=_eigen(U,Jcontact)
    je=_eigen(U,JE)
    den=spec.omega+e[None,:]-e[:,None]+1j*spec.eta
    rhoE=(f[:,None]-f[None,:])*pb/den
    response=lambda op:complex(np.sum(op*rhoE.T))
    econtact=complex(np.dot(f,np.diag(je)))
    dev=response(jd)+econtact
    cont=response(jc)
    full=dev+cont
    dipole_response=complex(np.sum(ja*rhoE.T))
    dynamic_expected=-1j*(spec.omega+1j*spec.eta)*dipole_response

    d=e[None,:]-e[:,None]
    diff=f[:,None]-f[None,:]
    ff=_fermi((e[:,None]+e[None,:])/2,spec.mu,spec.temperature)
    slope=ff*(1-ff)/spec.temperature
    tol=1e-9*max(1.,float(np.max(np.abs(e))))
    kernel=np.divide(diff,d,out=np.zeros_like(d),where=np.abs(d)>tol)
    kernel=np.where(np.abs(d)>tol,kernel,slope)
    rhoE_static=kernel*pb
    isothermal=lambda op:complex(np.sum(op*rhoE_static.T))
    static_dev=isothermal(jd)+econtact
    static_cont=isothermal(jc)
    static_full=static_dev+static_cont
    return dict(dynamic_device=dev,dynamic_contact=cont,dynamic_full=full,
      dynamic_full_ward_prediction=dynamic_expected,
      dynamic_ward_residual=float(abs(full-dynamic_expected)),
      static_isothermal_device=static_dev,static_isothermal_contact=static_cont,
      static_isothermal_full=static_full,
      contact_operator_norm=float(np.linalg.norm(Jcontact)),
      current_decomposition_residual=float(np.linalg.norm(Jfull-Jdev-Jcontact)),
      electric_current_contact=econtact,
      frame=frame,
      status='VERIFIED_FINITE_BOUNDED_POLARIZATION_CONTACT_NOT_BULK_TRANSPORT')


def equilibrium_full_current(spec=Spec(length=20,omega=0.),E=0.,frame='symmetric',K=None):
    """Independent full E rediagonalization, never uses Frechet response."""
    if not np.isfinite(E):raise ValueError('E must be finite')
    b=build_finite_system(spec,frame,K,electric_gate=E)
    H,Pa=b['H'],b['Pa']
    eig,U=eigh(H,check_finite=False,driver='evd')
    occ=_fermi(eig,spec.mu,spec.temperature)
    Jfull=1j*_comm(H,Pa)
    Jdev=b['J0']
    Jcont=Jfull-Jdev
    def trace(J):return complex(np.dot(occ,np.diag(_eigen(U,J))))
    return dict(full=trace(Jfull),device=trace(Jdev),contact=trace(Jcont))


def source_eq20_eq31_sign(energies,position,zeeman_spin,df,tau=1.):
    """Literal printed SM Eqs 20 and 31 with documented index convention.

    epsilon_nm=e_n-e_m; v_nm=i*epsilon_nm*r_nm.
    SM20: rho_mn=-tau*S_mn*(df_n-df_m)/epsilon_mn for m != n.
    This contracts to +tau sum_nm Z_nm*df_n where
    Z_nm=i*(r_nm S_mn-r_mn S_nm).
    SM31 *prints* -tau sum Z_nm*df_n. This function NEVER chooses a
    'corrected' author branch or adjusts physical Bbar=-g mu_B B.
    """
    e=np.asarray(energies,float);r=np.asarray(position,complex);s=np.asarray(zeeman_spin,complex)
    f=np.asarray(df,float)
    n=len(e)
    if e.ndim!=1 or n<2 or any(x.shape!=(n,n) for x in (r,s)) or f.shape!=(n,):
        raise ValueError('matching finite band dimensions')
    if not np.isfinite(e).all() or not np.isfinite(f).all() or not np.isfinite(r).all() or not np.isfinite(s).all() or not np.isfinite(tau):
        raise ValueError('finite inputs required')
    if np.min(np.abs(e[:,None]-e[None,:])+np.eye(n)*1e9)<1e-8:
        raise ValueError('this literal interband formula requires nondegenerate bands')
    if not (np.allclose(r,r.conj().T,rtol=0,atol=1e-11) and np.allclose(s,s.conj().T,rtol=0,atol=1e-11)):
        raise ValueError('Hermitian position and spin required')
    eps=e[:,None]-e[None,:]
    v=1j*eps*r
    z=1j*(r*s.T-s*r.T)
    np.fill_diagonal(z,0.)
    rho=np.zeros((n,n),complex)
    for nidx in range(n):
        for midx in range(n):
            if nidx==midx:continue
            # rho_mn with f_n-f_m and epsilon_mn
            rho[midx,nidx]=-tau*s[midx,nidx]*(f[nidx]-f[midx])/(e[midx]-e[nidx])
    direct=complex(np.sum(v*rho.T))
    inferred=complex(tau*np.sum(z*f[:,None]))
    printed=complex(-tau*np.sum(z*f[:,None]))
    return dict(direct_eq20=direct,inferred_Z=inferred,printed_eq31=printed,
                algebra_residual=float(abs(direct-inferred)),
                sign_discrepancy=complex(direct-printed),
                status='LITERAL_SOURCE_INDEX_DISCREPANCY_NOT_AUTHOR_ERRATUM')


def pauli_helicity_gate(kx=.37,ky=-.26,tilt=.13,velocity=1.2):
    """Original documented main Eq10 vs SM Eq62 helicity as a basis test."""
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1.,-1.]).astype(complex)
    h_main=tilt*kx*np.eye(2)+velocity*(kx*sy-ky*sx)
    h_sm=tilt*kx*np.eye(2)+velocity*(ky*sx-kx*sy)
    return dict(hamiltonian_equivalence_error=float(np.linalg.norm(sz@h_main@sz-h_sm)),
                transformed_spin_error=float(np.linalg.norm(sz@sx@sz+sx)),
                holding_spin_fixed_error=float(np.linalg.norm(sz@sx@sz-sx)),
                source_status='MAIN_SM_HELICITY_SIGN_WITH_PHYSICAL_SPIN_FIXED_OPEN')


REQUIRED_BY_CHANNEL={
 'finite_device':frozenset({'H_device','Pa','Pb','lead_contacts','field_protocol'}),
 'spin_phe':frozenset({'H_k','velocity_k','spin_zeeman','fermi_protocol','magnetic_units','two_color_poles','tensor_definition'}),
 'orbital_mnhe':frozenset({'H_k','velocity_k','full_wannier_position','magnetic_positional_shift','orbital_moment','interband_connections','BZ_measure','fermi_protocol','magnetic_units','tensor_definition'}),
}


def capability_gate(channel,provided:Mapping[str,object],complete_embedding=False):
    """Fail-closed source-observable input inventory, NOT solver certification.

    Even when all named pieces supplied, this only establishes an INPUT
    contract; gauge-consistency, SI and source-equation fidelity remain open.
    """
    if channel not in REQUIRED_BY_CHANNEL:raise ValueError('unknown source channel')
    if not isinstance(provided,Mapping):raise TypeError('provided must be a mapping')
    req=REQUIRED_BY_CHANNEL[channel]
    missing=sorted(k for k in req if k not in provided or provided[k] is None)
    if channel=='orbital_mnhe' and not complete_embedding:missing.append('verified_full_magnetic_Wannier_embedding')
    if missing:raise ValueError('incomplete '+channel+' input: '+', '.join(missing))
    return {'channel':channel,'input_contract':'SATISFIED_NOT_SCIENTIFIC_REPRODUCTION',
            'ready_for_author_claim':False,'proof_of_bulk_mapping':False}
