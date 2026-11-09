"""Stage 9: explicit stacking coordinates and B-prepared electric response.

All operators and units follow Stage 8. No sliding dynamics or total structural
energy is supplied by this electronic Hamiltonian.
"""
from dataclasses import replace
import numpy as np
from nmee_lab.clean_core.core import MU_B
from .bilayer_audit import AuditBilayer


def fraction120_to60(uv):
    """A1=a1, A2=a2-a1; rows a1,a2 in the code meet at 60 degrees."""
    u,v=np.asarray(uv,float)
    return np.array([u-v,v])


def ab_ba_shift(lam):
    """One shortest AB→BA segment, length one nearest-neighbour bond.

    AB=(-1/3,2/3) and BA=(1/3,1/3) in the 60-degree basis. These are the
    literal images of (1/3,2/3)→(2/3,1/3) in the declared 120-degree basis.
    AB/BA labels depend on the top/bottom and sublattice convention.
    """
    return tuple(fraction120_to60(((1+lam)/3,(2-lam)/3)))


def half_filled_spec(spec,k,batch=1024):
    """Choose mu inside the sampled four-band gap, never change filling."""
    model=AuditBilayer(spec);vmax=-np.inf;cmin=np.inf
    for start in range(0,len(k),batch):
        e=np.linalg.eigvalsh(model.matrices(k[start:start+batch]))
        vmax=max(vmax,float(e[:,3].max()));cmin=min(cmin,float(e[:,4].min()))
    if cmin<=vmax:
        raise ValueError('no sampled half-filling insulating gap')
    return replace(spec,mu_eV=(vmax+cmin)/2)


def prepared_B_response(model,k,weights,energies_eV,*,eta_eV=1e-6,step_T=.0625,batch=1024,
                        observation='total'):
    """B equilibrated first; linear E at complex energy hbar*omega+i*eta.

    Differentiate the complete gapped linear electric susceptibility at +/-B.
    This defines an operational response, not an identification with a paper's
    intrinsic/extrinsic split. Return total in-plane spin Cij(omega).
    """
    frequencies=np.asarray(energies_eV,float)
    if frequencies.ndim!=1 or not np.isfinite(frequencies).all() or eta_eV<=0 or step_T<=0:
        raise ValueError('finite real frequency vector and positive eta/step required')
    omega=frequencies+1j*eta_eV
    if observation not in ('total','top','bottom'):
        raise ValueError('observation must be total/top/bottom; B always couples to total spin')
    total=np.zeros((len(omega),2,2),complex)
    for start in range(0,len(k),batch):
        ks=k[start:start+batch];ws=weights[start:start+batch]
        h=model.matrices(ks);vel=[model.matrices(ks,axis=j) for j in range(2)]
        for sign in (1,-1):
            e,u=np.linalg.eigh(h-sign*model.mono.g_factor*MU_B*step_T*model.spins['total_z'])
            if not np.all((e[:,3]<model.spec.mu_eV)&(model.spec.mu_eV<e[:,4])):
                raise ValueError('finite B changed the insulating occupations')
            ud=u.conj().swapaxes(1,2);delta=e[:,:4,None]-e[:,None,4:]
            for i,axis in enumerate(('x','y')):
                spin=(ud@model.spins[observation+'_'+axis]@u)[:,:4,4:]
                for j,v in enumerate(vel):
                    vv=(ud@v@u)[:,4:,:4].swapaxes(1,2);a=spin*vv
                    for index,z in enumerate(omega):
                        chi=(1j/delta*(a.conj()/(z-delta)+a/(z+delta))).sum(axis=(1,2))
                        total[index,i,j]+=sign*(ws@chi)/(2*step_T)
    return {'energies_eV':frequencies.tolist(),'eta_eV':eta_eV,'step_T':step_T,
            'real':total.real.tolist(),'imag':total.imag.tolist()}
