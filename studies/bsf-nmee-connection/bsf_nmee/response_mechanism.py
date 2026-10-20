"""Finite-B decomposition of occupied-to-empty mixed (spin, momentum) curvature.

Only the well-defined T=0, gapped, frozen-occupation STATIC susceptibility.
The observable Gamma is the finite-B derivative of chi(s_x,E_x); unlike the
reported driven-DC NMEE coefficient, it is not an author-verified transport
response. An optional *diagnostic* split into numerators versus energy
separations uses indexed instantaneous energy eigenvectors and is NOT a
basis-independent decomposition across degenerate same-occupation bands.
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from .sliding_tight_binding import SlidingBilayer, SlidingSpec
from nmee_lab.clean_core.core import MU_B


def evaluate_mechanism(spec: SlidingSpec, k, weights, *, step_T=0.25,
                       hopping_channel='full', cross_amplitude_sign=1, batch=2048):
    """Return total and optional diagnostic contributions, per supplied region.

    Finite-field kernel identity: chi = -2 Im sum_{v occupied,c empty}
    Sx_{vc} vx_{cv} / (E_v-E_c)^2.

    At each k we *index* sorted band eigenpairs at ±B (nonunique if same-side
    degeneracies); the symmetric product identity splits d_B chi exactly.
    The SUM is physical/gauge-invariant, but each diagnostic contribution may
    depend on the choice of eigenvectors within degenerate band manifolds.
    """
    if hopping_channel not in ('full','same_only','cross_only'):
        raise ValueError('bad hopping channel')
    if step_T <= 0: raise ValueError('B step must be positive')
    if cross_amplitude_sign not in (1,-1): raise ValueError('cross sign must be +1 or -1')
    obj=SlidingBilayer(spec)
    if hopping_channel!='full':
        obj.hops=[hop for hop in obj.hops if (hop[0]==hop[1])==(hopping_channel=='same_only')]
    if cross_amplitude_sign==-1:
        obj.hops=[(ia,ib,disp,(-amp if ia!=ib else amp)) for ia,ib,disp,amp in obj.hops]
    k=np.asarray(k,float);w=np.asarray(weights,float)
    if k.shape != (len(w),2):raise ValueError('k/weight shapes mismatch')
    so=obj.spins['total_x'];sz=obj.spins['total_z']
    gamma=0.;wave=0.;spectrum=0.;baregap=np.inf;mu_gap=np.inf
    for i in range(0,len(k),batch):
        kp=k[i:i+batch];wp=w[i:i+batch]
        H=obj.matrices(kp);V=obj.matrices(kp,axis=0)
        h0=np.linalg.eigvalsh(H)
        if not np.all((h0[:,3]<spec.mu_eV)&(spec.mu_eV<h0[:,4])):
            raise ValueError('base Hamiltonian not insulating at fixed mu')
        baregap=min(baregap,float(np.min(h0[:,4]-h0[:,3])))
        mu_gap=min(mu_gap,float(np.min(h0[:,4]-spec.mu_eV)),float(np.min(spec.mu_eV-h0[:,3])))
        numer=[];denom=[]
        for sign in (1,-1):
            energies,U=np.linalg.eigh(H-sign*obj.mono.g_factor*MU_B*step_T*sz)
            if not np.all((energies[:,3]<spec.mu_eV)&(spec.mu_eV<energies[:,4])):
                raise ValueError('field closes the insulating Fermi gap')
            Uh=U.conj().swapaxes(1,2)
            S=Uh@so@U;Vm=Uh@V@U
            n=-2*np.imag(S[:,:4,4:]*np.swapaxes(Vm[:,4:,:4],1,2))
            d=(energies[:,:4,None]-energies[:,None,4:])**-2
            numer.append(n);denom.append(d)
        dp=denom[0];dm=denom[1];np_=numer[0];nm=numer[1]
        total=(np_*dp-nm*dm)/(2*step_T)
        part_wave=(np_-nm)*(dp+dm)/(4*step_T)
        part_spec=(np_+nm)*(dp-dm)/(4*step_T)
        gamma+=float(np.dot(wp,np.sum(total,axis=(1,2))))
        wave+=float(np.dot(wp,np.sum(part_wave,axis=(1,2))))
        spectrum+=float(np.dot(wp,np.sum(part_spec,axis=(1,2))))
    return {
      'gamma':gamma,'matrix_element_diagnostic':wave,
      'energy_denominator_diagnostic':spectrum,
      'decomposition_closure_error':gamma-wave-spectrum,
      'sampled_min_middle_gap_eV':baregap,
      'sampled_min_mu_distance_eV':mu_gap,
      'method':'exact static mixed-spin/momentum Berry-curvature finite-field difference; split is nonunique under intrablock degeneracy',
      'hopping_channel':hopping_channel,'cross_amplitude_sign':cross_amplitude_sign}
