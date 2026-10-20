"""Independent two-channel amplitude extension, WITHOUT modifying frozen stage-4 code.

This code evaluates the occupied-projector static spin-electric susceptibility
at Zeeman fields ±B.  It avoids nondegenerate-band sums inside the occupied
and empty subspaces and checks the fixed 0 K Fermi gap explicitly.

s, c are real dimensionless scales of AA/BB and AB/BA interlayer hopping in a
fixed orbital/Wannier gauge.  Not a gauge-invariant splitting of an observable.
"""
from __future__ import annotations
import numpy as np
from .sliding_tight_binding import SlidingBilayer, SlidingSpec
from nmee_lab.clean_core.core import MU_B


def evaluate_scaled_channels(spec: SlidingSpec, k, weights, *, same_scale=1.0,
                             cross_scale=1.0, step_T=.25, batch=2048):
    if not np.isfinite([same_scale,cross_scale,step_T]).all() or step_T<=0:
        raise ValueError('nonfinite channel scales or nonpositive B step')
    k=np.asarray(k,float);w=np.asarray(weights,float)
    if k.ndim!=2 or k.shape[1]!=2 or w.shape!=(len(k),):
        raise ValueError('k and weights must be (N,2) and (N,)')
    obj=SlidingBilayer(spec)
    obj.hops=[(a,b,dr,h*(same_scale if a==b else cross_scale))
              for a,b,dr,h in obj.hops]
    sx=obj.spins['total_x'];sz=obj.spins['total_z']
    result=0.;gap=float('inf');mu_gap=float('inf')
    for first in range(0,len(k),batch):
        sub=k[first:first+batch];wp=w[first:first+batch]
        H=obj.matrices(sub);vx=obj.matrices(sub,axis=0)
        e0=np.linalg.eigvalsh(H)
        if not np.all((e0[:,3]<spec.mu_eV)&(spec.mu_eV<e0[:,4])):
            raise ValueError('base state not insulating at fixed mu')
        gap=min(gap,float(np.min(e0[:,4]-e0[:,3])))
        mu_gap=min(mu_gap,float(np.min(e0[:,4]-spec.mu_eV)),
                   float(np.min(spec.mu_eV-e0[:,3])))
        chis=[]
        for sign in (+1,-1):
            e,U=np.linalg.eigh(H-sign*obj.mono.g_factor*MU_B*step_T*sz)
            if not np.all((e[:,3]<spec.mu_eV)&(spec.mu_eV<e[:,4])):
                raise ValueError('Zeeman perturbed state is not insulating')
            udag=U.conj().swapaxes(1,2)
            sigma=udag@sx@U;vel=udag@vx@U
            num=sigma[:,:4,4:]*np.swapaxes(vel[:,4:,:4],1,2)
            den=(e[:,:4,None]-e[:,None,4:])**2
            chis.append(-2*np.sum(np.imag(num)/den,axis=(1,2)))
        result+=float(wp@((chis[0]-chis[1])/(2*step_T)))
    return {'gamma':result,'sampled_min_middle_gap_eV':gap,
            'sampled_min_mu_distance_eV':mu_gap,'same_scale':same_scale,
            'cross_scale':cross_scale,
            'protocol':'T=0 gapped static mixed curvature d_B; frozen occupations'}
