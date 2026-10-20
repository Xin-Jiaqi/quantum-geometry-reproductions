"""Static in-plane spin-electric NMEE tensor at B_z, for symmetry diagnostics.

Outputs [spin=x/y,electric=x/y]. This is a gapped *equilibrium static*
finite-Zeeman difference, not the PRL's unresolved driven DC definition.
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from .sliding_tight_binding import SlidingBilayer, SlidingSpec
from nmee_lab.clean_core.core import MU_B


def static_inplane_matrix(spec: SlidingSpec, k, weights, *, same_scale=1.,
                          cross_scale=1., step_T=.25, batch=2048):
    k=np.asarray(k,float);weights=np.asarray(weights,float)
    if k.ndim!=2 or k.shape[1]!=2 or weights.shape!=(len(k),):
        raise ValueError('k/weights shape mismatch')
    if not np.isfinite(weights).all() or not np.isfinite(k).all() or step_T<=0:
        raise ValueError('invalid finite inputs')
    obj=SlidingBilayer(spec)
    obj.hops=[(ia,ib,disp,amp*(same_scale if ia==ib else cross_scale))
              for ia,ib,disp,amp in obj.hops]
    S=[obj.spins['total_x'],obj.spins['total_y']]
    Sz=obj.spins['total_z']
    accumulated=np.zeros((2,2),float)
    min_distance=float('inf')
    for i in range(0,len(k),batch):
        kp=k[i:i+batch];w=weights[i:i+batch]
        H=obj.matrices(kp)
        e0=np.linalg.eigvalsh(H)
        if not np.all((e0[:,3]<spec.mu_eV)&(spec.mu_eV<e0[:,4])):
            raise ValueError('not an insulator at chemical potential')
        min_distance=min(min_distance,float(np.min(e0[:,4]-spec.mu_eV)),
                         float(np.min(spec.mu_eV-e0[:,3])))
        V=[obj.matrices(kp,axis=j) for j in range(2)]
        cur=[]
        for sign in (+1,-1):
            e,U=np.linalg.eigh(H-sign*obj.mono.g_factor*MU_B*step_T*Sz)
            if not np.all((e[:,3]<spec.mu_eV)&(spec.mu_eV<e[:,4])):
                raise ValueError('field changes occupations')
            Udag=U.conj().swapaxes(1,2)
            Sm=[(Udag@s@U)[:,:4,4:] for s in S]
            Vm=[np.swapaxes((Udag@v@U)[:,4:,:4],1,2) for v in V]
            den=(e[:,:4,None]-e[:,None,4:])**2
            cur.append(np.array([[-2*np.imag(sm*vm/den).sum(axis=(1,2))
                                  for vm in Vm] for sm in Sm]))
        derivative=(cur[0]-cur[1])/(2*step_T)
        accumulated+=np.einsum('an,n->a', derivative.reshape(4,-1),w).reshape(2,2)
    return {'gamma':accumulated,'mu_distance':min_distance,
            'protocol':'gapped T=0 static d_B mixed curvature, Bz; alpha=x/y beta=x/y'}


def split_parity(spec, k, weights, **kwargs):
    plus=static_inplane_matrix(spec,k,weights,**kwargs)
    minus=static_inplane_matrix(spec,k,weights,cross_scale=-1.,**kwargs)
    return {'full':plus['gamma'],
            'even':(plus['gamma']+minus['gamma'])/2,
            'mixed_odd':(plus['gamma']-minus['gamma'])/2,
            'mu_distance':min(plus['mu_distance'],minus['mu_distance'])}


def project_interference_tensor(tensor, rotations, characters_same=None,
                                characters_cross=None):
    """Reynolds projector for a *fixed stacking* and 1D hopping irreps.

    Group must consist of a complete finite unitary point co-group with the
    supplied physically verified characters chi_s, chi_c = +/-1; rotated
    orbitals must not mix the two declared channels. For a channel-mixed
    representation or for an operation mapping to a different stacking, do
    NOT call this function. The two axial transformation determinants cancel.
    """
    T=np.asarray(tensor,float)
    Rs=[np.asarray(R,float) for R in rotations]
    if T.shape!=(3,3,3) or not Rs:
        raise ValueError('tensor (3,3,3) and nonempty group required')
    cs=[1]*len(Rs) if characters_same is None else list(characters_same)
    cc=[1]*len(Rs) if characters_cross is None else list(characters_cross)
    if len(cs)!=len(Rs) or len(cc)!=len(Rs) or any(x not in (-1,1) for x in cs+cc):
        raise ValueError('characters must be +/-1 for every operation')
    out=np.zeros((3,3,3),float)
    for R,a,b in zip(Rs,cs,cc):
        if R.shape!=(3,3) or not np.allclose(R@R.T,np.eye(3),atol=1e-9):
            raise ValueError('operations must be Cartesian orthogonal')
        out+=a*b*np.einsum('ia,jb,kc,abc->ijk',R,R,R,T,optimize=True)
    return out/len(Rs)
