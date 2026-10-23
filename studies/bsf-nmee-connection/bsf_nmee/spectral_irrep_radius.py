"""Stage7g: exact eight-band spectral-sector tensor on two K/K' disks.

Four occupied-to-empty groups: LL, remote-occupied/low-empty (RL),
low-occupied/remote-empty (LR), and remote-to-remote (RR).
Return ALL 2x2 spin-x/y, velocity-x/y components; groups are sums of
whole 2x2 band blocks, not arbitrary individual degenerate eigenstates.

Protocol: synthetic 8-band, T=0 insulating *equilibrium static* d_B NMEE.
No 4-band downfolding, no driven-DC or real-material prediction.
"""
from dataclasses import replace
import numpy as np
from nmee_lab.clean_core.core import MU_B
from .sliding_tight_binding import SlidingBilayer, SlidingSpec

CATEGORIES=("low_to_low", "remote_occ_to_low_empty",
            "low_occ_to_remote_empty", "remote_to_remote")


def _disk(radius, nr, nt):
    if radius <= 0 or nr < 2 or nt < 4:
        raise ValueError("Bad valley disk integration inputs")
    x, gx=np.polynomial.legendre.leggauss(nr)
    rr=radius*(x+1)/2
    wr=radius*gx/2
    phi=2*np.pi*(np.arange(nt)+.5)/nt
    displacement=(rr[:,None,None]*np.stack(
        (np.cos(phi),np.sin(phi)),axis=-1)[None,:,:]).reshape(-1,2)
    weights=np.repeat(rr*wr,nt)/(2*np.pi*nt)
    return displacement, weights


def spectral_tensor(spec: SlidingSpec, *, same_scale=1., cross_scale=1.,
                    radius=.75, nr=26, nt=72, step_T=.25):
    """(4,2,2) signed spectral-sector d_B response, K+K' disks.

    A finite isolated band-group gap must exist. For sectors LL/RL/LR/RR,
    the group projection is covariant under rotations of degenerate states
    INSIDE each group. Not a standalone experimental observable.
    """
    if step_T <= 0:
        raise ValueError("positive Zeeman step required")
    model=SlidingBilayer(spec)
    model.hops=[(ia,ib,d,am*(same_scale if ia==ib else cross_scale))
                for ia,ib,d,am in model.hops]
    g=model.model().reciprocal
    centers=np.array([[1/3,2/3],[2/3,1/3]])@g
    dk,weights=_disk(radius,nr,nt)
    total=np.zeros((4,2,2))
    min_mu_gap=float("inf")
    min_group_gap=float("inf")
    spins=(model.spins['total_x'],model.spins['total_y'])
    Sz=model.spins['total_z']
    valley_parts=[]
    for center in centers:
        k=dk+center
        H=model.matrices(k)
        vel=[model.matrices(k,axis=j) for j in range(2)]
        en=np.linalg.eigvalsh(H)
        mu_gap=min(float(np.min(en[:,4]-spec.mu_eV)),
                   float(np.min(spec.mu_eV-en[:,3])))
        group_gap=min(float(np.min(en[:,2]-en[:,1])),
                      float(np.min(en[:,6]-en[:,5])))
        if mu_gap <= 0 or group_gap <= 0:
            raise ValueError("An insulating gap or chosen remote/low band group closes")
        min_mu_gap=min(min_mu_gap,mu_gap)
        min_group_gap=min(min_group_gap,group_gap)
        field=[]
        for sign in (1,-1):
            es,U=np.linalg.eigh(H-sign*model.mono.g_factor*MU_B*step_T*Sz)
            if not np.all((es[:,3]<spec.mu_eV)&(spec.mu_eV<es[:,4])):
                raise ValueError("Finite B closes chemical-potential gap")
            Ud=U.conj().swapaxes(1,2)
            Sm=[Ud@Si@U for Si in spins]
            Vm=[Ud@vj@U for vj in vel]
            den=(es[:,:4,None]-es[:,None,4:])**2
            sectors=np.empty((4,2,2,len(k)),float)
            for i in range(2):
                for j in range(2):
                    pair=-2*np.imag(Sm[i][:,:4,4:] *
                          np.swapaxes(Vm[j][:,4:,:4],1,2)/den)
                    sectors[0,i,j]=pair[:,2:4,0:2].sum((1,2))
                    sectors[1,i,j]=pair[:,0:2,0:2].sum((1,2))
                    sectors[2,i,j]=pair[:,2:4,2:4].sum((1,2))
                    sectors[3,i,j]=pair[:,0:2,2:4].sum((1,2))
            field.append(sectors)
        response=np.einsum('aijn,n->aij',(field[0]-field[1])/(2*step_T),weights)
        total+=response
        valley_parts.append(response.tolist())
    return {'sectors':total,'valley_sectors':valley_parts,
            'min_mu_gap_eV':min_mu_gap,
            'min_remote_group_gap_eV':min_group_gap}


def _irrep(m):
    return {'matrix':np.asarray(m).tolist(),
            'scalar':float(np.trace(m)/2),
            'nematic_xx':float((m[0,0]-m[1,1])/2),
            'nematic_xy':float((m[0,1]+m[1,0])/2),
            'antisymmetric_xy':float((m[0,1]-m[1,0])/2)}


def channel_irrep_sectors(spec: SlidingSpec, *, radius=.75,
                          nr=26, nt=72, step_T=.25):
    """Finite-t estimates of Gamma=A*s^2+B*c^2+C*s*c+O(t^4).

    Each coefficient is partitioned by occupied-to-empty band GROUPS.
    t must be small for the second-order coefficient approximation.
    """
    t=float(spec.tunneling_eV)
    if not np.isfinite(t) or t<=0:
        raise ValueError("positive hopping reference amplitude required")
    quadrature=dict(radius=radius,nr=nr,nt=nt,step_T=step_T)
    settings={'full':(1,1),'cross_reversed':(1,-1),
              'same_only':(1,0),'cross_only':(0,1),'zero':(0,0)}
    cases={name:spectral_tensor(spec,same_scale=s,cross_scale=c,**quadrature)
           for name,(s,c) in settings.items()}
    sect={name:data['sectors'] for name,data in cases.items()}
    coeff={'A':(sect['same_only']-sect['zero'])/t**2,
           'B':(sect['cross_only']-sect['zero'])/t**2,
           'C':(sect['full']-sect['cross_reversed'])/(2*t**2)}
    return {
        'scope':'2 K/Kprime disks, 8-band synthetic, insulating T=0 STATIC NMEE',
        'u_shift_fractional':list(spec.shift),
        't_eV':t,'radius_inv_nm':radius,'nr':nr,'nt':nt,'step_T':step_T,
        'min_mu_gap_eV':min(z['min_mu_gap_eV'] for z in cases.values()),
        'min_remote_group_gap_eV':min(z['min_remote_group_gap_eV'] for z in cases.values()),
        'categories':list(CATEGORIES),
        'channels':{name:{
            'per_sector':{cat:_irrep(m) for cat,m in zip(CATEGORIES,arr)},
            'total':_irrep(arr.sum(axis=0))}
            for name,arr in coeff.items()},
        'caveat':'spectral-sector and fixed orbital-channel decomposition, not dressed SW or driven DC',
    }
