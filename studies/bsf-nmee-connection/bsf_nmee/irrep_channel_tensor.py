"""Stage 7f: analytic weak-hopping *full* in-plane NMEE tensors.

A/B/C refer to same^2, cross^2 and same*cross hopping amplitudes in the
fixed source-orbital embedding. Unlike the velocity-path split, their total
2x2 tensors are directly the derivatives of the physical static response.
The high-symmetry reference has C3 and layer-exchange mirror mz.

Protocol is the synthetic 8-band, insulating, T=0, equilibrium/static
finite-Bz response, NOT the disputed driven-DC NMEE or a material.
"""
from dataclasses import replace
import numpy as np
from .sliding_tight_binding import SlidingBilayer, SlidingSpec
from .projector_jet_response import _polar_nodes, projector_coefficient


def analytic_tensor_ABC(spec=SlidingSpec(), *, radius=.75, nr=14, nt=48):
    reference=SlidingBilayer(replace(spec,tunneling_eV=0.))
    same=SlidingBilayer(replace(spec,tunneling_eV=1.,
                                cross_sublattice=False))
    both=SlidingBilayer(replace(spec,tunneling_eV=1.,
                                cross_sublattice=True))
    G=reference.model().reciprocal
    centers=np.array([[1/3,2/3],[2/3,1/3]])@G
    q,weights=_polar_nodes(radius,nr,nt)
    result={name:np.zeros((2,2),float) for name in ('A','B','C')}
    per_valley=[]
    for center in centers:
        k=q+center
        H0=reference.matrices(k)
        Hs=same.matrices(k)
        Hf=both.matrices(k)
        Vs=Hs-H0
        Vc=Hf-Hs
        vel0=[reference.matrices(k,axis=j) for j in range(2)]
        Ws=[same.matrices(k,axis=j)-vel0[j] for j in range(2)]
        Wc=[both.matrices(k,axis=j)-same.matrices(k,axis=j)
             for j in range(2)]
        local={name:np.zeros((2,2),float) for name in result}
        for i,spin in enumerate((reference.spins['total_x'],
                                 reference.spins['total_y'])):
            for j in range(2):
                for name,X,Y,Wx,Wy,multiplier in (
                    ('A',Vs,Vs,Ws[j],Ws[j],.5),
                    ('B',Vc,Vc,Wc[j],Wc[j],.5),
                    ('C',Vs,Vc,Ws[j],Wc[j],1.)):
                    kernel=projector_coefficient(
                        H0,X,Y,vel0[j],Wx,Wy,spin,
                        reference.spins['total_z'],spec.mu_eV,
                        g_factor=reference.mono.g_factor,part='all')
                    local[name][i,j]=multiplier*float(weights@kernel)
                    result[name][i,j]+=local[name][i,j]
        per_valley.append({k:v.tolist() for k,v in local.items()})
    tensors={}
    for name,M in result.items():
        scalar=np.trace(M)/2
        nematic=.5*(M+M.T)-scalar*np.eye(2)
        anti=.5*(M-M.T)
        tensors[name]={
            'matrix':M.tolist(),
            'trace':float(2*scalar),
            'scalar':float(scalar),
            'nematic':nematic.tolist(),
            'nematic_norm':float(np.linalg.norm(nematic)),
            'antisymmetric':anti.tolist(),
        }
    return {
        'scope':'two K/Kprime disks, synthetic 8-band, 0 K gapped static d_B',
        'shift_fractional':list(spec.shift),
        'radius_inv_nm':radius,'nr':nr,'nt':nt,
        'tensors':tensors,'per_valley':per_valley,
        'caveat':'channels are defined in the declared A/B orbital basis',
    }


def rotate_slide(spec: SlidingSpec, degrees):
    """Rotate the actual Cartesian stacking displacement, not fractional axes."""
    a=SlidingBilayer(spec).a
    theta=np.deg2rad(degrees)
    R=np.array([[np.cos(theta),-np.sin(theta)],
                [np.sin(theta), np.cos(theta)]])
    cart=np.asarray(spec.shift)@a
    frac=(R@cart)@np.linalg.inv(a)
    return replace(spec,shift=tuple(float(v) for v in frac))


def tensor_covariance_residual(base, target, degrees, *, odd=False):
    theta=np.deg2rad(degrees)
    R=np.array([[np.cos(theta),-np.sin(theta)],
                [np.sin(theta), np.cos(theta)]])
    e=0
    for name in ('A','B','C'):
        M=np.asarray(base['tensors'][name]['matrix'])
        N=np.asarray(target['tensors'][name]['matrix'])
        expectation=(-1 if odd else 1)*(R@M@R.T)
        e=max(e,float(np.max(np.abs(N-expectation))))
    return e
