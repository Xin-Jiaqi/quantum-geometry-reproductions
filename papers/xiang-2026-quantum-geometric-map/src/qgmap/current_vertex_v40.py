"""v40: same-disorder RA/AR/RR/AA current ladder and its Hall zero controls.

A physically separated, noncrossing, q=0, fixed-energy Bethe-Salpeter
*method* for the independent QWZ model (not the paper's EB DC response).

At fixed complex-energy propagators Gs(k), Gt(k) and scalar white-noise
variance W, solve exactly in the four-Pauli basis:
  Gamma_a^{st}(k)=v_a(k)+C_a^{st},
  C_a^{st}=W <Gs(k) [v_a(k)+C_a^{st}] Gt(k)>_BZ.

The RR zero-current correction follows periodic momentum covariance with
momentum-independent Sigma, whereas the RA current ladder may be nonzero.
The *isolated* antisymmetric RA bubble can be nonzero even when the fully
symmetrized spectral RA+AR-RR-AA sector is antisymmetry-free. Do NOT read
bare or dressed RA antisymmetry as a physical Hall conductivity.

Normalized BZ means, dimensionless e=hbar=1. External eta is a regulator;
no energy-dependent Streda I/surface term, crossed diagrams, or nonlinear
mixed EB DC response is implemented by this module.
"""
from __future__ import annotations
import numpy as np
from .streda_clean_sea_v38 import qwz_model, I, SX, SY, SZ
from .scba_streda_v39 import SCBAConfig, solve_self_energy

PAULI=np.stack((I,SX,SY,SZ))


def pauli_vertex(G_left,G_right,v,W):
    """Independent 4x4 Bethe-Salpeter solve; returns matrices and diagnostics."""
    G_left=np.asarray(G_left,complex);G_right=np.asarray(G_right,complex)
    v=np.asarray(v,complex)
    if G_left.ndim!=4 or G_left.shape[-2:]!=(2,2) or G_right.shape!=G_left.shape or v.shape!=G_left.shape:
        raise ValueError('G_left,G_right,v need matching (Nx,Ny,2,2) shapes')
    if not np.isfinite(G_left).all() or not np.isfinite(G_right).all() or not np.isfinite(v).all():
        raise ValueError('All propagators and vertices must be finite')
    if W<0 or not np.isfinite(W):raise ValueError('W nonnegative and finite')
    phi=np.mean(G_left@v@G_right,axis=(0,1))
    basis_rhs=np.array([np.trace(t@phi)*.5 for t in PAULI])
    M=np.column_stack([
        np.array([np.trace(t@(G_left@s@G_right).mean((0,1)))*.5 for t in PAULI])
        for s in PAULI
    ])
    linear=np.eye(4,dtype=complex)-W*M
    cond=float(np.linalg.cond(linear))
    if not np.isfinite(cond) or cond>1e12:
        raise ValueError('Ladder near a singular diffuson; explicit limit analysis needed')
    coeff=np.linalg.solve(linear,W*basis_rhs)
    constant=np.einsum('a,aij->ij',coeff,PAULI)
    dressed=v+constant
    residual=float(np.max(np.abs(constant-W*np.mean(G_left@dressed@G_right,axis=(0,1)))))
    return dict(dressed=dressed,constant=constant,coeff=coeff,residual=residual,
                condition=cond,max_constant=float(np.max(abs(constant))))


def _bubble(vx,gl,dressed_y,gt):
    return complex(np.mean(np.trace(vx@gl@dressed_y@gt,axis1=-2,axis2=-1)))


def _as_pair(v):
    return {'real':float(np.real(v)),'imag':float(np.imag(v))}


def current_audit(mesh=21,mass=-1.,W=.12,energy=1.1,eta=.12,offset=(0.,0.)):
    """SCBA+current ladder in four channels and separate density Ward checks."""
    cfg=SCBAConfig(mesh=mesh,mass=mass,disorder_variance=W,eta=eta)
    H,vx,vy=qwz_model(mesh,mass,offset=offset)
    solved=solve_self_energy(H,energy,cfg)
    GR=solved['G'];GA=GR.conj().swapaxes(-1,-2)
    conditions=[];residuals=[];channels={};constants={}
    for name,L,R in (('ra',GR,GA),('ar',GA,GR),('rr',GR,GR),('aa',GA,GA)):
        cx=pauli_vertex(L,R,vx,W)
        cy=pauli_vertex(L,R,vy,W)
        conditions.extend([cx['condition'],cy['condition']]);residuals.extend([cx['residual'],cy['residual']])
        channels[name]={
            'xy':_bubble(vx,L,cy['dressed'],R),
            'yx':_bubble(vy,L,cx['dressed'],R),
            'bare_xy':_bubble(vx,L,vy,R),
            'bare_yx':_bubble(vy,L,vx,R)
        }
        constants[name]={'x':cx['max_constant'],'y':cy['max_constant']}
    assembled={a:(channels['ra'][a]+channels['ar'][a]-channels['rr'][a]-channels['aa'][a])/4
               for a in ('xy','yx')}
    spectral=(GA-GR)/(2j)
    direct_xy=_bubble(vx,spectral,vy,spectral)
    direct_yx=_bubble(vy,spectral,vx,spectral)
    # q=0 RA density Ward using exactly the same SCBA W and eta.
    density=pauli_vertex(GR,GA,np.broadcast_to(I,GR.shape),W)
    gamma_from_sigma=I-(solved['sigma']-solved['sigma'].conj().T)/(2j*eta)
    return dict(protocol='noncrossing_QWZ_current_ladder_four_channel_Hall_zero_control',
                mesh=mesh,mass=mass,W=W,energy=energy,eta=eta,offset=list(offset),
                selfenergy_residual=solved['residual'],selfenergy_iterations=solved['iterations'],
                max_bethe_salpeter_residual=max(residuals),max_linear_condition=max(conditions),
                current_correction_norms=constants,
                channel_bubbles={c:{key:_as_pair(val) for key,val in v.items()} for c,v in channels.items()},
                ra_antisymmetric_dressed=float(np.real((channels['ra']['xy']-channels['ra']['yx'])/2)),
                ra_antisymmetric_bare=float(np.real((channels['ra']['bare_xy']-channels['ra']['bare_yx'])/2)),
                ra_vertex_effect=float(np.real(((channels['ra']['xy']-channels['ra']['yx'])-
                    (channels['ra']['bare_xy']-channels['ra']['bare_yx']))/2)),
                four_channel_spectral_xy=_as_pair(assembled['xy']),
                four_channel_spectral_yx=_as_pair(assembled['yx']),
                four_channel_hall_antisym_abs=float(abs((assembled['xy']-assembled['yx'])/2)),
                bare_spectral_hall_antisym_abs=float(abs((direct_xy-direct_yx)/2)),
                ra_density_ward_residual=float(np.max(abs(density['constant']+I-gamma_from_sigma))),
                ra_density_vertex_norm=float(np.max(abs(gamma_from_sigma))),
                warnings=['RA antisymmetric bubble is NOT a physical Hall conductivity',
                          'four-channel spectral combination is symmetric in current indices',
                          'RR current Ward uses periodic momentum covariance and momentum-independent Sigma',
                          'full disorder Hall Streda/Bastin and original nonlinear EB DC remain OPEN'])


def exact_disorder_ra_torus(L=4,mass=-1.,W=.12,mu=1.1,eta=.12,seeds=8):
    """Exact finite disordered torus bare RA cross bubble, eigenbasis computed.

    This finite-size/finite-sample Kubo integrand has NO averaged vertex;
    compare to clean QWZ with offset=(0,0) at W=0, but NEVER claim
    generic agreement with infinite-system SCBA at nonzero W.
    """
    if L<3 or L>9 or seeds<1 or W<0 or eta<=0:raise ValueError('invalid torus settings')
    k=-np.pi+2*np.pi*np.arange(L)/L
    x,y=np.meshgrid(k,k,indexing='ij')
    d=np.stack((np.sin(x),np.sin(y),mass+np.cos(x)+np.cos(y)))
    dx=np.stack((np.cos(x),np.zeros_like(x),-np.sin(x)))
    dy=np.stack((np.zeros_like(y),np.cos(y),-np.sin(y)))
    spin=np.stack((SX,SY,SZ))
    H=np.einsum('a...,aij->...ij',d,spin)
    vx=np.einsum('a...,aij->...ij',dx,spin)
    vy=np.einsum('a...,aij->...ij',dy,spin)
    nk=L*L
    pts=np.stack(np.meshgrid(-np.pi+2*np.pi*np.arange(L)/L,
        -np.pi+2*np.pi*np.arange(L)/L,indexing='ij'),axis=-1).reshape((-1,2))
    sites=np.stack(np.meshgrid(np.arange(L),np.arange(L),indexing='ij'),axis=-1).reshape((-1,2))
    U=np.kron(np.exp(1j*sites@pts.T)/np.sqrt(nk),I)
    def site_operator(t):
        block=np.zeros((2*nk,2*nk),complex)
        for j,val in enumerate(t.reshape((nk,2,2))):block[2*j:2*j+2,2*j:2*j+2]=val
        return U@block@U.conj().T
    bareH=site_operator(H); jx=site_operator(vx);jy=site_operator(vy)
    values=[]
    for seed in range(seeds):
        rng=np.random.default_rng(38200+seed)
        onsites=np.repeat(rng.normal(scale=np.sqrt(W),size=nk),2)
        h=bareH+np.diag(onsites)
        es,evec=np.linalg.eigh(h)
        rot=evec.conj().T
        xx=rot@jx@evec;yy=rot@jy@evec
        gr=1/(mu+1j*eta-es)
        ga=gr.conj()
        xy=np.sum(xx*yy.T*(gr[None,:]*ga[:,None]))/nk
        yx=np.sum(yy*xx.T*(gr[None,:]*ga[:,None]))/nk
        values.append(float(np.real((xy-yx)/2)))
    return dict(protocol='finite_torus_exact_disorder_bare_RA_antisymmetric_bubble_only',
                L=L,mass=mass,W=W,mu=mu,eta=eta,seeds=seeds,
                mean=float(np.mean(values)),standard_deviation=float(np.std(values)),
                samples=values,
                caveat='W=0 agrees with same-mesh clean QWZ; W>0 finite ensembles are not SCBA and do not include a transport current vertex')
