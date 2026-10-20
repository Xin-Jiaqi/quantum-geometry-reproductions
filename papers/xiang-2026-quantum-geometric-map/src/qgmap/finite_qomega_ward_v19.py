"""v19 independent SCBA finite-(q,omega) resolvent Ward check.

Independent periodic two-band Gaussian scalar disorder model from v17.
R/A propagators use energies mu +- omega/2 and regulators +-i*eta,
with external momentum k +- q/2. Self-energies are determined separately
by self-consistent Born equations at both complex energies. The density
and longitudinal finite-difference current vertices obey noncrossing
Bethe-Salpeter ladders with the SAME disorder variance W.

The audited matrix identity is
 (omega + 2i eta) Gamma^0 - q Gamma^L
   = GR(k+q/2)^(-1) - GA(k-q/2)^(-1)
   = (omega+2i eta)I - (H_+-H_-) -(Sigma_R-Sigma_A).

The 2i eta term is an external resolvent regulator. This finite-eta identity
is NOT by itself the real-frequency physical continuity equation, nor a
calculation of an experimental dc PHE. The eta->0, q->0, omega->0 limits
need further collision/lead treatment.
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from .warped_disorder_v17 import Lattice,k_grid,h_and_vel
from .microscopic_ward_v11 import I,scba,ladder_vertex


def two_energy_greens(p:Lattice,*,qx=0.,qy=0.,omega=.04,b=0.,protocol='onsite'):
    if not np.isfinite([qx,qy,omega,b]).all():
        raise ValueError('External perturbations must be finite')
    k=k_grid(p)
    q=np.array([qx,qy],float)
    kp=k+q/2
    km=k-q/2
    hp,_,_=h_and_vel(kp,p,b,protocol)
    hm,_,_=h_and_vel(km,p,b,protocol)
    sr,gr,_,nr=scba(hp,replace(p,chemical_potential=p.chemical_potential+omega/2))
    sa_ret,_,ga,na=scba(hm,replace(p,chemical_potential=p.chemical_potential-omega/2))
    sa=sa_ret.conj().T
    return dict(hp=hp,hm=hm,gr=gr,ga=ga,sr=sr,sa=sa,
                regulator=complex(omega,2*p.eta), iterations=(nr,na),q=q)


def finite_ward(p:Lattice, *, q=.08, omega=.04, direction=1,
                b=0.,protocol='onsite',calc_bare=True):
    if direction not in (0,1):raise ValueError('direction must be 0 or 1')
    if not np.isfinite(q) or abs(q)>=np.pi:raise ValueError('q finite, |q|<pi')
    vals=dict(qx=q if direction==0 else 0.,qy=q if direction==1 else 0.,omega=omega,b=b,protocol=protocol)
    x=two_energy_greens(p,**vals)
    gr,ga=x['gr'],x['ga']
    bare_density=np.broadcast_to(I,x['hp'].shape)
    gam0,c0,r0=ladder_vertex(gr,ga,bare_density,p.disorder_variance)
    if q:
        bare_L=(x['hp']-x['hm'])/q
        gamL,cL,rL=ladder_vertex(gr,ga,bare_L,p.disorder_variance)
    else:
        gamL=np.zeros_like(gr);cL=np.zeros((2,2),complex);rL=0.
        bare_L=np.zeros_like(gr)
    z=x['regulator'];dsig=x['sr']-x['sa']
    left=z*gam0-q*gamL
    right=z*I-(x['hp']-x['hm'])-dsig
    gr_inv=np.linalg.inv(gr)
    ga_inv=np.linalg.inv(ga)
    independent=gr_inv-ga_inv
    defects=dict(
        physical_frequency_without_regulator=float(np.max(abs(omega*gam0-q*gamL-right))),
        missing_regulator_expected=float(2*p.eta*np.max(abs(gam0))),
        wti=float(np.max(abs(left-right))),
        resolvent=float(np.max(abs(independent-right))),
        ladder=float(max(r0,rL)),
        selfenergy=float(max(np.max(abs(x['sr']-p.disorder_variance*gr.mean((0,1)))),
                              np.max(abs(x['sa']-p.disorder_variance*ga.mean((0,1)))))),
        bare=float(np.max(abs((z*bare_density-q*bare_L)-right)) if calc_bare else None,
    )
    out=dict(protocol='independent_periodic_SCBA_finite_qomega_resolvent_WTI',
             q=float(q),omega=float(omega),eta=float(p.eta),W=float(p.disorder_variance),
             direction='x' if direction==0 else 'y',
             field_protocol=protocol,b=float(b),
             residuals=defects,
             density_ladder_correction_norm=float(np.max(abs(c0))),
             longitudinal_ladder_correction_norm=float(np.max(abs(cL))),
             density_vertex_norm=float(np.max(abs(gam0))),
             mean_density_vertex_trace=(np.trace(gam0,axis1=-2,axis2=-1).mean()/2),
             iterations=x['iterations'],
             warnings=['finite eta identity includes 2i eta external regulator',
                       'noncrossing scalar disorder, NOT physical current Ward at eta>0',
                       'not original Xiang Liouville/RTA PHE'])
    return out
