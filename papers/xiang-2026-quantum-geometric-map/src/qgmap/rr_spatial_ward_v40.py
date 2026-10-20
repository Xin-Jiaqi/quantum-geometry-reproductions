"""v40b: independent uniform-twist finite-difference RR spatial Ward audit.

Use the SAME QWZ scalar-disorder self-consistency, but obtain dSigma/dA and
dGR/dA by solving SCBA again for H(k + A). Compare to independently solved
four-Pauli current ladder C_RR, Gamma_a = v_a+C_RR.

Identity from differentiating Sigma = W <(E+i eta-H(k+A)-Sigma)^-1>:
    dSigma/dA_a = W < GR (v_a+dSigma/dA_a) GR > = C_a^RR.
This is a *spatial* q->0 same-sector (RR) Ward check at finite eta,
not a physical disordered Hall transport or the original mixed EB DC PHE.
"""
from __future__ import annotations
import numpy as np
from .streda_clean_sea_v38 import qwz_model
from .scba_streda_v39 import SCBAConfig,solve_self_energy
from .current_vertex_v40 import pauli_vertex


def twisted_qwz(mesh,mass,ax=0.,ay=0.,offset=(0.,0.)):
    """Uniform lattice vector potential via a full periodic momentum shift."""
    if len(offset)!=2 or not np.isfinite([ax,ay,*offset]).all():
        raise ValueError('finite vector-potential and 2D momentum-offset required')
    return qwz_model(mesh,mass,offset=(
        offset[0]+mesh*ax/(2*np.pi),offset[1]+mesh*ay/(2*np.pi)))


def audit(mesh=21,mass=-1.,W=.12,energy=1.1,eta=.12,step=2e-4,offset=(0.,0.)):
    if not np.isfinite(step) or step<=0 or step>=.1:
        raise ValueError('Positive small finite-difference step is required')
    cfg=SCBAConfig(mesh=mesh,mass=mass,disorder_variance=W,eta=eta)
    H,vx,vy=twisted_qwz(mesh,mass,offset=offset)
    center=solve_self_energy(H,energy,cfg)
    GR=center['G']
    GA=GR.conj().swapaxes(-1,-2)
    details={}
    for axis,v in (('x',vx),('y',vy)):
        rr=pauli_vertex(GR,GR,v,W)
        ra=pauli_vertex(GR,GA,v,W)
        hplus=twisted_qwz(mesh,mass,step if axis=='x' else 0.,
                          step if axis=='y' else 0.,offset=offset)[0]
        hminus=twisted_qwz(mesh,mass,-step if axis=='x' else 0.,
                           -step if axis=='y' else 0.,offset=offset)[0]
        plus=solve_self_energy(hplus,energy,cfg,initial=center['sigma'])
        minus=solve_self_energy(hminus,energy,cfg,initial=center['sigma'])
        dsigma=(plus['sigma']-minus['sigma'])/(2*step)
        dG=(plus['G']-minus['G'])/(2*step)
        dH=(hplus-hminus)/(2*step)
        rr_g=GR@rr['dressed']@GR
        details[axis]={
            'rr_ladder_residual':rr['residual'],
            'ra_ladder_residual':ra['residual'],
            'twisted_hamiltonian_velocity_error':float(np.max(abs(dH-v))),
            'dSigma_twist_vs_RR_ladder_error':float(np.max(abs(dsigma-rr['constant']))),
            'dGR_twist_vs_RR_vertex_error':float(np.max(abs(dG-rr_g))),
            'RR_current_correction_norm':rr['max_constant'],
            'RA_current_correction_norm':ra['max_constant'],
            'RR_condition':rr['condition'],
            'RA_condition':ra['condition'],
            'selfenergy_residual_plus':plus['residual'],
            'selfenergy_residual_minus':minus['residual']
        }
    return {
        'schema':'qgmap-v40b-rr-spatial-ward',
        'science_status':'PARTIAL',
        'protocol':'independent_twisted_SCBA_vs_RR_spatial_current_Ward',
        'parameters':{'mesh':mesh,'mass':mass,'W':W,'energy':energy,'eta':eta,'step':step,'offset':list(offset)},
        'selfenergy_residual_center':center['residual'],
        'directions':details,
        'limitations':['same-sector RR spatial identity, NOT the RA Hall current',
                       'finite mesh and finite external eta',
                       'full physical Kubo-Streda and original mixed EB DC not implemented']
    }
