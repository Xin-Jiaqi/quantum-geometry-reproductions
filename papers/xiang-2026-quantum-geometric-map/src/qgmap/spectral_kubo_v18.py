"""v18: zero-temperature symmetric spectral Kubo-Greenwood sector for a periodic two-band disordered model.

This is a fixed-energy, finite eta calculation within SCBA; no claim of the
published paper's finite-T/two-color response or complete disorder Kubo-Streda
(Hall sea term/crossed diagrams). Symmetric conductivity can be assembled from
RA, AR, RR, AA channels; including only RA is NOT a complete dissipative
conductivity. Consistent noncrossing ladder is applied separately to each
channel, in the same Gaussian white-noise model as v17.

Normalize integral by BZ average; common physical e^2*hbar/pi and k measure
omitted for protocol contrast. Field b enters as an energy (not Tesla).
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from .warped_disorder_v17 import Lattice, k_grid, h_and_vel
from .microscopic_ward_v11 import scba, ladder_vertex


def spectral_bubble(p: Lattice, b:float=0.0, protocol:str='onsite', vertex:bool=True):
    k=k_grid(p)
    h,vx,vy=h_and_vel(k,p,b,protocol)
    sigma,gr,ga,it=scba(h,p)
    imgr=(gr-ga)/(2j)
    bare_xy=np.mean(np.trace(vx@imgr@vy@imgr,axis1=-2,axis2=-1))
    bare_yx=np.mean(np.trace(vy@imgr@vx@imgr,axis1=-2,axis2=-1))
    if vertex:
        s={'r':gr,'a':ga}; terms={};res=[]
        for left,right in [('r','a'),('a','r'),('r','r'),('a','a')]:
            gl,gt=s[left],s[right]
            gy,_,ry=ladder_vertex(gl,gt,vy,p.disorder_variance)
            gx,_,rx=ladder_vertex(gl,gt,vx,p.disorder_variance)
            res.extend([ry,rx])
            terms[left+right]=(
              np.mean(np.trace(vx@gl@gy@gt,axis1=-2,axis2=-1)),
              np.mean(np.trace(vy@gl@gx@gt,axis1=-2,axis2=-1))
            )
        xy=(terms['ra'][0]+terms['ar'][0]-terms['rr'][0]-terms['aa'][0])/4
        yx=(terms['ra'][1]+terms['ar'][1]-terms['rr'][1]-terms['aa'][1])/4
        channels={key:[float(np.real(v)) for v in pair] for key,pair in terms.items()}
    else:
        xy,yx=bare_xy,bare_yx
        channels=None;res=[]
    out={
      'protocol':'SCBA_symmetric_spectral_Kubo_with_four_channel_non_crossing_ladders' if vertex else 'SCBA_symmetric_spectral_Kubo_bare_bubble',
      'symmetric':float(np.real((xy+yx)/2)),
      'xy':float(np.real(xy)), 'yx':float(np.real(yx)),
      'bare_symmetric':float(np.real((bare_xy+bare_yx)/2)),
      'hermiticity_error':float(abs(xy-yx.conjugate())),
      'vertex_residual':float(max(res,default=0.)),
      'scba_residual':float(np.max(abs(sigma-p.disorder_variance*gr.mean((0,1))))),
      'scba_iterations':it,
      'channels':channels,
      'imaginary_symmetric_part':float(np.imag((xy+yx)/2)),
    }
    return out


def slope(p:Lattice, step:float=0.02, protocol:str='onsite',vertex:bool=True):
    if step<=0: raise ValueError('step must be >0')
    plus=spectral_bubble(p,step,protocol,vertex)
    minus=spectral_bubble(p,-step,protocol,vertex)
    return {'value':(plus['symmetric']-minus['symmetric'])/(2*step),
            'plus':plus['symmetric'],'minus':minus['symmetric'],
            'vertex_residual':max(plus['vertex_residual'],minus['vertex_residual']),
            'scba_residual':max(plus['scba_residual'],minus['scba_residual'])}
