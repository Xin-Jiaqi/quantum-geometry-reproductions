#!/usr/bin/env python3
"""v39 model-only SCBA/Streda-sea benchmark in normalized units; not full dc conductivity."""
import argparse,json
from qgmap.scba_streda_v39 import derivative_audit,sea_energy_scan,exact_finite_disorder_sea

def reproduce(mesh=21,mass=-1.,W=.12,eta=.12,quadrature=256):
    return dict(schema='qgmap-scba-sea-v39-1',scientific_status='PARTIAL',
     derivative_W=derivative_audit(mesh=mesh,mass=mass,W=W,eta=eta,h=3e-5),
     clean_insulator=sea_energy_scan(mesh=mesh,mass=mass,mu=0,eta=eta,W=W,quadrature=quadrature),
     metallic_control=sea_energy_scan(mesh=mesh,mass=mass,mu=1.1,eta=eta,W=W,quadrature=quadrature),
     exact_finite_sample=exact_finite_disorder_sea(L=4,mass=mass,W=W,eta=eta,seeds=6),
     warning='Bare-current SCBA sea only; RR density vertex does not replace full transverse current ladder or Streda I, nor mixed EB dc')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--mesh',type=int,default=21)
    p.add_argument('--mass',type=float,default=-1.)
    p.add_argument('--W',type=float,default=.12)
    p.add_argument('--eta',type=float,default=.12)
    p.add_argument('--quadrature',type=int,default=256)
    c=p.parse_args()
    print(json.dumps(reproduce(c.mesh,c.mass,c.W,c.eta,c.quadrature),indent=2,allow_nan=False))
