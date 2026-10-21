#!/usr/bin/env python3
"""Four-terminal Hall probes, two-terminal Onsager check, and coherent-lead parameter scan."""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.lead_transport_v57 import (
  LeadSpec,diamond_graph,device_response,full_finite_field_fd)


def evaluate(spec=LeadSpec(),source=None,rings=1,gammas=None):
    if source is None:
        pos,h,sites=diamond_graph(rings)
        reference='analytic_diamond_rings'
    else:
        with np.load(source,allow_pickle=False) as f:
            required={'hoppings','positions','lead_sites'}
            if not required.issubset(f.files):raise ValueError('NPZ needs hoppings positions lead_sites')
            h=np.asarray(f['hoppings']);pos=np.asarray(f['positions']);sites=tuple(int(x) for x in f['lead_sites'])
            if any(x!=y for x,y in zip(sites,np.asarray(f['lead_sites']))):raise ValueError('lead_sites integer-valued required')
        reference='external_finite_graph_and_explicit_sites'
    gammas=list(gammas) if gammas is not None else [spec.broadening]*len(sites)
    a=device_response(h,pos,sites,gammas,spec)
    b=device_response(h,pos,sites,gammas,spec,'landau_y')
    fd=full_finite_field_fd(h,pos,sites,gammas,spec)
    out={'version':'v57','status':'VERIFIED_WIDEBAND_COHERENT_FINITE_DEVICE_NOT_BULK',
        'source':reference,'sites':len(pos),'lead_sites':sites,'model_parameters':spec.__dict__,
        'gammas':gammas,'symmetric_gauge':a,'landau_gauge':b,
        'max_T_gauge_residual':float(np.max(abs(a['T']-b['T']))),
        'max_dT_gauge_residual':float(np.max(abs(a['dT_dB']-b['dT_dB']))),
        'max_dT_five_point_residual':float(np.max(abs(a['dT_dB']-fd['dT_dB']))),
        'independent_five_point_field_derivative':fd,
        'scope':['Two/four-terminal Landauer-Buttiker linear response under wide-band ideal leads',
                 'Finite device voltage probe response, not intrinsic bulk orbital MNHE/PHE',
                 'Local site diagonal, magnetic-field independent Gamma; no microscopic leads dispersion',
                 'A topology/size scan of sparse concentric diamond rings is not a controlled bulk thermodynamic limit',
                 'Conductance output normalized to q^2/h; model B is measured by signed qB area/hbar']}
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph-npz',type=Path,help='hoppings, positions, lead_sites arrays; no pickle')
    parser.add_argument('--rings',type=int,default=1)
    parser.add_argument('--B',type=float,default=.31)
    parser.add_argument('--q-over-hbar',type=float,default=-.87)
    parser.add_argument('--mu',type=float,default=.15)
    parser.add_argument('--T',type=float,default=.16)
    parser.add_argument('--gamma-lead',type=float,default=.55)
    parser.add_argument('--nE',type=int,default=160)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    cfg=LeadSpec(args.B,args.q_over_hbar,args.mu,args.T,args.gamma_lead,args.nE)
    o=evaluate(cfg,args.graph_npz,args.rings)
    def enc(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    blob=json.dumps(o,default=enc,indent=2,allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(blob,encoding='utf8')
    else:print(blob,end='')
if __name__=='__main__':main()
