"""Linear gate-sensitivity diagnostic: interlayer step versus residual field.
V_linear = V_step + V_residual, with layer-center separation fixed to Mo planes.
The split is explicitly defined, not unique. Additivity is only the tangent
derivative in the control potential, never an additive decomposition of NMEE.
"""
from pathlib import Path
import json,os,socket
import numpy as np
from analyze import Response,embedding,adj
from mixing_channels import split_channels
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute node required')
root=Path.cwd();source=Path('/work/sunq/xjq/mosse-strain-20261108/e000')
d=np.load(source/'A/cache.npz');ops=np.load(root/'operators.npz')
bottom=d['orb_bottom'];layers=[np.flatnonzero(bottom),np.flatnonzero(~bottom)]
off=bottom[:,None]!=bottom[None,:]
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'scope':__doc__,
     'step':'d_Mo/20 * (I-2 P_lower), eV per model V/nm','rows':[]}
for nb in (96,112):
    r=Response(d,nb);v=ops['A_gate'][:,:nb,:nb]
    step=float(ops['Mo_separation_A'])/20*(np.eye(nb)-2*d['P'][:,:nb,:nb])
    fields={'linear':v,'layer_step':step,'residual':v-step}
    configurations=[('full',None,r.h,fields)]
    for frozen in (4,8,12):
        lift,_=embedding(d,frozen,nb);local=adj(lift)@r.h@lift
        _,parts,p,_=split_channels(local,layers,14);q=np.eye(len(bottom))-p
        def vc(matrix):
            t=(adj(lift)@matrix@lift)*off
            return t-p@t@p-q@t@q
        h=r.h-lift@parts['VC']@adj(lift)
        controls={name:vv-lift@vc(vv)@adj(lift) for name,vv in fields.items()}
        configurations.append(('no_VC',frozen,h,controls))
    for name,frozen,h,controls in configurations:
        for delta in (.01,.02):
            slopes={}
            for channel,pert in controls.items():
                plus=r.response(h+delta*pert);minus=r.response(h-delta*pert)
                slopes[channel]=((np.array(plus['C_total_lower_upper'])-minus['C_total_lower_upper'])/(2*delta)).tolist()
            error=np.array(slopes['linear'])-slopes['layer_step']-slopes['residual']
            out['rows'].append({'nb':nb,'control':name,'frozen_empty':frozen,'delta':delta,
                                'G':slopes,'tangent_additivity_residual':error.tolist()})
            (root/'gate_direction.json').write_text(json.dumps(out,indent=2))
            print(nb,name,frozen,delta,slopes,flush=True)
out['complete']=True;(root/'gate_direction.json').write_text(json.dumps(out,indent=2))
