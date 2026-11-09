"""Check existing controls across empty-state cutoffs; no new fitting or scans."""
from pathlib import Path
import os,socket,json,numpy as np
from analyze import Response
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute required')
p=Path.cwd();d=np.load(p/'mono/cache.npz');v=np.load(p/'electrostatic_operators.npz');out={}
for nb in (40,48,56,64):
    r=Response(d,nb);out[str(nb)]={'baseline':r.response(),'controls':{}}
    for name,op in zip(v['names'],v['operators']):
        out[str(nb)]['controls'][str(name)]=r.response(r.h+op[:,:nb,:nb])
(p/'electrostatic_audit.json').write_text(json.dumps(out,indent=2))
print('ELECTROSTATIC_AUDIT',[(nb,v['controls']['lower_3.8_full']['C_total_lower_upper'][0],v['controls']['uniform_F+0.25_3.8']['C_total_lower_upper'][0]) for nb,v in out.items()])
