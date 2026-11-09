"""Joint zero-gate and screened-gate channel intervention, full KS response."""
from pathlib import Path
import os,socket,json
import numpy as np
from analyze import Response,embedding,adj,polar
from mixing_channels import split_channels
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute node required')
root=Path.cwd();source=Path('/work/sunq/xjq/mosse-strain-20261108/e000/A')
gates=Path('/work/sunq/xjq/mosse-gate-scf-20261109')
d0=np.load(source/'cache.npz');ov=np.load(root/'overlaps.npz')
data={'zero':d0,**{case:np.load(gates/case/'A/cache.npz') for case in ('fm010','fp010')}}
fields={'zero':0.,'fm010':-.1,'fp010':.1}
bottom=d0['orb_bottom'];layers=[np.flatnonzero(bottom),np.flatnonzero(~bottom)]
lb=np.diag(bottom.astype(float));lt=np.eye(len(bottom))-lb
off=bottom[:,None]!=bottom[None,:]
result={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'rows':[],'checks':[],
    'units':'C:1e-9 (nm V T)^-1; G:1e-9 (V^2 T)^-1; field: applied V/nm',
    'scope':'Frozen original gate-SCF potentials; fixed zero-gate channel projectors transported by PAW overlap. No new SCF.',
    'order':'total, lower, upper','complete':False}
def save():(root/'result.json').write_text(json.dumps(result,indent=2))
for nb in (96,112):
    r0=Response(d0,nb)
    for frozen in (4,8,12):
        l0,minimum=embedding(d0,frozen,nb)
        _,_,p,checks=split_channels(adj(l0)@r0.h@l0,layers,14)
        q=np.eye(len(bottom))-p;pb=p@lb;pt=p@lt;qb=q@lb;qt=q@lt
        for case,d in data.items():
            r=Response(d,nb)
            if case=='zero':lift=l0;loss=0.;smin=1.
            else:
                projected=ov[case][:,:nb,:nb]@l0
                lift,s=polar(projected);smin=float(s.min());loss=float(np.max(abs(s-1)))
                if loss>.02:raise ValueError(f'Poor common-subspace representation {case} {nb} {frozen}: {loss}')
            hl=adj(lift)@r.h@lift;t=hl*off
            bc=pb@t@qt;bc=bc+adj(bc);tc=pt@t@qb;tc=tc+adj(tc)
            pieces={'VC':bc+tc,'Vb_Ct':bc,'Vt_Cb':tc,'VV_CC':p@t@p+q@t@q}
            recon=float(abs(pieces['VC']+pieces['VV_CC']-t).max())
            assert recon<1e-8
            controls=[('full',1.,np.zeros_like(t))]
            controls += [('VC',a,(1-a)*pieces['VC']) for a in (0.,.5,.75)]
            controls += [(name,0.,pieces[name]) for name in ('VV_CC','Vb_Ct','Vt_Cb')]
            for name,alpha,removed in controls:
                h=r.h-lift@removed@adj(lift);ans=r.response(h)
                if ans['gap_eV']<.05:raise ValueError('Gap below .05 eV')
                row={'nb':nb,'frozen_empty':frozen,'case':case,'F':fields[case],
                    'channel':name,'alpha':alpha,**ans}
                if nb==112 and frozen==8:
                    half=r.response(h,b=.001)
                    row['half_probe_max_change']=float(abs(np.array(ans['C_total_lower_upper'])-half['C_total_lower_upper']).max())
                result['rows'].append(row)
            result['checks'].append({'nb':nb,'frozen_empty':frozen,'case':case,
                'representation_max_singular_deviation':loss,'min_transport_singular':smin,
                'min_trial_singular':minimum,'channel_reconstruction_eV':recon,**checks})
            save();print('SCREENED_CHANNELS',nb,frozen,case,flush=True)
result['complete']=True;save();print('SCREENED_CHANNELS_COMPLETE',flush=True)
