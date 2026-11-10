"""Fixed-channel low-order hypothesis test; all scientific work in PBS."""
from pathlib import Path
import os,socket,json,time,hashlib
import numpy as np
from response_core import Response,embedding,adj,polar
from mixing_channels import split_channels
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('Compute node required')
root=Path.cwd();stage=os.getenv('STAGE','train');assert stage in ('train','validate','refine')
alpha={'train':[-.25,0,.25],'validate':[-.5,.5,1.],'refine':[-.5,.5,1.]}[stage]
b=.001 if stage=='refine' else .002;du=.0005 if stage=='refine' else .001
specs={
'MoSSe':{'zero':'/work/sunq/xjq/mosse-strain-20261108/e000/A/cache.npz','nbs':[96,112],
 'minus':'/work/sunq/xjq/mosse-gate-scf-20261109/fm010/A/cache.npz','plus':'/work/sunq/xjq/mosse-gate-scf-20261109/fp010/A/cache.npz',
 'overlaps':'/work/sunq/xjq/mosse-screened-channels-20261109/overlaps.npz','keys':['fm010','fp010']},
'InSe':{'zero':'/work/sunq/xjq/inse-nmee-20261109/A0/A/cache.npz','nbs':[144,160],
 'minus':'/work/sunq/xjq/inse-nmee-20261109/Am/A/cache.npz','plus':'/work/sunq/xjq/inse-nmee-20261109/Ap/A/cache.npz',
 'overlaps':'/work/sunq/xjq/inse-nmee-20261109/overlaps.npz','keys':['Am','Ap']},
'MoS2':{'zero':'/work/sunq/xjq/mos2-interface-20261108/A/cache.npz','nbs':[96,112]}}
res={'complete':False,'stage':stage,'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'order':['net','lower','upper'],'rows':[],'checks':[],'inputs':{},'b_eV':b,'step_halfwidth_eV':du,'alpha':alpha,'new_DFT':False};start=time.monotonic()
def save():
 res['elapsed_seconds']=time.monotonic()-start;(root/f'{stage}.json').write_text(json.dumps(res,indent=2)+'\n')
for mat,s in specs.items():
 data={k:np.load(s[k]) for k in ('zero','minus','plus') if k in s};d0=data['zero'];over=np.load(s['overlaps']) if 'overlaps' in s else None
 res['inputs'][mat]={k:{'path':s[k],'sha256':hashlib.sha256(Path(s[k]).read_bytes()).hexdigest()} for k in ('zero','minus','plus','overlaps') if k in s}
 bottom=d0['orb_bottom'];lb=np.diag(bottom.astype(float));lt=np.eye(len(bottom))-lb;off=bottom[:,None]!=bottom[None,:];layers=[np.flatnonzero(bottom),np.flatnonzero(~bottom)]
 for nb in s['nbs'][-1:] if stage=='refine' else s['nbs']:
  r0=Response(d0,nb)
  for frozen in ([8] if stage=='refine' else [4,8]):
   l0,minimum=embedding(d0,frozen,nb);h0=adj(l0)@r0.h@l0
   _,_,p,diag=split_channels(h0,layers,14);q=np.eye(len(bottom))-p;pb=p@lb;qt=q@lt
   def selected(local):
    z=pb@(local*off)@qt;return z+adj(z)
   w=(np.eye(nb)-2*d0['P'][:,:nb,:nb])/2;assert np.max(abs(w-adj(w)))<1e-8
   bc0=selected(h0);bcw=selected(adj(l0)@w@l0)
   if stage=='refine':
    # Gauge-invariant singular values within exactly degenerate occupied/empty sectors.
    ids_b,ids_t=layers
    eb,ub=np.linalg.eigh(h0[:,ids_b[:,None],ids_b]);et,ut=np.linalg.eigh(h0[:,ids_t[:,None],ids_t])
    z=adj(ub[:,:,:14])@h0[:,ids_b[:,None],ids_t]@ut[:,:,14:]
    delta=et[:,None,14:]-eb[:,:14,None]
    bare_gap=float(delta.min())
    bare_eta=float(np.linalg.svd(z/delta,compute_uv=False).max()) if bare_gap>1e-6 else None
    kfull=l0@bc0@adj(l0);background=r0.h-kfull
    ee,uu=np.linalg.eigh(background);nv=r0.nocc;oc=uu[:,:,:nv];em=uu[:,:,nv:]
    vv=adj(oc)@kfull@em;gap=ee[:,None,nv:]-ee[:,:nv,None]
    eta=np.linalg.svd(vv/gap,compute_uv=False)
    diagnostic={'material':mat,'nb':nb,'frozen_empty':frozen,'bare_cross_gap_eV':bare_gap,'bare_eta_max':bare_eta,
      'background_direct_gap_eV':float(gap.min()),'background_eta_max':float(eta.max()),
      'background_eta_rms':float(np.sqrt(np.mean(abs(vv/gap)**2))),
      'selected_K_operator_norm_eV':float(np.linalg.norm(kfull,ord=2,axis=(-2,-1)).max()),
      'background_gap_perturbation_bound':float(2*np.linalg.norm(kfull,ord=2,axis=(-2,-1)).max()/gap.min())}
    if over is not None:
     means={}
     for cname,key in zip(['minus','plus'],s['keys']):
      liftf,_=polar(over[key][:,:nb,:nb]@l0);rf=Response(data[cname],nb);hf=adj(liftf)@rf.h@liftf
      means[cname]=float((np.trace(hf[:,ids_t[:,None],ids_t],axis1=-2,axis2=-1).real/len(ids_t)-np.trace(hf[:,ids_b[:,None],ids_b],axis1=-2,axis2=-1).real/len(ids_b)).mean())
     diagnostic['screened_upper_minus_lower_trace_change_eV']=means['plus']-means['minus']
     diagnostic['screened_trace_slope_eV_per_Vnm']=(means['plus']-means['minus'])/.2
    res.setdefault('mixing_diagnostics',[]).append(diagnostic)
   # Complete physical vertices remain in Response; only the selected H block and its U dependence scale.
   cases=[('internal_step',name,sgn*du,r0,l0,r0.h+sgn*du*w,bc0+sgn*du*bcw) for name,sgn in [('minus',-1),('zero',0),('plus',1)]]
   if over is not None:
    for name in ('minus','zero','plus'):
     r=Response(data[name],nb)
     if name=='zero':lift=l0;loss=0.
     else:
      ovkey=s['keys'][0 if name=='minus' else 1];lift,sv=polar(over[ovkey][:,:nb,:nb]@l0);loss=float(abs(sv-1).max());assert loss<.02
     cases.append(('screened_external',name,{'minus':-.1,'zero':0.,'plus':.1}[name],r,lift,r.h,selected(adj(lift)@r.h@lift)))
     res['checks'].append({'material':mat,'nb':nb,'frozen_empty':frozen,'case':name,'transport_loss':loss})
   for mode,name,control,r,lift,h,bc in cases:
    bcfull=lift@bc@adj(lift)
    for a in alpha:
     hh=h-(1-a)*bcfull
     row={'material':mat,'nb':nb,'frozen_empty':frozen,'mode':mode,'case':name,'control':control,'alpha':a}
     try:
      ans=r.response(hh,b=b);row.update(ans,valid=True)
      if ans['gap_eV']<.05:row.update(valid=False,reason='gap below .05eV reliability floor')
     except (ValueError,np.linalg.LinAlgError) as exc:row.update(valid=False,reason=str(exc))
     res['rows'].append(row)
   res['checks'].append({'material':mat,'nb':nb,'frozen_empty':frozen,'min_trial_singular':minimum,**diag})
   save();print('DONE',mat,nb,frozen,stage,flush=True)
res['complete']=True;save();print('ALL_COMPLETE',flush=True)
