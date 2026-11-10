"""Reuse existing independently screened SCFs, PBS compute nodes only."""
from pathlib import Path
import json,os,socket,time,hashlib
import numpy as np
from split_kernel import adj,eig,match,reorder,factors,split
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute node required')
root=Path.cwd();start=time.monotonic()
specs={
 'MoSSe':{'zero':'/work/sunq/xjq/mosse-strain-20261108/e000/A/cache.npz',
          'minus':'/work/sunq/xjq/mosse-gate-scf-20261109/fm010/A/cache.npz',
          'plus':'/work/sunq/xjq/mosse-gate-scf-20261109/fp010/A/cache.npz',
          'overlaps':'/work/sunq/xjq/mosse-screened-channels-20261109/overlaps.npz',
          'keys':['fm010','fp010'],'cutoffs':[96,112]},
 'InSe':{'zero':'/work/sunq/xjq/inse-nmee-20261109/A0/A/cache.npz',
         'minus':'/work/sunq/xjq/inse-nmee-20261109/Am/A/cache.npz',
         'plus':'/work/sunq/xjq/inse-nmee-20261109/Ap/A/cache.npz',
         'overlaps':'/work/sunq/xjq/inse-nmee-20261109/overlaps.npz',
         'keys':['Am','Ap'],'cutoffs':[128,144,160]}}
result={'complete':False,'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),
        'order':['net','lower','upper'],'units_G':'1e-9 (V^2 T)^-1',
        'gate_interval_V_nm':.2,'convention':'-2 Im(Svc Vcv) / Delta^2, Vcv=-i Delta rcv',
        'scope':'Independent gate SCFs, fixed ions, spin-only frozen-potential B probe. Finite gate secant, not exact zero-field derivative.',
        'rows':[],'inputs':{}}
def save():
 result['elapsed_seconds']=time.monotonic()-start
 (root/('result_refined.json' if os.getenv('REFINE') else 'result.json')).write_text(json.dumps(result,indent=2)+'\n')
for material,spec in specs.items():
 data={s:np.load(spec[s]) for s in ('zero','minus','plus')};over=np.load(spec['overlaps'])
 result['inputs'][material]={s:{'path':spec[s],'sha256':hashlib.sha256(Path(spec[s]).read_bytes()).hexdigest()} for s in ('zero','minus','plus','overlaps')}
 nv=int(data['zero']['nocc']);nk=len(data['zero']['k'])
 for s in ('minus','plus'):
  assert int(data[s]['nocc'])==nv and np.max(abs(data[s]['cell']-data['zero']['cell']))<1e-9
  assert np.max(abs(data[s]['k']-data['zero']['k']))<1e-9
 for nb in spec['cutoffs'][-1:] if os.getenv('REFINE') else spec['cutoffs']:
  om=over[spec['keys'][0]][:,:nb,:nb];op=over[spec['keys'][1]][:,:nb,:nb]
  sv=np.concatenate([np.linalg.svd(om,compute_uv=False).ravel(),np.linalg.svd(op,compute_uv=False).ravel()])
  for b in ((.0005,.00025,.000125,.0000625) if os.getenv('REFINE') else (.002,.001)):
   accum={name:{q:np.zeros(3) for q in ('energy','matrix','total','energy_single_den','matrix_single_den','Cplus','Cminus')} for name in ('anchor','direct','energy_order')}
   details={name:[] for name in accum};maps={name:np.zeros((3,nk)) for name in accum}
   for sign in (1,-1):
    es={s:eig(d,nb,sign*b) for s,d in data.items()}
    e0,u0=es['zero'];em,um=es['minus'];ep,up=es['plus']
    pm,wm=match(u0,um,om,nv);pp,wp=match(u0,up,op,nv)
    pd,wd=match(um,up,op@adj(om),nv)
    identity=np.tile(np.arange(nb),(nk,1))
    for name,(im,ip) in {'anchor':(pm,pp),'direct':(identity,pd),'energy_order':(identity,identity)}.items():
     enm,unm=reorder(em,um,im);enp,unp=reorder(ep,up,ip)
     nm,dm,am=factors(data['minus'],nb,enm,unm);np_,dp,ap=factors(data['plus'],nb,enp,unp)
     ee,mm,tt=split(np_,dp,nm,dm,.2)
     deltap=enp[:,:nv,None]-enp[:,None,nv:]
     deltam=enm[:,:nv,None]-enm[:,None,nv:]
     e1,m1,t1=split(np_/deltap,1/deltap,nm/deltam,1/deltam,.2)
     assert np.max(abs(t1-tt))<1e-5, 'Kernel-convention identity failed'
     closure=float(np.max(abs(ee+mm-tt)))
     for q,values in [('energy',ee),('matrix',mm),('total',tt),('energy_single_den',e1),('matrix_single_den',m1)]:
      ck=values.sum(axis=(2,3));accum[name][q]+=sign*ck.mean(axis=1)/(2*b)
      if q=='energy':maps[name]+=sign*ck/(2*b)
     for q,n,d in [('Cplus',np_,dp),('Cminus',nm,dm)]:accum[name][q]+=sign*(n*d).sum(axis=(2,3)).mean(axis=1)/(2*b)
     details[name].append(dict(B_sign=sign,product_closure_max=closure,minus=am,plus=ap))
   for name,values in accum.items():
    row={'material':material,'nb':nb,'b_eV':b,'matching':name,**{q:v.tolist() for q,v in values.items()},'details':details[name],
         'min_full_subspace_singular':float(sv.min()),'anchor_min_pair_power':float(min(wm.min(),wp.min())),
         'anchor_pair_power_p01':float(np.quantile(np.r_[wm.ravel(),wp.ravel()],.01)),
         'direct_min_pair_power':float(wd.min()),'energy_k_net_lower_upper':maps[name].T.tolist()}
    result['rows'].append(row)
   save();print('COMPLETE',material,nb,b,'anchor',accum['anchor']['energy'],accum['anchor']['matrix'],accum['anchor']['total'],flush=True)
 result.setdefault('materials_complete',[]).append(material);save()
result['complete']=True;save();print('ALL_COMPLETE',flush=True)
