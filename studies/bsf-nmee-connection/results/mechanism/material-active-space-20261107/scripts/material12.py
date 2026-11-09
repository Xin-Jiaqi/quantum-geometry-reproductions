"""Response-independent active-sector tests using a true zero-field material cache."""
import os,socket,json,sys
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
if not os.environ.get('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute required')
sys.path.insert(0,str(Path.cwd().parent))
from model import localized_model,polar_cols,BOTTOM,NAMES
from material_subspaces import active_layer_basis,orthonormalize,adj
D=np.load('../zero-field-refined-20261107/full-controls-cache.npz')
nk=len(D['k']);mesh=round(nk**.5);cell=D['cell'];neighbors=D['neighbors']
layers=[np.where(BOTTOM)[0],np.where(~BOTTOM)[0]]
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'scope':'No response fit; valence/low-empty/coupled-high-empty subspaces from zero-field Hamiltonian only; no new DFT','units':'1e-9 (nm V T)^-1','runs':{}}

def chi(h,spin,links,nocc):
 e,u=np.linalg.eigh(h);occ=u[:,:,:nocc];emp=u[:,:,nocc:];der=np.zeros_like(occ)
 min_overlap=1.
 for a in range(4):
  shifted=links[:,a]@occ[neighbors[:,a]];metric=adj(occ)@shifted
  min_overlap=min(min_overlap,float(np.linalg.svd(metric,compute_uv=False).min()))
  der+=(-1)**a*mesh*cell[a//2,0]/(4*np.pi)*shifted@np.linalg.inv(metric)
 pairs=2*((adj(occ)@spin@emp)*(adj(emp)@(1j*der)).swapaxes(1,2)).real/(e[:,:nocc,None]-e[:,None,nocc:])
 return float(pairs.sum(axis=(1,2)).mean()*6.331927853335898e4),float(e[:,nocc].min()-e[:,nocc-1].max()),min_overlap

def response(h,spin,z,links,nocc,b=.002):
 a,ga,oa=chi(h+b*z,spin,links,nocc);bval,gb,ob=chi(h-b*z,spin,links,nocc)
 return {'C':(a-bval)/(2*b) if min(ga,gb)>0 else None,'gap':min(ga,gb),'even_chi':(a+bval)/2,'minimum_occupied_overlap_singular':min(oa,ob)}

def norm(a):return float(np.sqrt(np.mean(np.sum(abs(a)**2,axis=(-2,-1)))))

for window in [int(x) for x in os.environ.get('WINDOWS','4,8,12').split(',')]:
 grid,_=localized_model(D,window);ix=np.rint(D['k'][:,:2]*mesh).astype(int)%mesh
 h=np.array([grid[tuple(i)] for i in ix]);nf=28+window;lifts=[]
 for cf,e in zip(D['coeff'],D['energy']):
  f,_=polar_cols(cf[:,:nf].astype(complex));rest=null_space(adj(f));u,_=polar_cols(((adj(rest)@cf[:,nf:])/(1+((e[nf:]-e[nf])/4)**2)).conj().T)
  lifts.append(np.vstack([np.zeros((20,44)),adj(f),u@adj(rest)]))
 lift=np.array(lifts);inter=BOTTOM[:,None]!=BOTTOM[None,:];t=h*inter
 for method in ['amplitude','mixing']:
  basis,selection=active_layer_basis(h,layers,weighting=method)
  ta=adj(basis)@t@basis;ha=adj(basis)@h@basis
  vrmask=np.zeros((12,12));ccmask=np.zeros((12,12))
  for off,other in [(0,6),(6,0)]:
   vrmask[off:off+2,other+4:other+6]=1;vrmask[other+4:other+6,off:off+2]=1
   ccmask[off+2:off+6,other+4:other+6]=1;ccmask[other+4:other+6,off+2:off+6]=1
  vr=ta*vrmask;cc=ta*ccmask
  norms={};weights={};orb=np.repeat(NAMES,2);element=np.repeat(['Mo','Mo','S','S','Se','Se'],[10,10,6,6,6,6])
  for i,name in enumerate(['bottom','top']):
   off=6*i;other=6*(1-i)
   cv=ta[:,off+4:off+6,other:other+2];ce=ta[:,off+4:off+6,other+2:other+4];cr=ta[:,off+4:off+6,other+4:other+6]
   norms[name]={'r_from_v_eV':norm(cv),'r_from_c_eV':norm(ce),'r_from_r_eV':norm(cr),'c_over_v_amplitude':norm(ce)/norm(cv),'rr_over_v_amplitude':norm(cr)/norm(cv)}
   for label,start in [('v',0),('c',2),('r',4)]:
    b=basis[:,:,off+start:off+start+2]
    weights[name+'_'+label]={el+'_'+orbital:float(np.mean(np.sum(abs(b[:,(element==el)&(orb==orbital)])**2,axis=(1,2)))/2) for el in ['Mo','S','Se'] for orbital in sorted(set(orb)) if np.any((element==el)&(orb==orbital))}
  embedded={'full':lift@t@adj(lift),'active':lift@(basis@ta@adj(basis))@adj(lift),'vr':lift@(basis@vr@adj(basis))@adj(lift),'remote_cc':lift@(basis@cc@adj(basis))@adj(lift)}
  for nb in [int(x) for x in os.environ.get('BANDS','96').split(',')]:
   z=D['Z'][:,:nb,:nb];spin=D['spin'][:,:nb,:nb];links=D['links'][:,:,:nb,:nb];h0=D['diag'][:,:nb,:nb]
   blocks={key:value[:,:nb,:nb] for key,value in embedded.items()};weak=h0-blocks['full']
   q,metric=orthonormalize((lift@basis)[:,:nb]);z12=adj(q)@z@q;s12=adj(q)@spin@q
   links12=np.array([adj(q)@links[:,a]@q[neighbors[:,a]] for a in range(4)]).transpose(1,0,2,3)
   rows={}
   for alpha in [0.,.5,1.]:
    hf=weak+alpha*blocks['full'];hs=weak+alpha*blocks['active']
    rows['full_'+str(alpha)]=response(hf,spin,z,links,48)
    rows['active_in_background_'+str(alpha)]=response(hs,spin,z,links,48)
    rows['bare12_'+str(alpha)]=response(adj(q)@hf@q,s12,z12,links12,4)
   for name,hc in [('full_remove_vr',h0-blocks['vr']),('full_remove_remote_cc',h0-blocks['remote_cc']),('active_remove_vr',weak+blocks['active']-blocks['vr'])]:
    rows[name]=response(hc,spin,z,links,48)
   for name,hc in [('bare12_remove_vr',h0-blocks['vr']),('bare12_remove_remote_cc',h0-blocks['remote_cc'])]:
    rows[name]=response(adj(q)@hc@q,s12,z12,links12,4)
   rows['bare12_half_step']=response(adj(q)@h0@q,s12,z12,links12,4,.001)
   bg=rows['full_0.0']['C']-rows['bare12_0.0']['C']
   predictions={str(a):rows['bare12_'+str(a)]['C']+bg for a in [0.,.5,1.]}
   key=f'window{window}_{method}_bands{nb}'
   out['runs'][key]={'selection':selection,'cutoff_metric_minimum':metric,'coupling_norms':norms,'orbital_weights':weights,'rows':rows,'constant_background_approximation':predictions,'active_hamiltonian_projection_error_eV':float(np.max(abs(adj(lift@basis)@D['diag']@(lift@basis)-ha)))}
   Path('material12.json').write_text(json.dumps(out,indent=2)+'\n')
   print(key,{k:None if v['C'] is None else round(v['C'],5) for k,v in rows.items()},flush=True)
