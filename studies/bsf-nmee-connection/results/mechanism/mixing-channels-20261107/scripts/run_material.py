"""Occupied/empty channel controls using the existing zero-B material cache."""
import os,socket,json,sys
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
if not os.environ.get('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute required')
sys.path.insert(0,str(Path.cwd().parent))
from model import localized_model,polar_cols,BOTTOM
from mixing_channels import split_channels,adj
D=np.load('../zero-field-refined-20261107/full-controls-cache.npz')
nk=len(D['k']);n=round(nk**.5);cell=D['cell'];neighbors=D['neighbors']
combinations={'none':(), 'VV':('VV',),'CC':('CC',),'VC':('VC',),'VV_CC':('VV','CC'),'VV_VC':('VV','VC'),'CC_VC':('CC','VC'),'all':('VV','CC','VC')}
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'scope':'zero-field fixed layer valence/empty projectors; 14 of 22 states occupied in each d+p layer; artificial hopping interventions, no strain/gate DFT','units':'1e-9 (nm V T)^-1','runs':{}}
for frozen in (4,8,12):
 grid,_=localized_model(D,frozen);ix=np.rint(D['k'][:,:2]*n).astype(int)%n
 h=np.array([grid[tuple(i)] for i in ix]);lifts=[];nf=28+frozen
 for cf,e in zip(D['coeff'],D['energy']):
  f,_=polar_cols(cf[:,:nf].astype(complex));rest=null_space(adj(f));u,_=polar_cols(((adj(rest)@cf[:,nf:])/(1+((e[nf:]-e[nf])/4)**2)).conj().T)
  lifts.append(np.vstack([np.zeros((20,44)),adj(f),u@adj(rest)]))
 lift=np.array(lifts)
 _,parts,p,checks=split_channels(h,[np.where(BOTTOM)[0],np.where(~BOTTOM)[0]],14)
 embedded={name:lift@mat@adj(lift) for name,mat in parts.items()}
 checks['lift_isometry']=float(np.max(abs(adj(lift)@lift-np.eye(44))))
 for nb in (96,112):
  Z=D['Z'][:,:nb,:nb];spin=D['spin'][:,:nb,:nb];links=D['links'][:,:,:nb,:nb]
  h0=D['diag'][:,:nb,:nb];blocks={key:v[:,:nb,:nb] for key,v in embedded.items()};weak=h0-sum(blocks.values())
  ew,uw=np.linalg.eigh(weak);pw=uw[:,:,:48]@adj(uw[:,:,:48])
  def chi(h):
   e,u=np.linalg.eigh(h);occ=u[:,:,:48];emp=u[:,:,48:];v=np.zeros_like(occ,dtype=complex)
   for a in range(4):
    shifted=links[:,a]@occ[neighbors[:,a]]
    v+=(-1)**a*n*cell[a//2,0]/(4*np.pi)*shifted@np.linalg.inv(adj(occ)@shifted)
   pairs=2*((adj(occ)@spin@emp)*(adj(emp)@(1j*v)).swapaxes(1,2)).real/(e[:,:48,None]-e[:,None,48:])
   return float(pairs.sum(axis=(1,2)).mean()*6.331927853335898e4),float(e[:,48].min()-e[:,47].max())
  def response(h,b):
   plus,gp=chi(h+b*Z);minus,gm=chi(h-b*Z)
   return {'C':(plus-minus)/(2*b),'gap':min(gp,gm),'even_chi':(plus+minus)/2,'step_eV':b}
  rows={}
  for label,include in combinations.items():
   hc=weak+sum((blocks[t] for t in include),np.zeros_like(weak))
   e,u=np.linalg.eigh(hc);pc=u[:,:,:48]@adj(u[:,:,:48])
   rows[label]={'response':response(hc,.002),'half_step':response(hc,.001),'occupied_projector_change':float(np.sqrt(np.mean(np.sum(abs(pc-pw)**2,axis=(1,2))))),'hermiticity':float(np.max(abs(hc-adj(hc))))}
  key=f'window{frozen}_bands{nb}'
  out['runs'][key]={'checks':checks,'rows':rows,'channel_norms':{key:float(np.sqrt(np.mean(np.sum(abs(mat)**2,axis=(1,2))))) for key,mat in blocks.items()}}
  Path('material-channels.json').write_text(json.dumps(out,indent=2)+'\n')
  print(key,{k:round(v['response']['C'],5) for k,v in rows.items()},flush=True)
