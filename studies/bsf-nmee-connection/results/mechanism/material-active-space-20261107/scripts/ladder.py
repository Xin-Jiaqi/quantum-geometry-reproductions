"""Diagnose omitted occupied/empty states without fitting response signs."""
from pathlib import Path
exec(Path('material12.py').read_text().split('for window in ')[0])
window=8;nb=96
hgrid,_=localized_model(D,window);ix=np.rint(D['k'][:,:2]*mesh).astype(int)%mesh
h=np.array([hgrid[tuple(i)] for i in ix]);nf=28+window;lifts=[]
for cf,e in zip(D['coeff'],D['energy']):
 f,_=polar_cols(cf[:,:nf].astype(complex));rest=null_space(adj(f));u,_=polar_cols(((adj(rest)@cf[:,nf:])/(1+((e[nf:]-e[nf])/4)**2)).conj().T)
 lifts.append(np.vstack([np.zeros((20,44)),adj(f),u@adj(rest)]))
lift=np.array(lifts);t=h*(BOTTOM[:,None]!=BOTTOM[None,:]);inter=(lift@t@adj(lift))[:,:nb,:nb]
h0=D['diag'][:,:nb,:nb];weak=h0-inter;spin=D['spin'][:,:nb,:nb];z=D['Z'][:,:nb,:nb];links=D['links'][:,:,:nb,:nb]
base,selection=active_layer_basis(h,layers,weighting='amplitude')
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'scope':'No-fit occupied/empty truncation ladder, fixed original material hoppings and projected observables','units':'1e-9 (nm V T)^-1','ladder':{},'groups':{}}

def occupied_groups(hc,b=.002):
 vals=[];boundaries=[]
 for sign in [1,-1]:
  e,u=np.linalg.eigh(hc+sign*b*z);occ=u[:,:,:48];emp=u[:,:,48:];der=np.zeros_like(occ)
  for a in range(4):
   shifted=links[:,a]@occ[neighbors[:,a]]
   der+=(-1)**a*mesh*cell[a//2,0]/(4*np.pi)*shifted@np.linalg.inv(adj(occ)@shifted)
  pairs=2*((adj(occ)@spin@emp)*(adj(emp)@(1j*der)).swapaxes(1,2)).real/(e[:,:48,None]-e[:,None,48:])
  vals.append(pairs.mean(axis=0).sum(axis=1)*6.331927853335898e4)
  boundaries.append([float(np.min(e[:,j]-e[:,j-1])) for j in [20,36,40,44]])
 c=(vals[0]-vals[1])/(2*b)
 return {'groups':{f'{a+1}-{b}':float(c[a:b].sum()) for a,b in [(0,20),(20,36),(36,40),(40,44),(44,48)]},'sum':float(c.sum()),'per_occupied_band':c.tolist(),'minimum_group_boundary_gaps_eV':np.min(boundaries,axis=0).tolist()}

for alpha in [0.,.5,1.]:out['groups'][str(alpha)]=occupied_groups(weak+alpha*inter)
for occupied in [2,4,6,8,10,14]:
 for empty in ['coupled','all']:
  per=occupied+(4 if empty=='coupled' else 8);basis=np.zeros((nk,44,2*per),complex)
  for i,ids in enumerate(layers):
   e,u=np.linalg.eigh(h[:,ids[:,None],ids]);v=u[:,:,14-occupied:14]
   c=np.concatenate([u[:,:,14:16],base[:,ids,6*i+4:6*i+6]],axis=2) if empty=='coupled' else u[:,:,14:]
   basis[:,ids[:,None],np.arange(i*per,(i+1)*per)]=np.concatenate([v,c],axis=2)
  q,metric=orthonormalize((lift@basis)[:,:nb]);s=adj(q)@spin@q;zz=adj(q)@z@q
  ll=np.array([adj(q)@links[:,a]@q[neighbors[:,a]] for a in range(4)]).transpose(1,0,2,3)
  rows={str(alpha):response(adj(q)@(weak+alpha*inter)@q,s,zz,ll,2*occupied) for alpha in [0.,1.]}
  key=f'v{occupied}_{empty}'
  out['ladder'][key]={'dimension':2*per,'metric_minimum':metric,'rows':rows}
  Path('ladder.json').write_text(json.dumps(out,indent=2)+'\n')
  print(key,{a:round(r['C'],5) if r['C'] is not None else None for a,r in rows.items()},flush=True)
q,_=orthonormalize((lift@base)[:,:nb]);p=q@adj(q)
out['spin_readout_controls']={'full':response(h0,spin,z,links,48),'active_PSP':response(h0,p@spin@p,z,links,48),'remainder':response(h0,spin-p@spin@p,z,links,48)}
Path('ladder.json').write_text(json.dumps(out,indent=2)+'\n')
print('GROUPS',out['groups'],flush=True)
