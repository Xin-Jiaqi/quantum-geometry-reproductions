"""Window/cutoff checks of the top-valence plus complete d+p empty sector."""
from pathlib import Path
exec(Path('material12.py').read_text().split('for window in ')[0])
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'scope':'Unfitted top-valence plus all d+p empty states, 20 states total; no new DFT','units':'1e-9 (nm V T)^-1','runs':{}}
for window in [4,8,12]:
 grid,_=localized_model(D,window);ix=np.rint(D['k'][:,:2]*mesh).astype(int)%mesh;h=np.array([grid[tuple(i)] for i in ix]);nf=28+window;lifts=[]
 for cf,e in zip(D['coeff'],D['energy']):
  f,_=polar_cols(cf[:,:nf].astype(complex));rest=null_space(adj(f));u,_=polar_cols(((adj(rest)@cf[:,nf:])/(1+((e[nf:]-e[nf])/4)**2)).conj().T);lifts.append(np.vstack([np.zeros((20,44)),adj(f),u@adj(rest)]))
 lift=np.array(lifts);basis=np.zeros((nk,44,20),complex)
 for i,ids in enumerate(layers):
  e,u=np.linalg.eigh(h[:,ids[:,None],ids]);basis[:,ids[:,None],np.arange(10*i,10*i+10)]=u[:,:,12:]
 inter=lift@(h*(BOTTOM[:,None]!=BOTTOM[None,:]))@adj(lift)
 for nb in [96,112]:
  q,metric=orthonormalize((lift@basis)[:,:nb]);h0=D['diag'][:,:nb,:nb];t=inter[:,:nb,:nb];spin=adj(q)@D['spin'][:,:nb,:nb]@q;z=adj(q)@D['Z'][:,:nb,:nb]@q;links=np.array([adj(q)@D['links'][:,a,:nb,:nb]@q[neighbors[:,a]] for a in range(4)]).transpose(1,0,2,3)
  rows={str(alpha):response(adj(q)@(h0+(alpha-1)*t)@q,spin,z,links,4) for alpha in [0.,.5,1.]}
  rows['half_step']=response(adj(q)@h0@q,spin,z,links,4,.001)
  key=f'window{window}_bands{nb}';out['runs'][key]={'metric_minimum':metric,'rows':rows};print(key,{k:v['C'] for k,v in rows.items()},flush=True)
Path('validate20.json').write_text(json.dumps(out,indent=2)+'\n')
