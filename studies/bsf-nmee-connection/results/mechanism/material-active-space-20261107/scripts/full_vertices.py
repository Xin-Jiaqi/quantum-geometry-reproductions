"""Probe-connection controls in the complete material space, without deleting states."""
from pathlib import Path
exec(Path('material12.py').read_text().split('for window in ')[0])
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'scope':'Fixed zero-field active projector; block probe vertices in full 96/112 material space, no state deletion','units':'1e-9 (nm V T)^-1','runs':{}}
for window in [4,8,12]:
 grid,_=localized_model(D,window);ix=np.rint(D['k'][:,:2]*mesh).astype(int)%mesh;h=np.array([grid[tuple(i)] for i in ix]);nf=28+window;lifts=[]
 for cf,e in zip(D['coeff'],D['energy']):
  f,_=polar_cols(cf[:,:nf].astype(complex));rest=null_space(adj(f));u,_=polar_cols(((adj(rest)@cf[:,nf:])/(1+((e[nf:]-e[nf])/4)**2)).conj().T);lifts.append(np.vstack([np.zeros((20,44)),adj(f),u@adj(rest)]))
 lift=np.array(lifts);basis,selection=active_layer_basis(h,layers,weighting='amplitude')
 for nb in [96,112]:
  a,_=orthonormalize((lift@basis)[:,:nb]);p=a@adj(a);q=np.eye(nb)-p
  h0=D['diag'][:,:nb,:nb];s=D['spin'][:,:nb,:nb];z=D['Z'][:,:nb,:nb];ll=D['links'][:,:,:nb,:nb]
  block=lambda x:p@x@p+q@x@q
  hb,sb,zb=[block(x) for x in [h0,s,z]]
  lb=np.array([p@ll[:,j]@p[neighbors[:,j]]+q@ll[:,j]@q[neighbors[:,j]] for j in range(4)]).transpose(1,0,2,3)
  cases={'all':(h0,s,z,ll),'H_block':(hb,s,z,ll),'Z_block':(h0,s,zb,ll),'S_block':(h0,sb,z,ll),'E_block':(h0,s,z,lb),'all_block':(hb,sb,zb,lb),'Z_offdiag_only':(h0,s,z-zb,ll),'S_offdiag_only':(h0,s-sb,z,ll)}
  rows={name:response(hc,sc,zc,lc,48) for name,(hc,sc,zc,lc) in cases.items()}
  key=f'window{window}_bands{nb}';out['runs'][key]={'rows':rows,'Zeeman_linearity_residual':rows['all']['C']-rows['Z_block']['C']-rows['Z_offdiag_only']['C'],'readout_linearity_residual':rows['all']['C']-rows['S_block']['C']-rows['S_offdiag_only']['C']}
  Path('full_vertices.json').write_text(json.dumps(out,indent=2)+'\n');print(key,{k:None if r['C'] is None else round(r['C'],5) for k,r in rows.items()},flush=True)
