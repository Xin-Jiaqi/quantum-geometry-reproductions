"""Identify which empty states and probe vertices are lost by direct truncation."""
from pathlib import Path
exec(Path('ladder.py').read_text().split('for alpha in [0.,.5,1.]')[0])
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'scope':'Predefined energy-pair empty-sector tests and zero-field/probe block controls; no fitted model parameters','units':'1e-9 (nm V T)^-1','pairs':{}}
eig=[np.linalg.eigh(h[:,ids[:,None],ids]) for ids in layers]
for selected in [(0,),(1,),(2,),(0,1),(0,2),(1,2)]:
 per=4+2*len(selected);basis=np.zeros((nk,44,2*per),complex)
 for i,ids in enumerate(layers):
  e,u=eig[i];columns=list(range(12,16))+[b for pair in selected for b in range(16+2*pair,18+2*pair)]
  basis[:,ids[:,None],np.arange(per*i,per*(i+1))]=u[:,:,columns]
 q,metric=orthonormalize((lift@basis)[:,:nb]);s=adj(q)@spin@q;zz=adj(q)@z@q;ll=np.array([adj(q)@links[:,a]@q[neighbors[:,a]] for a in range(4)]).transpose(1,0,2,3)
 rows={str(alpha):response(adj(q)@(weak+alpha*inter)@q,s,zz,ll,4) for alpha in [0.,1.]}
 out['pairs'][str(selected)]={'dimension':2*per,'rows':rows};print('PAIR',selected,{k:v['C'] for k,v in rows.items()},flush=True)
# Complete the active 12-space with its high-empty orthogonal complement.
qd=np.zeros((nk,44,8),complex)
for i,ids in enumerate(layers):
 high=eig[i][1][:,:,16:];active=base[:,ids,6*i+4:6*i+6];coord=adj(high)@active
 dark=np.array([null_space(adj(c)) for c in coord])
 qd[:,ids[:,None],np.arange(4*i,4*i+4)]=high@dark
qa,_=orthonormalize((lift@base)[:,:nb]);raw=(lift@qd)[:,:nb];raw-=qa@(adj(qa)@raw);qdKS,_=orthonormalize(raw);q=np.concatenate([qa,qdKS],axis=2)
ss=adj(q)@spin@q;zz=adj(q)@z@q;hh=adj(q)@h0@q;ll=np.array([adj(q)@links[:,a]@q[neighbors[:,a]] for a in range(4)]).transpose(1,0,2,3)
mask=np.zeros((20,20));mask[:12,:12]=1;mask[12:,12:]=1
variants={'all':(hh,ss,zz,ll),'H_block':(hh*mask,ss,zz,ll),'Z_block':(hh,ss,zz*mask,ll),'S_block':(hh,ss*mask,zz,ll),'E_block':(hh,ss,zz,ll*mask),'H_Z_block':(hh*mask,ss,zz*mask,ll),'H_Z_E_block':(hh*mask,ss,zz*mask,ll*mask),'all_block':(hh*mask,ss*mask,zz*mask,ll*mask)}
out['vertices']={name:response(hc,sc,zc,lc,4) for name,(hc,sc,zc,lc) in variants.items()}
out['checks']={'isometry':float(np.max(abs(adj(q)@q-np.eye(20)))),'H_active_dark_eV':norm(hh[:,:12,12:]),'Zeeman_active_dark':norm(zz[:,:12,12:]),'spin_active_dark':norm(ss[:,:12,12:])}
Path('closure.json').write_text(json.dumps(out,indent=2)+'\n');print('VERTICES',out['vertices'],flush=True)
