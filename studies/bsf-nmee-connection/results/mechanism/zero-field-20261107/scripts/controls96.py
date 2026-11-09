"""112-state frozen-potential diagnostic. Reject a claimed zero-field mechanism if baseline is magnetic."""
import os,socket,json
from pathlib import Path
import numpy as np
if not os.environ.get('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute required')
import sys
sys.path.insert(0,str(Path.cwd().parent))
from model import localized_model,polar_cols,BOTTOM
from scipy.linalg import null_space
D=np.load('full-controls-cache.npz'); nk=len(D['k']); n=round(nk**.5); links=D['links']; neighbors=D['neighbors'];cell=D['cell'];Z=D['Z'];spin=D['spin']; P=D['P']; diag=D['diag']
links=links[:,:,:96,:96];Z=Z[:,:96,:96];spin=spin[:,:96,:96];P=P[:,:96,:96];diag=diag[:,:96,:96]
adj=lambda a:a.conj().swapaxes(-1,-2)
def chi(h,spin_op=spin):
 e,u=np.linalg.eigh(h);occ=u[:,:,:48];emp=u[:,:,48:];deriv=np.zeros_like(occ,dtype=complex)
 for a in range(4):
  shifted=links[:,a]@occ[neighbors[:,a]]
  deriv+=(-1)**a*n*cell[a//2,0]/(4*np.pi)*shifted@np.linalg.inv(adj(occ)@shifted)
 pairs=2*((adj(occ)@spin_op@emp)*(adj(emp)@(1j*deriv)).swapaxes(1,2)).real/(e[:,:48,None]-e[:,None,48:])
 return float(pairs.sum(axis=(1,2)).mean()*6.331927853335898e4),float(e[:,48].min()-e[:,47].max())
def response(h,step=.002,spin_op=spin):
 a,ga=chi(h+step*Z,spin_op);b,gb=chi(h-step*Z,spin_op)
 return {'C':(a-b)/(2*step),'gap':min(ga,gb),'even_chi':(a+b)/2,'step':step}
H0=diag
out={'scope':'96-state projection of genuinely zero-B self-consistent KS reference, frozen-potential bare Zeeman derivative; no gate self-consistency','node':socket.gethostname(),'job':os.environ['PBS_JOBID'],'mesh':n}
zero,gap=chi(H0);out['baseline']={'chi_zero':zero,'gap':gap,'steps':[response(H0,b) for b in (.001,.002,.008,.016)],'zero_reference':chi(diag)[0]}
out['zero_chi_ratio']=abs(zero)/(abs(out['baseline']['steps'][1]['C'])*.002+1e-15)
out['total_baseline']=response(H0,spin_op=D['spin_total'][:,:96,:96])
out['total_chi_zero']=chi(H0,D['spin_total'][:,:96,:96])[0]
out['windows']={}
for frozen in (4,8,12):
 grid,_=localized_model(D,frozen); ix=np.rint(D['k'][:,:2]*n).astype(int)%n
 h=np.array([grid[tuple(i)] for i in ix]);lifts=[];nf=28+frozen
 for cf,e in zip(D['coeff'],D['energy']):
  f,_=polar_cols(cf[:,:nf].astype(complex));rest=null_space(adj(f));u,_=polar_cols(((adj(rest)@cf[:,nf:])/(1+((e[nf:]-e[nf])/4)**2)).conj().T)
  lifts.append(np.vstack([np.zeros((20,44)),adj(f),u@adj(rest)]))
 lift=np.array(lifts);HT=lift@(h*(BOTTOM[:,None]!=BOTTOM[None,:]))@adj(lift)
 HT=HT[:,:96,:96]
 weak=H0-HT;strong=H0;sources=[np.linalg.eigh(x) for x in (weak,strong)]
 rows={}
 for p in (0,1):
  for q in (0,1):
   ep,up=sources[p];eq=sources[q][0];de=ep[:,:,None]-ep[:,None,:];dq=eq[:,:,None]-eq[:,None,:]
   if np.any((abs(de)<1e-9)&(abs(dq)>1e-7)):raise ValueError('degenerate spectral ambiguity')
   hc=(up*eq[:,None,:])@adj(up);rows[f'states{p}_energy{q}']=response(hc)
 gate=[dict(U=float(U),total=response(H0+U*(.5*np.eye(96)-P),spin_op=D['spin_total'][:,:96,:96]),**response(H0+U*(.5*np.eye(96)-P))) for U in (-.6,-.45,-.3,-.15,0,.15,.3,.45,.6)] if frozen==8 else []
 out['windows'][str(frozen)]={'spectral':rows,'gate':gate,'hopping_scan':[dict(alpha=t,**response(H0+(t-1)*HT)) for t in (0,.25,.5,.75,1)]}
 Path('controls96-analysis.json').write_text(json.dumps(out,indent=2)+'\n')
Path('controls96-analysis.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
