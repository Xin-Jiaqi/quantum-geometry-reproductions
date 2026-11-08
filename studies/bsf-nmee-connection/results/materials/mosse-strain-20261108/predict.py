"""Before-DFT strain hypothesis: geometric scaling of existing cross-layer blocks.
No fitted deformation potentials, no claim to predict relaxation or experimental thresholds.
"""
from pathlib import Path
import os,socket,json,hashlib,shutil
import numpy as np
from scipy.linalg import null_space
from ase.io import read
from analyze import Response,polar,adj
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute only')
root=Path.cwd();old=Path('/work/sunq/xjq/mosse-7269-mechanism-20261104/zero-field-refined-20261107')
D=np.load(old/'full-controls-cache.npz');at=read(old/'seed/POSCAR');mesh=round(len(D['k'])**.5)
# Keep original KS vertices and full material Hilbert space; change selected interlayer H only.
energy=np.real(np.diagonal(D['diag'],axis1=-2,axis2=-1))
d={k:D[k] for k in ('Z','spin','spin_total','P','links','neighbors','k','cell')};d.update(energies=energy,nocc=48,ncore=20)
xyz=np.repeat(at.positions,[10,10,6,6,6,6],axis=0)
frac=np.repeat(at.get_scaled_positions(),[10,10,6,6,6,6],axis=0)
bottom=xyz[:,2]<np.mean([at.positions[:,2].min(),at.positions[:,2].max()]);cross=bottom[:,None]!=bottom[None,:]
index=np.rint(D['k'][:,:2]*mesh).astype(int)%mesh
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'before_new_strained_DFT':True,
 'assumptions':'Fixed vertical positions/onsite/vertices. Cross-layer real-space hopping scales as (r0/r)^p. p=3,5 are uncalibrated sensitivity assumptions, not fitted material constants. In-plane metric/area changed consistently. Strain relaxation and intralayer deformation potentials omitted.',
 'units':'1e-9 (nm V T)^-1','source':str(old),'rows':[],'source_keys':list(D.keys()),'available_bytes':shutil.disk_usage(root).free}
for nb in (96,112):
 cf=D['coeff'][:,:,:nb-20];lifts=[]
 for coeff,en in zip(cf,energy[:,:nb]):
  nf=36;f,_=polar(coeff[:,:nf]);rest=null_space(adj(f));u,_=polar((adj(rest)@coeff[:,nf:]/(1+((en[20+nf:]-en[20+nf])/4)**2)).conj().T)
  lifts.append(np.vstack([np.zeros((20,44)),adj(f),u@adj(rest)]))
 lift=np.array(lifts);r0=Response(d,nb);hlocal=adj(lift)@r0.h@lift
 grid=np.zeros((mesh,mesh,44,44),complex)
 for ij,h in zip(index,hlocal):grid[tuple(ij)]=h*cross
 hr=np.fft.fft2(grid,axes=(0,1))/mesh**2
 diff=frac[None,:,:2]-frac[:,None,:2];dz=xyz[None,:,2]-xyz[:,None,2]
 for power in (3,5):
  for strain in (-.02,-.01,0,.01,.02):
   scaled=np.zeros_like(hr)
   for i,j in np.ndindex(mesh,mesh):
    images=np.array([(i+a,j+b) for a in (-mesh,0,mesh) for b in (-mesh,0,mesh)])
    dr=(images[:,None,None,:]+diff)@at.cell.array[:2,:2];rho2=(dr**2).sum(-1)
    mask=np.isclose(rho2,rho2.min(0),atol=1e-7,rtol=0);radius0=np.sqrt(rho2+dz**2);radius=np.sqrt((1+strain)**2*rho2+dz**2)
    ratio=np.divide(radius0,radius,out=np.ones_like(radius),where=radius>1e-9)**power
    scaled[i,j]=hr[i,j]*(ratio*mask).sum(0)/mask.sum(0)
   g=np.fft.ifft2(scaled,axes=(0,1))*mesh**2
   hs=np.array([g[tuple(ij)] for ij in index]);delta=lift@(hs-hlocal*cross)@adj(lift)
   herm=float(abs(delta-adj(delta)).max());assert herm<1e-7,herm
   dd=dict(d);dd['cell']=np.array(d['cell']).copy();dd['cell'][:2]*=1+strain;r=Response(dd,nb)
   ans=r.response(r.h+(delta+adj(delta))/2)
   row=dict(nb=nb,power=power,strain=strain,selected_cross_norm=float(np.linalg.norm(hs)),hermiticity_error=herm,**ans);out['rows'].append(row)
   (root/'prediction.json').write_text(json.dumps(out,indent=2))
print('PREDICTION_COMPLETE',json.dumps(out),flush=True)
