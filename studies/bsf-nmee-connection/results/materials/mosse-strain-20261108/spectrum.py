"""Post-DFT spectral counterfactual, not an additive energy/geometry decomposition.
Swap energy-ordered spectra between independently strained states only when the
within-k degeneracy blocks are compatible. Keep physical spin and position
operators belonging to the state reference; never substitute old velocities.
"""
from pathlib import Path
import os,socket,json
import numpy as np
from analyze import Response
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute only')
root=Path.cwd();labels=tuple(os.getenv('CASES','em02:e000:ep02').split(':'));assert 'e000' in labels and set(labels)<=set(('em02','e000','ep02'))
datasets={c:np.load(root/c/'A/cache.npz') for c in labels}
output=root/('spectrum.json' if len(labels)==3 else 'spectrum_partial.json')
strain={'em02':-.02,'e000':0.,'ep02':.02}
out={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'units':'1e-9 (nm V T)^-1','scope':__doc__,'pairing':'energy order at identical fractional k; reject if a <1e-5 eV degeneracy in either spectrum is split by >1e-4 eV in the other','rows':[],'layer_electrons':{}}
for state,d in datasets.items():
 assert np.max(abs(d['k']-datasets['e000']['k']))<1e-8
 out['layer_electrons'][state]=float(np.diagonal(d['P'],axis1=-2,axis2=-1)[:,:48].real.sum(1).mean())
 for nb in (96,112):
  r=Response(d,nb);es=d['energies'][:,:nb];ds=np.diff(es,axis=1)
  for spectrum,other in datasets.items():
   en=other['energies'][:,:nb];de=np.diff(en,axis=1)
   bad=((ds<1e-5)&(de>1e-4))|((de<1e-5)&(ds>1e-4))
   row={'states':state,'spectrum':spectrum,'nb':nb,'degeneracy_mismatch_count':int(bad.sum())}
   if bad.any(): row['invalid_control']='Incompatible degeneracy blocks; not used for mechanism inference'
   else:
    paired=en.copy()
    if state!=spectrum:
     for ik in range(len(en)):
      boundaries=[0]+(np.flatnonzero((ds[ik]>=1e-5)&(de[ik]>=1e-5))+1).tolist()+[nb]
      for lo,hi in zip(boundaries[:-1],boundaries[1:]):
       if hi-lo>1: paired[ik,lo:hi]=paired[ik,lo:hi].mean()
    row['max_degenerate_energy_adjustment_eV']=float(np.max(abs(paired-en)))
    h=paired[:,:,None]*np.eye(nb);ans=r.response(h);row.update(ans)
    row['metric_normalized_lower']=ans['C_total_lower_upper'][1]*(1+strain[state])
   out['rows'].append(row)
   output.write_text(json.dumps(out,indent=2));print(state,spectrum,nb,row.get('C_total_lower_upper',row.get('invalid_control')),flush=True)
