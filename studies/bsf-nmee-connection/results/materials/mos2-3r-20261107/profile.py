"""Real-space response profile and corrected seed-gap metadata. PBS only."""
from pathlib import Path
import json,os,socket
import numpy as np
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('compute allocation required')
root=Path.cwd();labels=os.environ.get('LABELS','mono,A,B').split(',')

def plane(d):
 with (d/'CHGCAR').open() as f:
  head=[next(f) for _ in range(8)];nat=sum(map(int,head[6].split()))
  for _ in range(nat):next(f)
  dims=next(s for s in f if s.strip()).split();nx,ny,nz=map(int,dims)
  for ax in range(2):
   if ax:
    for line in f:
     if line.split()==dims:break
    else:raise ValueError('missing spin grid')
   vals=[]
   while len(vals)<nx*ny*nz:vals.extend(map(float,next(f).split()))
   p=np.asarray(vals).reshape(nz,ny,nx).mean(axis=(1,2))/nz
   if ax==0 and abs(p.sum()-(24 if nat==3 else 48))>.01:raise ValueError('charge normalization')
  return p

for label in labels:
 for d in sorted(root.glob(f'{label}_F*_k*_e*_b*')):
  if not (d/'answer.json').exists():continue
  r=json.loads((d/'answer.json').read_text())
  gate=r['gate_physical_Vnm']
  if not (d/'answer_raw.json').exists():(d/'answer_raw.json').write_text(json.dumps(r,indent=2))
  signs=r.get('B_signs',[1,-1]);e=r['probe_VA'];b=r['zeeman_eV'];area=r['area_A2']
  seed_gaps={}
  for bs in signs:
   seed=d/f'seedB{bs:+d}';g=json.loads((seed/'finished.json').read_text())['gap_eV']
   if not 0<g<100:raise ValueError('invalid zero-probe reference gap')
   seed_gaps[f'B{bs:+d}']=g
  r['raw_finite_field_eigenvalue_gaps_not_physical']=r.pop('gaps_eV',r.get('raw_finite_field_eigenvalue_gaps_not_physical',{}))
  r['seed_gaps_eV']=seed_gaps
  r['gap_note']='Finite-PEAD empty EIGENVAL values are 10000-eV placeholders, not physical gaps. Only normal-diagonalized zero-Ex seed gaps are reported. Response densities and C are unchanged.'
  cell=np.array([list(map(float,s.split())) for s in (d/'start/POSCAR').read_text().splitlines()[2:5]])
  r['electric_probe_to_conservative_breakdown_ratio']=float(10*r['mesh']*max(abs(cell[:2,0]*e))/min(seed_gaps.values()))
  if r['electric_probe_to_conservative_breakdown_ratio']>=1:raise ValueError('probe exceeds conservative seed-gap bound')
  mu=2.00231930436256*5.7883818060e-5/2
  cp=sum(bs*(plane(d/f'E+1B{bs:+d}')-plane(d/f'E-1B{bs:+d}')) for bs in signs)/len(signs)*(1e10*mu/area)/(2*e*b)
  nz=len(cp);L=cell[2,2];z=np.arange(nz)*L/nz-L/2;weights=np.zeros(nz);weights[1:nz//2]=1;weights[[0,nz//2]]=.5
  intg=np.array([cp.sum(),weights@cp,(1-weights)@cp])
  np.testing.assert_allclose(intg,r['C_xyz_total_bottom_top_1e-9'][0],atol=1e-6,rtol=1e-7)
  positions=np.array([list(map(float,s.split()[:3])) for s in (d/'start/POSCAR').read_text().splitlines()[8:]])
  out={'structure':label,'gate_physical_Vnm':gate,'z_A':z.tolist(),'response_per_A':(cp/(L/nz)).tolist(),'slice_response':cp.tolist(),'total_bottom_top':intg.tolist(),'first_z_moment_C_A':float(z@cp),'horizontal_mirror_oddness_relative_L2':float(np.linalg.norm(cp+cp[(-np.arange(nz))%nz])/np.linalg.norm(cp)),'ion_z_A':(positions[:,2]*L-L/2).tolist(),'analysis_pbs':os.environ['PBS_JOBID'],'note':'Monolayer halves describe spatial cancellation within one layer, not two material layers. Profile is the planar average of the induced Pauli spin response.'}
  (d/'profile.json').write_text(json.dumps(out,indent=2));(d/'answer.json').write_text(json.dumps(r,indent=2))
  print(label,gate,'C',intg.tolist(),'first_z_moment',out['first_z_moment_C_A'],'oddness',out['horizontal_mirror_oddness_relative_L2'],flush=True)
