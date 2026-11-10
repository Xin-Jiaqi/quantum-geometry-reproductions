"""Read-only BiDB pilot inventory. No DFT, relaxation, remote access or jobs.
Export exact database structures; toolkit checks spatial constraints only.
"""
from pathlib import Path
import argparse, json, hashlib, csv, sys, subprocess
import numpy as np
import spglib
from ase.db import connect
from ase.io import write

p=argparse.ArgumentParser()
p.add_argument('--database',type=Path,required=True)
p.add_argument('--toolkit',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(a.toolkit/'src'))
from group_theory_operations import equivariant_map_basis
from group_theory_operations.stacking import partition_left_cosets
D=connect(str(a.database))
groups=[('Janus interface',7265,[7266,7268,7269,7270]),('TMD control',1929,[1930,1931]),('TMD control',2652,[2653,2654]),('III-VI control',862,[863,864]),('III-VI control',739,[740,741]),('anisotropic',783,[785]),('anisotropic',1082,[1084]),('centrosymmetric monolayer',2404,[2405,2406])]
old={7266,7269,7270,1930,864}

def dataset(atoms,tol):
 return spglib.get_symmetry_dataset((atoms.cell.array,atoms.get_scaled_positions(),atoms.numbers),symprec=tol)

def operations(atoms,tol):
 ds=dataset(atoms,tol);L=atoms.cell.array.T
 return ds,[L@r@np.linalg.inv(L) for r in ds.rotations]

def allowed(ops,layered):
 axial=[np.linalg.det(r)*r for r in ops]
 output=[np.kron(np.eye(2) if r[2,2]>0 else [[0,1],[1,0]],ax) if layered else ax for r,ax in zip(ops,axial)]
 inputs=[np.kron(r,ax) for r,ax in zip(ops,axial)]
 result={}
 for gate in [False,True]:
  ins=[m*(r[2,2] if gate else 1) for m,r in zip(inputs,ops)]
  basis=np.array(equivariant_map_basis([m.tolist() for m in output],[m.tolist() for m in ins]))
  for tensor in basis:
   for rout,rin in zip(output,ins):
    assert np.max(abs(rout@tensor-tensor@rin))<1e-8
  b=basis.reshape(-1,2 if layered else 1,3,3,3)
  for name,v in [('layer',b),('net',b.sum(axis=1,keepdims=True))]:
   result[('gate_derivative_' if gate else 'zero_gate_')+name]=[f'{"xyz"[i]}{ "xy"[j]}' for i in range(3) for j in range(2) if np.max(abs(v[:,:,i,j,2]),initial=0)>1e-8]
 return result

def connection(atoms):
 f=atoms.get_scaled_positions(wrap=False);cell=atoms.cell.array;n=len(f)//2
 order=np.argsort(atoms.positions[:,2]);lo=order[:n];hi=order[n:]
 names=atoms.numbers
 # Preferred representatives; nonunique modulo monolayer symmetry.
 choices=[('E',np.eye(3)),('m001',np.diag([1,1,-1])),('I',-np.eye(3)),('C2z',np.diag([-1,-1,1]))]
 parent=spglib.get_symmetry((cell,[[0,0,.5]],[1]),symprec=.001)
 choices += [('fractional_R',r) for r in parent['rotations']]
 for name,R in choices:
  if R[2,2]==0:continue
  for j in hi:
   if names[j]!=names[lo[0]]:continue
   t=f[j]-R@f[lo[0]];t[:2]%=1
   distances=[];assigned=[]
   for i in lo:
    d=f[hi]-(R@f[i]+t);d[:,:2]-=np.rint(d[:,:2]);cart=d@cell
    norm=np.linalg.norm(cart,axis=1);norm[names[hi]!=names[i]]=np.inf
    k=int(norm.argmin());distances.append(float(norm[k]));assigned.append(k)
   if max(distances)<.001 and len(set(assigned))==n:
    return {'representative':name,'fractional_rotation':R.tolist(),'inplane_translation_mod1':t[:2].tolist(),'vertical_affine_translation_A':float(t[2]*cell[2,2]),'center_separation_A':float(atoms.positions[hi,2].mean()-atoms.positions[lo,2].mean()),'max_mapping_error_A':max(distances),'nonunique_mod_monolayer_symmetry':True}
 raise ValueError('No tested lower-to-upper operation fits '+str(atoms.symbols))

rows=[];monos={}
for group,mid,ids in groups:
 m=D.get(id=mid);ma=m.toatoms();mds,mops=operations(ma,.001)
 monos[str(mid)]={'formula':m.formula,'uid':m.monolayer_uid,'metadata':m.key_value_pairs,'pointgroup_at_1e_3A':mds.pointgroup,'symmetry_by_tolerance':{str(t):dataset(ma,t).pointgroup for t in [1e-5,1e-4,1e-3,1e-2]},'spatial_allowance_Eip_Bz':allowed(mops,False)}
 sd=a.out/'structures';sd.mkdir(exist_ok=True)
 write(sd/f'mono-{mid}.vasp',ma,format='vasp',direct=True)
 for rid in ids:
  r=D.get(id=rid);at=r.toatoms();ds,ops=operations(at,.001);L=at.cell.array.T
  parent=spglib.get_symmetry((at.cell.array,[[0,0,.5]],[1]),symprec=.001)
  G=[L@x@np.linalg.inv(L) for x in parent['rotations']]
  row={'group':group,'parent_id':mid,'parent_formula':m.formula,'id':rid,'uid':r.bilayer_uid,'metadata':r.key_value_pairs,'pointgroup_at_1e_3A':ds.pointgroup,'symmetry_by_tolerance':{str(t):dataset(at,t).pointgroup for t in [1e-5,1e-4,1e-3,1e-2]},'layer_exchanging_ops':int(sum(x[2,2]<0 for x in ops)),'connection':connection(at),'spatial_allowance_Eip_Bz':allowed(ops,True),'existing_zero_gate_result':rid in old}
  try:row['point_orientation_coset_count']=len(partition_left_cosets(G,ops))
  except ValueError as e:row['point_orientation_coset_count']=None;row['coset_note']=str(e)
  assert r.slide_stability=='Stable' and not r.magnetic and r.gap_pbe>0
  write(sd/f'bilayer-{rid}.vasp',at,format='vasp',direct=True)
  rows.append(row)
out={'scope':'Proposed pilot only. No job submission or calculated new NMEE. Spatial allowance is necessary, not sufficient; symmetry at 1e-3A does not assert exact input symmetry.','database':{'local_path':str(a.database),'sha256':hashlib.sha256(a.database.read_bytes()).hexdigest(),'rows':D.count(),'source':'https://2dhub.org/bidb/bidb.html','license_note':'Structures/metadata from BiDB; cite Pakdel et al., Nat Commun 15,932 (2024).'},'toolkit_commit':subprocess.check_output(['git','-C',str(a.toolkit),'rev-parse','HEAD'],text=True).strip(),'versions':{'numpy':np.__version__,'spglib':spglib.__version__},'monolayers':monos,'bilayers':rows}
(a.out/'inventory.json').write_text(json.dumps(out,indent=2)+'\n')
with (a.out/'candidates.csv').open('w',newline='') as f:
 cols=['group','parent_formula','parent_id','id','uid','pointgroup','slide_stability','gap_pbe_eV','connection','tau1','tau2','mapping_error_A','zero_net','zero_layer','gate_net','existing_zero_gate_result','symmetry_caution']
 w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
 for r in rows:
  c=r['connection'];s=r['spatial_allowance_Eip_Bz'];tol=r['symmetry_by_tolerance']
  w.writerow({'group':r['group'],'parent_formula':r['parent_formula'],'parent_id':r['parent_id'],'id':r['id'],'uid':r['uid'],'pointgroup':r['pointgroup_at_1e_3A'],'slide_stability':r['metadata']['slide_stability'],'gap_pbe_eV':r['metadata']['gap_pbe'],'connection':c['representative'],'tau1':c['inplane_translation_mod1'][0],'tau2':c['inplane_translation_mod1'][1],'mapping_error_A':c['max_mapping_error_A'],'zero_net':','.join(s['zero_gate_net']),'zero_layer':','.join(s['zero_gate_layer']),'gate_net':','.join(s['gate_derivative_net']),'existing_zero_gate_result':r['existing_zero_gate_result'],'symmetry_caution':len(set(tol.values()))>1})
print('Exported',len(rows),'bilayers and',len(monos),'parents; no DFT or jobs.')
for r in rows:print(r['id'],r['parent_formula'],r['connection']['representative'],np.round(r['connection']['inplane_translation_mod1'],5),r['pointgroup_at_1e_3A'],r['spatial_allowance_Eip_Bz'])
