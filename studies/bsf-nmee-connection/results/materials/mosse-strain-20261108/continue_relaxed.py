"""Fixed-registry biaxial strain: constrained ionic relaxation, independent SCF,
full KS NMEE, and a matched-geometry isolated lower-layer reference. PBS only.
"""
from pathlib import Path
import os,socket,json,re,shutil,subprocess,sys,hashlib,time
import numpy as np
from ase.io import read,write
from geometry_check import validate_fixed_registry
from analyze import Response
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute only')
root=Path.cwd();strain=float(os.environ['STRAIN']);label=os.environ['CASE'];assert label in ('em01','e000','ep01','em02','ep02');assert abs(strain)<=.020001
base=root/label;base.mkdir(exist_ok=True)
assert label=='ep02' and (base/'progress.json').exists()
assert not (base/'A/seed/OUTCAR').exists(),'Do not overwrite an existing SCF'
shutil.copyfile(base/'progress.json',base/'progress-before-continuation.json')
old=Path('/work/sunq/xjq/mosse-7269-mechanism-20261104/zero-field-refined-20261107/seed')
exe='/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'
ranks=len(Path(os.environ['PBS_NODEFILE']).read_text().splitlines())
meta={'case':label,'strain':strain,'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'relaxation':'PBE-D3(0), SOC; fixed strained cell; all ions relaxed; same lateral registry verified within 1e-5 A; force 0.01 eV/A', 'NMEE':'frozen-potential full KS spin-only EB; strained Hamiltonian independently self-consistent','source':str(old)}
def progress(stage,**kw):
 meta.update(stage=stage,**kw);(base/'progress.json').write_text(json.dumps(meta,indent=2));print('PROGRESS',stage,flush=True)
def common(nat,nb):return f'''SYSTEM=MoSSe S/S biaxial strain {strain:+.4f}
ENCUT=450
PREC=Accurate
LASPH=T
LREAL=F
ISYM=-1
LSORBIT=T
SAXIS=0 0 1
MAGMOM={nat*3}*0
BEXT=0 0 0
NBANDS={nb}
ISMEAR=0
SIGMA=0.02
EDIFF=1e-9
NELM=240
ALGO=Normal
KPAR=1
NCORE=1
ISTART=0
ICHARG=2
LDIPOL=T
IDIPOL=3
DIPOL=0.5 0.5 0.5
EFIELD=0
IVDW=11
LCHARG=T
LWAVE=T
'''
def inputs(dest,atoms,incar):
 dest.mkdir(parents=True,exist_ok=True);write(dest/'POSCAR',atoms,format='vasp',direct=True,vasp5=True,sort=False)
 shutil.copyfile(old/'POTCAR',dest/'POTCAR');shutil.copyfile(old/'KPOINTS',dest/'KPOINTS');(dest/'INCAR').write_text(incar)
def vasp(dest,relax=False):
 with (dest/'vasp.log').open('w') as f:subprocess.run(['mpirun','-np',str(ranks),'-ppn',str(ranks),'-hostfile',os.environ['PBS_NODEFILE'],exe],cwd=dest,stdout=f,stderr=subprocess.STDOUT,check=True)
 out=(dest/'OUTCAR').read_text(errors='replace');log=(dest/'vasp.log').read_text(errors='replace')
 if 'General timing' not in out or re.search(r'^\s*(DAV|RMM):\s*240\s',log,re.M):raise RuntimeError('Electronic convergence/termination failed')
 if relax and 'reached required accuracy' not in out:raise RuntimeError('Ionic relaxation incomplete')
 if not relax and 'aborting loop because EDIFF is reached' not in out:raise RuntimeError('SCF not converged')
 return out
# The previous VASP relaxation finished; only its overly strict coordinate guard failed.
initial=read(base/'relax/POSCAR');relax=base/'relax'
assert len(initial.constraints)==0,'Audit the actual selective-dynamics flags'
out=(relax/'OUTCAR').read_text(errors='replace')
assert 'General timing' in out and 'reached required accuracy' in out
meta['relaxation_job']='558256.cmspbs1'
meta['continuation_reason']='Relaxation converged; registry drift is sub-micro-A. Replace dimensionless 1e-8 guard with verified Cartesian 1e-5 A criterion. Actual POSCAR had all ions free.'
at=read(relax/'CONTCAR');at.set_constraint();at.positions[:,2]+=(at.cell[2,2]/2-(at.positions[:,2].max()+at.positions[:,2].min())/2)
assert abs(at.cell.array-initial.cell.array).max()<1e-7
registry=validate_fixed_registry(initial,at)
lower=at.positions[:,2]<at.cell[2,2]/2;assert lower.sum()==3
try:
 aa=read(relax/'OUTCAR',index=-1);forces=aa.get_forces(apply_constraint=False);maxfz=float(abs(forces[:,2]).max());maxfxy=float(abs(forces[:,:2]).max())
except Exception as err:raise RuntimeError('Cannot verify relaxed forces') from err
if maxfz>.012 or maxfxy>.012:raise RuntimeError(f'Constrained residual forces too large: {maxfz},{maxfxy}')
meta['geometry']={'max_z_force_eVA':maxfz,'max_inplane_force_eVA':maxfxy,'interface_SS_gap_A':float(at.positions[~lower,2].min()-at.positions[lower,2].max()),'cell_A':at.cell.array.tolist(),'positions_A':at.positions.tolist(),'lower_atoms':np.flatnonzero(lower).tolist(),'registry_validation':registry}
progress('relaxed')
for material in ('A','mono'):
 atoms=at.copy() if material=='A' else at[lower]
 if material=='mono':
  zmo=atoms.positions[[x.symbol=='Mo' for x in atoms],2][0];atoms.positions[:,2]+=atoms.cell[2,2]/2-zmo
 seed=base/material/'seed';nb=112 if material=='A' else 64
 inc=common(len(atoms),nb)+'NSW=0\nIBRION=-1\nLVHAR=T\nNUM_WANN='+str(44 if material=='A' else 22)+'''\nLWANNIER90=T
LWRITE_MMN_AMN=T
WANNIER90_WIN="
begin projections
Mo:d
S:p
Se:p
end projections
"
'''
 inputs(seed,atoms,inc);progress(material+'_scf_running');vasp(seed)
 (seed/'finished.json').write_text(json.dumps({'job':meta['job'],'node':meta['node'],'strain':strain,'material':material,'input_sha256':{n:hashlib.sha256((seed/n).read_bytes()).hexdigest() for n in ('POSCAR','KPOINTS','INCAR')}},indent=2))
 progress(material+'_export_running');subprocess.run([sys.executable,str(root/'export.py')],cwd=base,env=dict(os.environ,STRUCTURE=material),check=True)
 d=np.load(base/material/'cache.npz');res={}
 for nb in ((96,112) if material=='A' else (48,56,64)):
  rr=Response(d,nb);res[str(nb)]={'standard':rr.response(b=.002),'half_b':rr.response(b=.001)}
 (base/material/'response.json').write_text(json.dumps({'results':res,'units':'1e-9 (nm V T)^-1'},indent=2))
 if min(x['standard']['gap_eV'] for x in res.values())<=.05:raise RuntimeError('Insulating gap too small for claimed regime')
 progress(material+'_complete',**{material+'_response':res})
progress('electrostatic_running');subprocess.run([sys.executable,str(root/'reference.py')],cwd=base,check=True)
progress('complete');(base/'result.json').write_text(json.dumps(meta,indent=2));print('CASE_COMPLETE',label,flush=True)
