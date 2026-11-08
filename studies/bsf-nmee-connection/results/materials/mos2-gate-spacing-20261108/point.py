"""A bounded fixed-registry material reference and NMEE response, PBS only."""
from pathlib import Path
import hashlib,json,os,re,shutil,socket,subprocess,sys
import numpy as np
from ase.io import read,write
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute allocation required')
root=Path.cwd(); case=os.environ['CASE']; spacing=float(os.environ['DELTA_D']); gate=float(os.environ['GATE_VNM'])
assert re.fullmatch(r'[a-zA-Z0-9_.+-]+',case)
assert -.3 <= spacing <= 1.0 and -0.6 <= gate <= .6
base=root/case; base.mkdir(exist_ok=True); seed=base/'seed'; seed.mkdir(exist_ok=True)
if (base/'result.json').exists() or (seed/'INCAR').exists():
    raise SystemExit('Existing point; inspect manually instead of overwriting')
old=Path('/work/sunq/xjq/mos2-interface-20261108/A/seed')
atoms=read(old/'POSCAR'); original=atoms.positions.copy(); lower=atoms.positions[:,2]<atoms.cell[2,2]/2
atoms.positions[:,2]+=np.where(lower,-spacing/2,spacing/2)
assert np.max(abs(atoms.positions[:,:2]-original[:,:2]))==0
write(seed/'POSCAR',atoms,format='vasp',direct=True,vasp5=True,sort=False)
for name in ('POTCAR','KPOINTS'):shutil.copyfile(old/name,seed/name)
inc=(old/'INCAR').read_text()
inc=re.sub(r'^ICHARG=.*$', 'ICHARG=2', inc, flags=re.M)
inc=re.sub(r'^EFIELD=.*$',f'EFIELD={-gate/10:.10f}',inc,flags=re.M)
inc=re.sub(r'^SYSTEM=.*$',f'SYSTEM=3R MoS2 fixed registry dd={spacing} Fz={gate}',inc,flags=re.M)
(seed/'INCAR').write_text(inc)
meta={'case':case,'delta_d_A':spacing,'physical_Fz_Vnm':gate,'VASP_EFIELD_VA':-gate/10,
      'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'lower_atoms_zero_based':np.where(lower)[0].tolist(),
      'interface_gap_A':float(atoms.positions[~lower,2].min()-atoms.positions[lower,2].max()),
      'model':'full KS static spin-only EB response; gate SCF, probes frozen potential',
      'input_sha256':{n:hashlib.sha256((seed/n).read_bytes()).hexdigest() for n in ('POSCAR','INCAR','KPOINTS')}}
(base/'input.json').write_text(json.dumps(meta,indent=2))
ranks=len(Path(os.environ['PBS_NODEFILE']).read_text().splitlines())
with (seed/'vasp.log').open('w') as f:
    subprocess.run(['mpirun','-np',str(ranks),'-ppn',str(ranks),'-hostfile',os.environ['PBS_NODEFILE'],
     '/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'],cwd=seed,stdout=f,stderr=subprocess.STDOUT,check=True)
out=(seed/'OUTCAR').read_text(errors='replace')
if 'General timing' not in out or 'aborting loop because EDIFF is reached' not in out:
    raise RuntimeError('SCF incomplete')
(seed/'finished.json').write_text(json.dumps(meta,indent=2))
env=dict(os.environ,STRUCTURE=case)
subprocess.run([sys.executable,'-u','export.py'],env=env,check=True)
from analyze import Response
d=np.load(base/'cache.npz'); results={}
for nb in (96,112):
    r=Response(d,nb)
    results[str(nb)]={'standard':r.response(b=.002),'half_b':r.response(b=.001)}
meta.update(results=results,units='1e-9 (nm V T)^-1, Pauli spin density')
(base/'result.json').write_text(json.dumps(meta,indent=2))
print('POINT_DONE',json.dumps(meta),flush=True)
