"""Zero-field references at the unchanged 3R geometry, executed on PBS only."""
from pathlib import Path
import json, os, re, shutil, socket, subprocess
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute allocation required')
root=Path.cwd(); label=os.environ['STRUCTURE']; assert label in ('A','mono')
old=Path('/work/sunq/xjq/mos2-nmee-20261107')
gate=float(os.environ.get('GATE_VNM','0'));case=label if gate==0 else f'{label}_F{gate:+.3f}'
base=root/case; base.mkdir(exist_ok=True); d=base/'seed'; d.mkdir(exist_ok=True)
if (d/'finished.json').exists():
    raise SystemExit('Reference already finished; never overwrite')
nb,nat,nwann,nocc=(112,6,44,48) if label=='A' else (64,3,22,24)
for name in ('POSCAR','POTCAR','KPOINTS','CHGCAR'):
    shutil.copyfile(old/f'{label}_F+0.000_k6_e0.004_b0.016/start'/name,d/name)
inc=f'''SYSTEM=MoS2 {label} zero-field interface control
ENCUT=450
PREC=Accurate
LASPH=T
LREAL=F
ISYM=-1
LSORBIT=T
SAXIS=0 0 1
MAGMOM={3*nat}*0
BEXT=0 0 0
NBANDS={nb}
ISMEAR=0
SIGMA=0.02
EDIFF=1e-9
NELM=240
NSW=0
IBRION=-1
LCHARG=T
LWAVE=T
ALGO=Normal
KPAR=1
NCORE=1
ISTART=0
ICHARG=1
LDIPOL=T
IDIPOL=3
DIPOL=0.5 0.5 0.5
EFIELD={-gate/10:.10f}
LVHAR=T
NUM_WANN={nwann}
LWANNIER90=T
LWRITE_MMN_AMN=T
WANNIER90_WIN="
begin projections
Mo:d
S:p
end projections
"
'''
(d/'INCAR').write_text(inc)
ranks=len(Path(os.environ['PBS_NODEFILE']).read_text().splitlines())
with (d/'vasp.log').open('w') as f:
    subprocess.run(['mpirun','-np',str(ranks),'-ppn',str(ranks),'-hostfile',os.environ['PBS_NODEFILE'],
      '/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'],cwd=d,stdout=f,stderr=subprocess.STDOUT,check=True)
out=(d/'OUTCAR').read_text(errors='replace');log=(d/'vasp.log').read_text(errors='replace')
if 'General timing' not in out or 'aborting loop because EDIFF is reached' not in out or re.search(r'^\s*DAV:\s*240\s',log,re.M):
    raise RuntimeError('Unconverged SCF')
for name in ('wannier90.mmn','WAVECAR','LOCPOT','CHGCAR'):
    assert (d/name).stat().st_size>1000,name
(d/'finished.json').write_text(json.dumps({'structure':label,'job':os.environ['PBS_JOBID'],
    'node':socket.gethostname(),'bands':nb,'occupied':nocc,'electrostatic_potential':'ionic+Hartree, LVHAR',
    'zero_magnetic_field':True,'gate_physical_Vnm':gate,'dipole_corrected':True},indent=2))
print('REFERENCE_DONE',label,flush=True)
