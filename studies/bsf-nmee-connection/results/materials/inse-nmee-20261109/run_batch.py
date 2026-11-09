"""Fixed-geometry InSe: gate-SCF then full-KS insulating spin response. PBS only."""
from pathlib import Path
import os,socket,json,re,shutil,subprocess,sys,time,hashlib
import numpy as np
from ase.io import read
from analyze import Response
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute node required')
root=Path.cwd();batch=os.environ['BATCH'];assert batch in ('mono','bilayer')
points=([('mono0','mono',0),('monom','mono',-.1),('monop','mono',.1)] if batch=='mono' else
        [('A0','A',0),('B0','B',0),('Am','A',-.1),('Ap','A',.1)])
exe='/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'
pots=Path('/work/sunq/mxy/bin/potpaw_PBE.64')
ranks=len(Path(os.environ['PBS_NODEFILE']).read_text().splitlines());assert ranks==4
for name,label,field in points:
    base=root/name;seed=base/label/'seed'
    if (base/'result.json').exists():
        prior=json.loads((base/'result.json').read_text());assert prior['complete'];continue
    assert not (seed/'OUTCAR').exists(),f'Existing unfinished SCF {name}; audit before restart'
    seed.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/f'POSCAR.{label}',seed/'POSCAR')
    atoms=read(seed/'POSCAR');single=label=='mono';nb=80 if single else 160
    nwan=28 if single else 56
    with (seed/'POTCAR').open('wb') as f:
        for element in ('In_d','Se'):f.write((pots/element/'POTCAR').read_bytes())
    (seed/'KPOINTS').write_text('full Gamma mesh\n0\nGamma\n6 6 1\n0 0 0\n')
    initial='ISTART=0\nICHARG=2'
    if field:
        prev=root/('mono0' if single else 'A0')/label/'seed'
        for file in ('WAVECAR','CHGCAR'):shutil.copyfile(prev/file,seed/file)
        initial='ISTART=1\nICHARG=1'
    incar=f'''SYSTEM=InSe {name} physical perpendicular F={field:+.3f} V/nm
ENCUT=450
PREC=Accurate
LASPH=T
LREAL=F
ISYM=-1
LSORBIT=T
SAXIS=0 0 1
MAGMOM={len(atoms)*3}*0
BEXT=0 0 0
NBANDS={nb}
ISMEAR=0
SIGMA=0.02
EDIFF=1e-9
NELM=240
ALGO=Normal
KPAR=4
NCORE=1
{initial}
LDIPOL=T
IDIPOL=3
DIPOL=0.5 0.5 0.5
EFIELD={-field/10:.8f}
IVDW=11
LCHARG=T
LWAVE=T
NSW=0
IBRION=-1
LVHAR=T
NUM_WANN={nwan}
LWANNIER90=T
LWRITE_MMN_AMN=T
WANNIER90_WIN="
begin projections
In:s;p
Se:p
end projections
"
'''
    (seed/'INCAR').write_text(incar)
    record={'name':name,'label':label,'F_external_V_nm':field,'job':os.environ['PBS_JOBID'],
        'node':socket.gethostname(),'complete':False,'scope':__doc__,'fixed_ions':True,
        'input_sha256':{f:hashlib.sha256((seed/f).read_bytes()).hexdigest() for f in ('POSCAR','KPOINTS','INCAR')},
        'geometry':{'positions_A':atoms.positions.tolist(),'cell_A':atoms.cell.array.tolist(),'symbols':atoms.get_chemical_symbols()}}
    def progress(stage):
        record['stage']=stage;(base/'progress.json').write_text(json.dumps(record,indent=2));print(name,stage,flush=True)
    progress('SCF_running');start=time.monotonic()
    with (seed/'vasp.log').open('w') as f:
        subprocess.run(['mpirun','-np',str(ranks),'-ppn',str(ranks),'-hostfile',os.environ['PBS_NODEFILE'],exe],cwd=seed,stdout=f,stderr=subprocess.STDOUT,check=True)
    out=(seed/'OUTCAR').read_text(errors='replace');log=(seed/'vasp.log').read_text(errors='replace')
    if 'General timing' not in out or 'aborting loop because EDIFF is reached' not in out or re.search(r'^\s*(DAV|RMM):\s*240\s',log,re.M):
        raise RuntimeError('SCF convergence/termination failed')
    record['scf_seconds']=time.monotonic()-start
    record['electronic_iterations']=len(re.findall(r'^\s*(DAV|RMM):',log,re.M))
    record['potcar_titles']=re.findall(r'TITEL\s*=\s*(.*)',(seed/'POTCAR').read_text())
    progress('export_running')
    subprocess.run([sys.executable,str(root/'export.py')],cwd=base,env=dict(os.environ,STRUCTURE=label),check=True)
    d=np.load(base/label/'cache.npz');record['export']=json.loads((base/label/'export.json').read_text())
    results={}
    for cutoff in ((64,72,80) if single else (128,144,160)):
        rr=Response(d,cutoff);results[str(cutoff)]={'standard':rr.response(b=.002),'half_b':rr.response(b=.001)}
    record['response']=results
    if min(x['standard']['gap_eV'] for x in results.values())<.05:raise RuntimeError('Insulating gap <.05 eV')
    record['complete']=True;progress('complete');(base/'result.json').write_text(json.dumps(record,indent=2))
print('BATCH_COMPLETE',batch,flush=True)
