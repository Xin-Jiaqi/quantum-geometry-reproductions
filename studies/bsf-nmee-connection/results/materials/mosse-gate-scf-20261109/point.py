"""Independent gate-SCF at fixed MoSSe S/S geometry; all work inside PBS.

Perpendicular gate is self-consistent. Subsequent spin-only EB susceptibility
uses the existing frozen-potential full-KS response definition, without ionic
motion or self-consistent exchange feedback to the probe magnetic field.
"""
from pathlib import Path
import os,socket,json,re,shutil,subprocess,sys,hashlib,time
import numpy as np
from ase.io import read,write
from analyze import Response

if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute node required')
root=Path.cwd();field=float(os.environ['FIELD']);case=os.environ['CASE']
assert (case,field) in [('fm010',-.1),('fp010',.1)]
base=root/case;base.mkdir(exist_ok=True)
assert not (base/'progress.json').exists(),'Do not overwrite a started case'
source=Path('/work/sunq/xjq/mosse-strain-20261108/e000')
exe='/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'
ranks=len(Path(os.environ['PBS_NODEFILE']).read_text().splitlines())
assert ranks==4
meta={'case':case,'F_external_V_nm':field,'VASP_EFIELD_eV_A':-field/10,
      'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'source':str(source),
      'geometry':'fixed existing relaxed e000 A; independently matched lower and upper isolated layers',
      'scope':__doc__,'materials':{}}
def progress(stage):
    meta['stage']=stage;(base/'progress.json').write_text(json.dumps(meta,indent=2))
    print('PROGRESS',case,stage,flush=True)

def planar(path):
    with path.open() as f:
        header=[next(f) for _ in range(8)]
        for _ in range(sum(map(int,header[6].split()))):next(f)
        nx,ny,nz=map(int,next(x for x in f if x.strip()).split())
        vals=[]
        while len(vals)<nx*ny*nz:vals.extend(map(float,next(f).split()))
    assert len(vals)==nx*ny*nz
    return np.asarray(vals).reshape(nz,ny,nx).mean(axis=(1,2))

at=read(source/'A/seed/POSCAR');at.set_constraint()
lower=at.positions[:,2]<at.cell[2,2]/2;assert lower.sum()==3
progress('preparing')
for material in ('A','mono','upper'):
    seed=base/material/'seed';seed.mkdir(parents=True,exist_ok=True)
    atoms=at.copy() if material=='A' else at[lower if material=='mono' else ~lower]
    if material!='A':
        zmo=atoms.positions[[a.symbol=='Mo' for a in atoms],2][0]
        atoms.positions[:,2]+=atoms.cell[2,2]/2-zmo
    if material=='A':shutil.copyfile(source/'A/seed/POSCAR',seed/'POSCAR')
    elif material=='mono':
        exact=read(source/'mono/seed/POSCAR')
        assert abs(atoms.positions-exact.positions).max()<1e-7
        shutil.copyfile(source/'mono/seed/POSCAR',seed/'POSCAR')
    else:write(seed/'POSCAR',atoms,format='vasp',direct=True,vasp5=True,sort=False)
    sr=source/('A' if material=='A' else 'mono')/'seed'
    for name in ('POTCAR','KPOINTS'):shutil.copyfile(sr/name,seed/name)
    text=(sr/'INCAR').read_text()
    text=re.sub(r'^SYSTEM=.*$',f'SYSTEM=MoSSe {material} fixed ions external F={field:+.3f} V/nm',text,flags=re.M)
    text=re.sub(r'^EFIELD=.*$',f'EFIELD={-field/10:.8f}',text,flags=re.M)
    if material!='upper':
        for name in ('WAVECAR','CHGCAR'):shutil.copyfile(sr/name,seed/name)
        text=text.replace('ISTART=0','ISTART=1').replace('ICHARG=2','ICHARG=1')
    # ICHARG=1 updates density self-consistently; no fixed-charge ICHARG=11.
    assert 'BEXT=0 0 0' in text and 'NSW=0' in text and 'EDIFF=1e-9' in text
    (seed/'INCAR').write_text(text)
    progress(material+'_scf_running');t0=time.monotonic()
    with (seed/'vasp.log').open('w') as log:
        subprocess.run(['mpirun','-np',str(ranks),'-ppn',str(ranks),'-hostfile',os.environ['PBS_NODEFILE'],exe],cwd=seed,stdout=log,stderr=subprocess.STDOUT,check=True)
    out=(seed/'OUTCAR').read_text(errors='replace');log=(seed/'vasp.log').read_text(errors='replace')
    if 'General timing' not in out or 'aborting loop because EDIFF is reached' not in out or re.search(r'^\s*(DAV|RMM):\s*240\s',log,re.M):
        raise RuntimeError('SCF did not converge/terminate normally')
    row={'scf_seconds':time.monotonic()-t0,'electronic_iterations':len(re.findall(r'^\s*(DAV|RMM):',log,re.M)),
         'input_sha256':{n:hashlib.sha256((seed/n).read_bytes()).hexdigest() for n in ('POSCAR','INCAR','KPOINTS')},
         'positions_A':atoms.positions.tolist(),'cell_A':atoms.cell.array.tolist(),
         'symbols':atoms.get_chemical_symbols(),'F_external_V_nm':field,'VASP_EFIELD_eV_A':-field/10,
         'scf_converged':True,'fixed_ions':True,'job':meta['job'],'node':meta['node']}
    (seed/'finished.json').write_text(json.dumps(row,indent=2));meta['materials'][material]=row
    # LOCPOT includes the saved scalar potential. Vacuum slopes are checked later,
    # rather than assuming its convention is identical to the input field.
    row['planar_potential_eV']=planar(seed/'LOCPOT').tolist()
    row['planar_charge_raw']=planar(seed/'CHGCAR').tolist()
    progress(material+'_export_running')
    subprocess.run([sys.executable,str(root/'export.py')],cwd=base,env=dict(os.environ,STRUCTURE=material),check=True)
    data=np.load(base/material/'cache.npz');res={}
    for nb in ((96,112) if material=='A' else (48,56,64)):
        r=Response(data,nb);res[str(nb)]={'standard':r.response(b=.002),'half_b':r.response(b=.001)}
    assert min(x['standard']['gap_eV'] for x in res.values())>.05,'Insulating response invalid'
    row.update(response=res,export=json.loads((base/material/'export.json').read_text()))
    (base/material/'response.json').write_text(json.dumps({'results':res,'units':'1e-9 (nm V T)^-1'},indent=2))
    progress(material+'_complete')
progress('complete');meta['complete']=True
(base/'result.json').write_text(json.dumps(meta,indent=2));print('CASE_COMPLETE',case,flush=True)
