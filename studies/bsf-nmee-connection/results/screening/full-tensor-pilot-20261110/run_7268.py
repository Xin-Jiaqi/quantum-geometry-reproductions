"""One bounded fixed-ion 7268 reference state per PBS job; no auto retry."""
from pathlib import Path
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import numpy as np
from tensor_response import TensorResponse
from validate_7269 import serial


def main():
    if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
        raise SystemExit('PBS compute required')
    root=Path.cwd()
    acceptance=json.loads((root/'final_anchor_acceptance.json').read_text())
    if not acceptance['passed']:raise RuntimeError('7269 acceptance required')
    case=os.environ['CASE']
    if case not in ('A0','lower0','upper0'):raise ValueError('Only authorized zero-gate pilot states')
    label=case[:-1]
    base=root/'pilot7268'/case
    seed=base/'seed'
    if base.exists():raise RuntimeError('Existing case requires audit; no automatic retry/overwrite')
    seed.mkdir(parents=True)
    source=Path('/work/sunq/xjq/mosse-strain-20261108/e000')/('A' if label=='A' else 'mono')/'seed'
    shutil.copyfile(root/'inputs7268'/f'POSCAR.{label}',seed/'POSCAR')
    for f in ('POTCAR','KPOINTS'):shutil.copyfile(source/f,seed/f)
    text=(source/'INCAR').read_text()
    for key,value in {'SYSTEM':f'MoSSe7268 {case} fixed ions Fz=0','ISTART':'0','ICHARG':'2','EFIELD':'0','KPAR':'1','NCORE':'1'}.items():
        pattern=rf'^{key}\s*=.*$'
        text=re.sub(pattern,f'{key}={value}',text,flags=re.M) if re.search(pattern,text,re.M) else text+f'\n{key}={value}\n'
    for required in ('LSORBIT=T','BEXT=0 0 0','NSW=0','EDIFF=1e-9','ENCUT=450'):
        if required not in text:raise ValueError('Inherited protocol mismatch: '+required)
    (seed/'INCAR').write_text(text)
    ranks=len(Path(os.environ['PBS_NODEFILE']).read_text().splitlines())
    if ranks!=4:raise ValueError('Budget requires4MPI ranks')
    rec={'case':case,'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'complete':False,'fixed_ions':True,'Fz_V_nm':0,'stage':'SCF'}
    path=base/'result.json'
    def save():path.write_text(json.dumps(rec,default=serial,indent=2))
    save();start=time.monotonic()
    with (seed/'vasp.log').open('w') as log:
        subprocess.run(['mpirun','-np','4','-ppn','4','-hostfile',os.environ['PBS_NODEFILE'],'/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'],cwd=seed,stdout=log,stderr=subprocess.STDOUT,check=True)
    rec['scf_seconds']=time.monotonic()-start
    output=(seed/'OUTCAR').read_text(errors='replace');log=(seed/'vasp.log').read_text(errors='replace')
    if 'General timing' not in output or 'aborting loop because EDIFF is reached' not in output or re.search(r'^\s*(DAV|RMM):\s*240\s',log,re.M):raise RuntimeError('SCF did not converge')
    rec['electronic_iterations']=len(re.findall(r'^\s*(DAV|RMM):',log,re.M))
    from ase.io import read
    state=read(seed/'OUTCAR',index=-1)
    rec['max_force_eV_A']=float(np.max(np.linalg.norm(state.get_forces(),axis=1)))
    rec['forces_eV_A']=state.get_forces().tolist()
    rec['stage']='export';save()
    subprocess.run([sys.executable,str(root/'export_tensor.py'),str(seed),str(base),'--layer-mode','bilayer' if label=='A' else 'mono','--mesh','6','6','--cut','0.5'],check=True)
    rec['export']=json.loads((base/'export.json').read_text())
    rec['stage']='response';save();t=time.monotonic()
    data=np.load(base/'tensor_cache.npz');rec['responses']={}
    for nb in ((96,112) if label=='A' else (48,56,64)):
        r=TensorResponse(data,nb)
        rec['responses'][str(nb)]={'standard':r.response(b=.002),'half_b':r.response(b=.001)}
    rec['response_seconds']=time.monotonic()-t
    rec['elapsed_seconds']=time.monotonic()-start
    rec['stage']='completed_requires_scientific_review';rec['complete']=True
    rec['physical_readout']='total,lower,upper' if label=='A' else 'total only; halfspaces are not material layers'
    save()


if __name__=='__main__':main()
