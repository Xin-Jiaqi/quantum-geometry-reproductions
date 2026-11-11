"""One zero-field WS2 bilayer. All VASP and response work runs in one PBS allocation."""
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import socket
import subprocess
import sys
import time

import numpy as np
from ase.io import read


def main():
    nodefile = Path(os.environ['PBS_NODEFILE'])
    node = socket.gethostname()
    if not os.environ.get('PBS_JOBID') or not node.startswith('csuper') or len(nodefile.read_text().splitlines()) != 4:
        raise RuntimeError('Four-rank PBS compute allocation required')
    here = Path.cwd()
    kernels = here.parent
    contract = json.loads((here / 'contract.json').read_text())
    if hashlib.sha256((here / 'POSCAR').read_bytes()).hexdigest() != contract['POSCAR_sha256']:
        raise ValueError('Unexpected bilayer structure')
    if shutil.disk_usage(here).free < 10 * 1024**3:
        raise RuntimeError('Insufficient disk headroom')
    resource.setrlimit(resource.RLIMIT_FSIZE, (6*1024**3, 6*1024**3))
    sys.path.insert(0, str(kernels))
    from tensor_response import TensorResponse
    from validate_7269 import serial

    seed = here / 'seed'
    seed.mkdir(exist_ok=False)
    for name in ('INCAR', 'POSCAR'):
        shutil.copy(here/name, seed/name)
    (seed/'KPOINTS').write_text('Full 2D mesh\n0\nGamma\n6 6 1\n0 0 0\n')
    lib = Path('/work/sunq/mxy/bin/potpaw_PBE.64')
    potential = b''.join((lib/name/'POTCAR').read_bytes() for name in contract['species_order'])
    (seed/'POTCAR').write_bytes(potential)
    titles = re.findall(r'TITEL\s*=\s*([^\n]+)', potential.decode('latin-1'))
    if [title.split()[1] for title in titles] != contract['species_order']:
        raise ValueError('PAW species order changed')
    populations = [float(x) for x in re.findall(r'ZVAL\s*=\s*([0-9.]+)',potential.decode('latin-1'))]
    nocc = int(np.dot([4, 2], populations))
    record = {'job': os.environ['PBS_JOBID'], 'node': node, 'stage':'SCF', 'complete':False,
              'potcar_titles': titles, 'potcar_sha256':hashlib.sha256(potential).hexdigest(), 'electrons': nocc}
    def save():
        tmp = here/'result.tmp'
        tmp.write_text(json.dumps(record, default=serial, allow_nan=False, indent=2)+'\n')
        tmp.replace(here/'result.json')
    def run(command, name, cwd):
        with (cwd/name).open('w') as log:
            subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, check=True)
    save()
    start = time.monotonic()
    try:
        run(['mpirun','-np','4','-ppn','4','-hostfile',str(nodefile),
             '/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'], 'vasp.log', seed)
        out = (seed/'OUTCAR').read_text(errors='replace')
        if 'General timing' not in out or 'aborting loop because EDIFF is reached' not in out:
            raise ValueError('Electronic convergence failed')
        forces = read(seed/'OUTCAR', index=-1).get_forces(apply_constraint=False)
        record['scf_seconds'] = time.monotonic()-start
        record['max_force_eV_A'] = float(np.max(np.linalg.norm(forces,axis=1)))
        if not np.isfinite(forces).all() or record['max_force_eV_A'] > 0.10:
            record['stage'] = 'force_limit_exceeded'
            save()
            return
        record['stage'] = 'tensor'; save()
        cache = here/'tensor_shell_cache.npz'
        audit = here/'gauge_audit.json'
        run([sys.executable,str(kernels/'export_tensor.py'),str(seed),str(here),
             '--layer-mode','bilayer','--mesh','6','6','--cut','.5'], 'export.log', here)
        run([sys.executable,str(kernels/'build_stencil.py'),str(seed),str(here/'tensor_cache.npz'),str(cache)], 'stencil.log',here)
        run([sys.executable,str(kernels/'audit_gauge.py'),str(seed),str(cache),'--mesh','6','6','--output',str(audit)],'gauge.log',here)
        record['gauge_passed'] = json.loads(audit.read_text())['passed']
        with np.load(cache) as data:
            if int(data['nocc']) != nocc or not record['gauge_passed']:
                raise ValueError('Occupation or overlap mismatch')
            record['responses'] = {str(nb):{key:TensorResponse(data,nb).response(b=b)
                                   for key,b in [('standard',.002),('half_b',.001)]} for nb in (96,112)}
        record.update(stage='zero_field_requires_scientific_review',complete=True)
    except Exception as exc:
        record.update(stage='error', reason=str(exc)); raise
    finally:
        record['elapsed_seconds'] = time.monotonic()-start
        save()


if __name__ == '__main__':
    main()
