"""Apply the same six-link kernel to the two published project gate seeds."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys

if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super', 'cms-ssh')):
    raise SystemExit('PBS compute node required')
here = Path.cwd()
root = here.parent
spec = json.loads((here/'inputs.json').read_text())
for name, expected in spec['sources'].items():
    if hashlib.sha256((root/name).read_bytes()).hexdigest() != expected:
        raise SystemExit('Kernel provenance differs: '+name)
sys.path.insert(0, str(root))
sys.path.insert(0, str(root/'overnight7268'))
from night_checks import response_acceptance
from tensor_response import TensorResponse
from validate_7269 import serial
from ase.io import read
import numpy as np

archive = Path('/work/sunq/xjq/mosse-gate-scf-20261109')
for case, field in [('fm010', -.10), ('fp010', .10)]:
    seed = archive/case/'A/seed'
    for name, expected in spec[case].items():
        if hashlib.sha256((seed/name).read_bytes()).hexdigest() != expected:
            raise SystemExit('Historical field input differs: '+case+'/'+name)
    if shutil.disk_usage(here).free < 10*1024**3:
        raise SystemExit('Insufficient filesystem space')
    dest = here/case
    dest.mkdir(exist_ok=False)
    def run(script, *args):
        subprocess.run([sys.executable, str(root/script), *map(str,args)], check=True, timeout=240)
    run('export_tensor.py', seed, dest, '--layer-mode', 'bilayer', '--mesh', 6, 6, '--cut', .5)
    cache = dest/'tensor_shell_cache.npz'
    run('build_stencil.py', seed, dest/'tensor_cache.npz', cache)
    run('audit_gauge.py', seed, cache, '--mesh', 6, 6, '--output', dest/'gauge.json')
    record = {'case':case, 'label':'A', 'Fz_V_nm':field, 'mesh':6, 'job':os.environ['PBS_JOBID'],
              'node':socket.gethostname(), 'source':str(seed), 'input_sha256':spec[case],
              'max_force_eV_A':float(np.linalg.norm(read(seed/'OUTCAR').get_forces(apply_constraint=False),axis=1).max()),
              'gauge_audit':json.loads((dest/'gauge.json').read_text()),
              'export':json.loads((dest/'export.json').read_text()), 'responses':{}}
    with np.load(cache) as f:
        for nb in (96,112):
            model = TensorResponse(f,nb)
            record['responses'][str(nb)] = {'standard':model.response(b=.002),'half_b':model.response(b=.001)}
    record['acceptance'] = response_acceptance(record,.1)
    (dest/'result.json').write_text(json.dumps(record,default=serial,indent=2,allow_nan=False)+'\n')
    if not record['acceptance']['passed']:
        raise SystemExit('Field reference did not pass: '+case)
