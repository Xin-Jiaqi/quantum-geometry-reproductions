"""Independent 8x8 SCF check for A0, initialized with the converged 6x6 density.

Triggered by the unexpectedly large, nearly compensated layer response.
No gate or geometry change. Run only after A0 completes, in a PBS allocation.
"""
from pathlib import Path
import os, socket, json, re, shutil, subprocess, sys, time, hashlib
import numpy as np
from analyze import Response
from kresolved import c_k

if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super', 'cms-ssh')):
    raise SystemExit('PBS compute allocation required')
root = Path.cwd(); base = root/'A0_mesh8'; seed = base/'A'/'seed'
assert not (seed/'OUTCAR').exists(), 'Do not overwrite an existing mesh-check run'
source = root/'A0'/'A'/'seed'
assert json.loads((root/'A0'/'result.json').read_text())['complete']
seed.mkdir(parents=True, exist_ok=True)
for name in ('POSCAR', 'POTCAR', 'CHGCAR'):
    shutil.copyfile(source/name, seed/name)
text = (source/'INCAR').read_text().replace('ICHARG=2', 'ICHARG=1')+'\nLMAXMIX=4\n'
assert 'ICHARG=1' in text and 'ISTART=0' in text
(seed/'INCAR').write_text(text)
(seed/'KPOINTS').write_text('Independent Gamma mesh check\n0\nGamma\n8 8 1\n0 0 0\n')
record = {'name': 'A0_mesh8', 'label': 'A', 'F_external_V_nm': 0., 'complete': False,
          'job': os.environ['PBS_JOBID'], 'node': socket.gethostname(), 'scope': __doc__,
          'input_sha256': {n: hashlib.sha256((seed/n).read_bytes()).hexdigest() for n in ('POSCAR', 'INCAR', 'KPOINTS')},
          'source_density_sha256': hashlib.sha256((seed/'CHGCAR').read_bytes()).hexdigest()}
(base/'progress.json').write_text(json.dumps(record, indent=2))
ranks = len(Path(os.environ['PBS_NODEFILE']).read_text().splitlines()); assert ranks == 4
start = time.monotonic()
with (seed/'vasp.log').open('w') as f:
    subprocess.run(['mpirun', '-np', str(ranks), '-ppn', str(ranks), '-hostfile', os.environ['PBS_NODEFILE'],
                    '/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'], cwd=seed, stdout=f, stderr=subprocess.STDOUT, check=True)
out = (seed/'OUTCAR').read_text(errors='replace'); log = (seed/'vasp.log').read_text(errors='replace')
assert 'General timing' in out and 'aborting loop because EDIFF is reached' in out
assert not re.search(r'^\s*(DAV|RMM):\s*240\s', log, re.M)
record['vasp_seconds'] = time.monotonic()-start
subprocess.run([sys.executable, str(root/'export.py')], cwd=base, env=dict(os.environ, STRUCTURE='A'), check=True)
d = np.load(base/'A'/'cache.npz')
record['export'] = json.loads((base/'A'/'export.json').read_text())
record['response'] = {str(nb): {'standard': Response(d, nb).response(), 'half_b': Response(d, nb).response(b=.001)} for nb in (128, 144, 160)}
rr=Response(d,160);ck=c_k(rr,rr.h)
assert max(abs(ck.mean(axis=0)-record['response']['160']['standard']['C_total_lower_upper']))<1e-7
record['k']=d['k'].tolist();record['cell']=d['cell'].tolist();record['C_k_total_lower_upper']=ck.tolist()
record['complete'] = True
(base/'result.json').write_text(json.dumps(record, indent=2))
print('MESH_CHECK_COMPLETE', flush=True)
