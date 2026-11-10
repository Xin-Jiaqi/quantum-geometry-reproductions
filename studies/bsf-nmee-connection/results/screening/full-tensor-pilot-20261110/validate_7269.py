"""Bounded 7269 cache benchmark; accept/reject before new 7268 DFT."""
from pathlib import Path
import importlib.util
import json
import os
import socket
import subprocess
import sys
import time
import numpy as np
from tensor_response import TensorResponse


def serial(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    raise TypeError(type(x).__name__)


def main():
    if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
        raise SystemExit('PBS compute only')
    root = Path.cwd()
    source = Path('/work/sunq/xjq/mosse-strain-20261108/e000/A')
    output = root/'anchor7269'
    start = time.monotonic()
    subprocess.run([sys.executable, str(root/'export_tensor.py'), str(source/'seed'), str(output), '--layer-mode','bilayer','--mesh','6','6','--cut','0.5'], check=True)
    spec = importlib.util.spec_from_file_location('legacy', root/'legacy_response.py')
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    old = np.load(source/'cache.npz')
    new = np.load(output/'tensor_cache.npz')
    results = []
    for nb in (96,112):
        t = time.monotonic()
        r = TensorResponse(new,nb)
        standard = r.response(b=.002)
        half = r.response(b=.001)
        previous = legacy.Response(old,nb).response(b=.002)
        c = np.asarray(standard['C'])
        ch = np.asarray(half['C'])
        reference = np.asarray(previous['C_total_lower_upper'])
        legacy_error = float(np.max(abs(c[:,0,0]-reference)))
        scale = max(float(np.max(abs(c[1:,0,0]))),1e-8)
        forbidden = float(np.max(abs(c[:,[0,1,2,2],[1,0,0,1]])))
        isotropy = float(np.max(abs(c[:,0,0]-c[:,1,1])))
        probe = float(np.max(abs(c-ch)))
        passed = bool(standard['gap_eV'] > .05 and half['gap_eV'] > .05 and legacy_error <= 1e-7*max(1,float(np.max(abs(reference)))) and
                      forbidden <= .02*scale and isotropy <= .05*scale and probe <= .01*scale)
        results.append(dict(nb=nb, standard=standard,half_b=half,legacy=previous,
                            legacy_max_abs_error=legacy_error, forbidden_max_abs=forbidden,
                            Cxx_Cyy_max_abs=isotropy, probe_max_abs=probe, layer_scale=scale,
                            response_and_checks_seconds=time.monotonic()-t,passed=passed))
    rec=dict(passed=all(r['passed'] for r in results),results=results,
             thresholds='Legacy relative1e-7; forbidden2%, Cxx-Cyy5%, probe1% of max layer Cxx; independent PAW layer sum checked at export. These are pilot acceptance, not final convergence.',
             elapsed_seconds=time.monotonic()-start,node=socket.gethostname(),job=os.environ['PBS_JOBID'])
    (root/'validation7269.json').write_text(json.dumps(rec,default=serial,indent=2))
    print('VALIDATION_PASSED',rec['passed'],flush=True)
    if not rec['passed']:raise SystemExit(2)


if __name__=='__main__':main()
