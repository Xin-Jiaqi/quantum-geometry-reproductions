"""Final bounded acceptance after six-direction interior/boundary gauge audit."""
from pathlib import Path
import json,os,socket,subprocess,sys,hashlib
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute only')
root=Path.cwd();seed=Path('/work/sunq/xjq/mosse-strain-20261108/e000/A/seed')
subprocess.run([sys.executable,str(root/'audit_gauge.py'),str(seed),str(root/'anchor7269/tensor_shell_cache.npz'),'--mesh','6','6','--output',str(root/'gauge_shell_audit.json')],check=True)
d=json.loads((root/'stencil_diagnosis.json').read_text());g=json.loads((root/'gauge_shell_audit.json').read_text())
result={'passed':bool(d['passed'] and g['passed']),'node':socket.gethostname(),'job':os.environ['PBS_JOBID'],'scope':'7269 full shell tensor pilot acceptance, not mesh-converged material amplitude','evidence_sha256':{f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in ['stencil_diagnosis.json','gauge_shell_audit.json']}}
(root/'final_anchor_acceptance.json').write_text(json.dumps(result,indent=2))
print('FINAL_ANCHOR_ACCEPTED',result['passed'],flush=True)
if not result['passed']:raise SystemExit(2)
