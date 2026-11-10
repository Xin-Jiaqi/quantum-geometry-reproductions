"""Predeclared four/six-link discrimination, preserving failed anchor evidence."""
from pathlib import Path
import json,os,socket,subprocess,sys,time
import numpy as np
from tensor_response import TensorResponse
from build_stencil import complete_stencil
from validate_7269 import serial


def main():
    if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute only')
    root=Path.cwd();seed=Path('/work/sunq/xjq/mosse-strain-20261108/e000/A/seed');start=time.monotonic()
    subprocess.run([sys.executable,str(root/'audit_gauge.py'),str(seed),str(root/'anchor7269/tensor_cache.npz'),'--mesh','6','6','--output',str(root/'gauge_audit.json')],check=True)
    stencil=complete_stencil(seed,root/'anchor7269/tensor_cache.npz',root/'anchor7269/tensor_shell_cache.npz')
    old=json.loads((root/'validation7269.json').read_text())
    rows=[]
    for nb in (96,112):
        t=time.monotonic();r=TensorResponse(np.load(root/'anchor7269/tensor_shell_cache.npz'),nb)
        standard=r.response(b=.002);half=r.response(b=.001)
        c=np.asarray(standard['C']);ch=np.asarray(half['C']);scale=max(float(np.max(abs(c[1:,0,0]))),1e-8)
        prev=next(x for x in old['results'] if x['nb']==nb)
        legacy=np.array(prev['legacy']['C_total_lower_upper'])
        four=np.array(prev['standard']['C'])
        forbidden=float(np.max(abs(c[:,[0,1,2,2],[1,0,0,1]])))
        isotropy=float(np.max(abs(c[:,0,0]-c[:,1,1])))
        probe=float(np.max(abs(c-ch)))
        cxx_change=float(np.max(abs(c[:,0,0]-legacy)))
        legacy_four_error=prev['legacy_max_abs_error']
        passed=bool(standard['gap_eV']>.05 and half['gap_eV']>.05 and legacy_four_error<1e-7*max(1,float(np.max(abs(legacy)))) and forbidden<.02*scale and isotropy<.05*scale and probe<.01*scale and cxx_change<.05*scale)
        rows.append(dict(nb=nb,standard=standard,half_b=half,passed=passed,forbidden_max_abs=forbidden,Cxx_Cyy_max_abs=isotropy,probe_max_abs=probe,Cxx_discretization_change=cxx_change,legacy_four_link_max_abs_error=legacy_four_error,layer_scale=scale,old_four_link_forbidden_max_abs=prev['forbidden_max_abs'],response_seconds=time.monotonic()-t))
    gauge=json.loads((root/'gauge_audit.json').read_text())
    rec={'passed':all(x['passed'] for x in rows),'rows':rows,'stencil':stencil,'gauge':gauge,'elapsed_seconds':time.monotonic()-start,'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'interpretation':'Full MMN shell changes derivative geometry, not spin tensor projection. Four-link failure remains archived. Legacy exact regression applies to unchanged four-link mode; corrected Cxx must remain within5%.'}
    (root/'stencil_diagnosis.json').write_text(json.dumps(rec,default=serial,indent=2))
    print('STENCIL_DIAGNOSIS_PASSED',rec['passed'],flush=True)
    if not rec['passed']:raise SystemExit(2)


if __name__=='__main__':main()
