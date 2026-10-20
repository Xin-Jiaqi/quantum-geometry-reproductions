#!/usr/bin/env python3
"""Strict checks for new Round7 metal evidence, no mutation of frozen products."""
import hashlib,json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def digest(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def check_numbers(exp,obs,path='root'):
    if isinstance(exp,dict):
        if not isinstance(obs,dict) or set(exp)!=set(obs):raise AssertionError('keys mismatch '+path)
        for k,v in exp.items():check_numbers(v,obs[k],path+'.'+k)
    elif isinstance(exp,list):
        if not isinstance(obs,list) or len(exp)!=len(obs):raise AssertionError('list mismatch '+path)
        for i,(a,b) in enumerate(zip(exp,obs)):check_numbers(a,b,path+f'[{i}]')
    elif isinstance(exp,(float,int)) and not isinstance(exp,bool):
        if not math.isclose(float(exp),float(obs),rel_tol=2e-9,abs_tol=1e-12):
            raise AssertionError(f'number mismatch {path}: {exp} vs {obs}')
    elif exp!=obs:raise AssertionError(f'metadata mismatch {path}: {exp} vs {obs}')

def main():
    ref=json.loads((ROOT/'results/frozen/round7_metal_reference_v1.json').read_text())
    for path,sha in ref['previous_frozen_files'].items():
        if digest(path)!=sha:raise AssertionError('changed earlier frozen '+path)
    for path,sha in ref['round7_evidence_files'].items():
        if digest(path)!=sha:raise AssertionError('changed round7 evidence '+path)
    for path,sha in ref['round7_source_files'].items():
        if digest(path)!=sha:raise AssertionError('changed round7 source '+path)
    r=subprocess.run([sys.executable,'scripts/audit_round7.py','--no-write'],cwd=ROOT,
         text=True,capture_output=True,env=__import__('os').environ|{'PYTHONPATH':str(ROOT/'src')})
    if r.returncode:raise AssertionError(r.stderr+'\n'+r.stdout)
    data=json.loads(r.stdout)
    check_numbers(ref['expected'],{k:data[k] for k in ref['expected']})
    if data['max_absolute_convergence_C_Ainv']>4e-8 or data['max_absolute_convergence_Dperomega_Ainv_eVinv']>2e-8:
        raise AssertionError('insufficient metal resolution')
    if data['max_absolute_occupied_vs_FS_eq13_C_Ainv']>2e-8:
        raise AssertionError('Eq13 Fermi line disagrees with sea')
    print('PASS_ROUND7_MODEL_METAL: metal Fermi contours, independent Eq13, frozen checksums; whole paper PARTIAL')
if __name__=='__main__':main()
