#!/usr/bin/env python3
"""Fail-closed Round8 scientific kernel and preservation gate."""
from __future__ import annotations
import hashlib,json,math,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def equivalent(a,b,tag='root'):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b):raise AssertionError('different metadata '+tag)
        for k,v in a.items():equivalent(v,b[k],tag+'.'+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise AssertionError('list mismatch '+tag)
        for i,(aa,bb) in enumerate(zip(a,b)):equivalent(aa,bb,tag+f'.{i}')
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        if not math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-12):raise AssertionError(f'value mismatch {tag}: {a} vs {b}')
    elif a!=b:raise AssertionError('metadata mismatch '+tag)

def main():
    ref=json.loads((ROOT/'results/frozen/round8_green_wavepacket_reference_v1.json').read_text())
    if sha('results/RELEASE_MANIFEST.json')!=ref['prior_release_manifest_SHA256']:
        raise AssertionError('Round7 release manifest was changed')
    for name,expect in ref['old_manifest_mutables_preserved'].items():
        archive=ROOT/'archive/round7_integrity_files'/name
        if hashlib.sha256(archive.read_bytes()).hexdigest()!=expect:
            raise AssertionError('Archived Round7 living file changed '+name)
    for name,expect in ref['prior_frozen'].items():
        if sha(name)!=expect:raise AssertionError('Prior frozen file changed '+name)
    for group in ('source_files','evidence_files'):
        for name,expect in ref[group].items():
            if sha(name)!=expect:raise AssertionError(f'Round8 {group} changed {name}')
    run=subprocess.run([sys.executable,'scripts/audit_round8.py','--no-write'],cwd=ROOT,
          env=os.environ|{'PYTHONPATH':str(ROOT/'src')},text=True,capture_output=True)
    if run.returncode:raise AssertionError(run.stdout+'\n'+run.stderr)
    observed=json.loads(run.stdout)
    equivalent(ref['expected'],observed)
    if observed['max_per_k_GF_vs_wave_abs']>=3e-12:
        raise AssertionError('Green/wavepacket independent complex check failed')
    if abs(observed['SM_eq46_FT_coefficient_positive_harmonic']-.5)>1e-12:
        raise AssertionError('Fourier cosine harmonic mismatch')
    if abs(observed['SM_eq34_metric_intermediate_over_main']-2)>1e-12:
        raise AssertionError('SM diagnostic lost provenance')
    print('PASS_ROUND8_MODEL_SEMICLASSICAL_LESSER_GF: scientific scope PARTIAL; Fig2(d) author curve not closed')
if __name__=='__main__':main()
