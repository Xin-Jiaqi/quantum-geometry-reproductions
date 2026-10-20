#!/usr/bin/env python3
"""Offline model regression; does NOT claim full-paper or SI physical validation."""
from __future__ import annotations
import hashlib,json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from igmc_lab.core import (dirac_response,rashba_analytic,rashba_integrate,
                           energy_from_frequency,physical_response_si)
from igmc_lab.independent_checks import (finite_field_wilson_derivative,
                                          orbital_liouville_kernel)

def main():
    manifest=json.loads((ROOT/'results/RELEASE_MANIFEST.json').read_text())
    for name,sha in manifest['files'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=sha:
            raise AssertionError(f'Release integrity mismatch: {name}')
    with (ROOT/'results/frozen/model_reference.json').open() as f:g=json.load(f)
    observed={
      'rashba_mu_minus0015':rashba_integrate(-.015,nrad=1200,nphi=96),
      'rashba_mu_plus0009':rashba_integrate(.009,nrad=1200,nphi=96),
      'dirac_gap_c_101':dirac_response(101,0,alpha='x',beta='y')['sigma_C'],
      'dirac_gap_dperomega_101':dirac_response(101,0,alpha='x',beta='x')['sigma_D_per_omega'],
      'dirac_metal_c_201':dirac_response(201,.85,alpha='x',beta='y')['sigma_C'],
      'dirac_metal_dperomega_201':dirac_response(201,.85,alpha='x',beta='x')['sigma_D_per_omega'],
    }
    for k,v in observed.items():
        expected=g['baseline'][k]
        if not math.isclose(v,expected,rel_tol=1e-9,abs_tol=1e-11):
            raise AssertionError(f'Frozen result drift {k}: observed={v}, expected={expected}')
    # Separate frozen SI/TR-regression checkpoint. Do not modify the old six baselines.
    h=json.loads((ROOT/'results/frozen/round2_si_reference.json').read_text())['expected']
    w=energy_from_frequency(1e12,protocol='cyclic_hz')
    new={
      'energy_1THz_eV':w,
      'si_factor_g2':physical_response_si(1.),
      'rashba_platform_si':physical_response_si(rashba_analytic(0.)),
      'dirac_insulator_D_si':physical_response_si(dirac_response(201,0.,omega_eV=w,alpha='x',beta='x')['sigma_D']),
      'dirac_insulator_D_fullKubo_si':physical_response_si(dirac_response(201,0.,omega_eV=w,alpha='x',beta='x')['Kubo_imag']),
      'dirac_metal_mu085_C_si':physical_response_si(dirac_response(201,.85,omega_eV=w,alpha='x',beta='y')['sigma_C']),
    }
    for k,v in new.items():
        if not math.isclose(v,h[k],rel_tol=1e-10,abs_tol=1e-11):
            raise AssertionError(f'Round2 SI frozen drift {k}: {v} vs {h[k]}')
    # Independent round-3 freeze: no changing old round-1 or round-2 standards.
    f3=json.loads((ROOT/'results/frozen/round3_independent_reference.json').read_text())['values']
    original=dirac_response(301,0.,omega_eV=w,alpha='x',beta='x')
    original_101=dirac_response(101,0.,omega_eV=w,alpha='x',beta='x')
    independent=orbital_liouville_kernel(101,omega_eV=w)
    latest={
       'omega_eV':w,
       'N301_wilson_b1e-4':finite_field_wilson_derivative(301,b_eV=1e-4),
       'N301_main_D_peromega':original['sigma_D_per_omega'],
       'N101_liouville_imag':independent.imag,
       'N101_kubo_imag':original_101['Kubo_imag'],
       'in_gap_D_SI':new['dirac_insulator_D_si'],
       'SM_literal_intermediate_D_SI_diagnostic':2*new['dirac_insulator_D_si']
    }
    for k,v in latest.items():
        if not math.isclose(v,f3[k],rel_tol=1e-9,abs_tol=1e-11):
            raise AssertionError(f'Round3 independent frozen drift {k}: {v} vs {f3[k]}')
    if abs(latest['N301_wilson_b1e-4']/latest['N301_main_D_peromega']-1) >= 8e-5:
        raise AssertionError('Independent finite-field Wilson error exceeds acceptance')
    if abs(latest['N101_liouville_imag']-latest['N101_kubo_imag']) >= 5e-11:
        raise AssertionError('Direct Liouville does not match full spectral Kubo')
    # Round 4 symmetry calculation: verify all 122 group closures + 244 tensors
    # against an independently generated, separately frozen all-MPG reference.
    rc4=subprocess.run([sys.executable,'scripts/audit_round4.py'],cwd=ROOT,
                      env=__import__('os').environ|{'PYTHONPATH':str(ROOT/'src')},
                      text=True,capture_output=True)
    if rc4.returncode:
        raise AssertionError('Round4 MPG frozen regression failed: '+rc4.stderr+'\n'+rc4.stdout)
    # Round 5 source transcription: 244 paper-table vs independent-MPG tensor
    # subspaces. One explicitly documented printed-parametric relation is
    # missing in SM Table II; fail on any *additional* mismatch.
    r5=subprocess.run([sys.executable,'scripts/audit_round5.py','--no-write'],cwd=ROOT,
                      env=__import__('os').environ|{'PYTHONPATH':str(ROOT/'src')},
                      text=True,capture_output=True)
    if r5.returncode:
        raise AssertionError('Round5 SM source table audit failed: '+r5.stderr+'\n'+r5.stdout)
    frozen5=json.loads((ROOT/'results/frozen/round5_source_table_reference_v1.json').read_text())
    audit5=json.loads((ROOT/'results/raw/round5_source_table_audit.json').read_text())
    for key,value in frozen5['expected'].items():
        if key=='exception':
            if any(audit5[key].get(k)!=v for k,v in value.items()):
                raise AssertionError('Round5 frozen exception mismatch')
        elif audit5[key]!=value:
            raise AssertionError(f'Round5 frozen audit mismatch {key}')
    for relative,sha in frozen5['sha256'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()!=sha:
            raise AssertionError(f'Round5 frozen SHA mismatch: {relative}')
    for relative,sha in frozen5['prior_round_frozen_sha256'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()!=sha:
            raise AssertionError(f'Prior frozen reference changed: {relative}')
    # Round 6: independent eigenvector-free d-vector formula and Gauss-Legendre
    # quadrature. Recompute the evidence without touching previously frozen CSVs.
    r6=subprocess.run([sys.executable,'scripts/audit_round6.py','--no-write'],cwd=ROOT,
                      env=__import__('os').environ|{'PYTHONPATH':str(ROOT/'src')},
                      text=True,capture_output=True)
    if r6.returncode:
        raise AssertionError('Round6 d-vector physics check failed: '+r6.stderr+'\n'+r6.stdout)
    observed6=json.loads(r6.stdout)
    ref6=json.loads((ROOT/'results/frozen/round6_dvector_reference_v1.json').read_text())
    for relative,sha in ref6['older_file_hashes_at_round5'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()!=sha:
            raise AssertionError(f'Older frozen checkpoint modified during Round6: {relative}')
    for relative,sha in ref6['new_evidence_file_hashes'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()!=sha:
            raise AssertionError(f'Round6 evidence SHA mismatch: {relative}')
    def _close_expected(a,b,key):
        if isinstance(a,dict):
            if not isinstance(b,dict) or set(a)!=set(b):raise AssertionError(f'Round6 frozen object mismatch {key}')
            for k in a:_close_expected(a[k],b[k],f'{key}.{k}')
        elif isinstance(a,list):
            if not isinstance(b,list) or len(a)!=len(b):raise AssertionError(f'Round6 frozen list mismatch {key}')
            for i,(aa,bb) in enumerate(zip(a,b)):_close_expected(aa,bb,f'{key}.{i}')
        elif isinstance(a,(float,int)) and not isinstance(a,bool):
            if not math.isclose(float(a),float(b),abs_tol=1e-12,rel_tol=1e-9):
                raise AssertionError(f'Round6 frozen value drift {key}: {a} != {b}')
        elif a!=b:raise AssertionError(f'Round6 frozen metadata mismatch {key}')
    for k,v in ref6['expected'].items():_close_expected(v,observed6[k],k)
    if observed6['max_complex_TR_conjugation_error']>1e-12:
        raise AssertionError('Round6 complex TR parity regression')
    if observed6['gauss_convergence_abs_Dperomega_last_step']>2e-8:
        raise AssertionError('Round6 independent BZ integration regression')

    rc7=subprocess.run([sys.executable,'scripts/verify_round7.py'],cwd=ROOT,
         env=__import__('os').environ|{'PYTHONPATH':str(ROOT/'src')},text=True,capture_output=True)
    if rc7.returncode:
        raise AssertionError('Round7 metal acceptance failed: '+rc7.stderr+'\n'+rc7.stdout)
    rc=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests','-v'],cwd=ROOT,env=__import__('os').environ|{'PYTHONPATH':str(ROOT/'src')})
    if rc.returncode:raise SystemExit(rc.returncode)
    print(json.dumps({'status':'PASS_MODEL_OFFLINE','scientific_status':'PARTIAL: Eq14 printed parity and Fig2(d) author graph unresolved; SM Eq34 intermediate double-count evidence documented','round1_frozen_cases':len(observed),'round2_si_frozen_cases':len(new),'round3_independent_frozen_cases':len(latest),'round4_mpg_groups':122,'round4_mpg_tensors':244,'round5_SM_tensor_shapes_compared':244,'round5_exact_source_matches':243,'round5_table_exception':"D -42'm'",'round6_eigenvector_free_gauss':True,'round6_full_complex_TR_cases':observed6['TR_exact_complex_checks'],'round7_metal_fermi_contour':True,'round7_publication':'local checkpoint only','test_suite':'see output'},indent=2))

if __name__=='__main__':main()
