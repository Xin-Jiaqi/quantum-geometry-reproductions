"""A single allocation, serial states, explicit quality gates; never submits jobs."""
from pathlib import Path
import hashlib
import json
import os
import re
import resource
import shutil
import signal
import socket
import subprocess
import sys
import time


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_tags(text, tags):
    for key, value in tags.items():
        pattern = rf'^\s*{key}\s*=.*$'
        text = re.sub(pattern, f'{key}={value}', text, flags=re.M) if re.search(pattern, text, re.M) else text+f'\n{key}={value}\n'
    return text


def main():
    nodefile = os.environ.get('PBS_NODEFILE', '')
    if not os.getenv('PBS_JOBID') or not nodefile or socket.gethostname().startswith(('super', 'cms-ssh')):
        raise SystemExit('PBS compute allocation required')
    if len(Path(nodefile).read_text().splitlines()) != 4:
        raise SystemExit('Exactly four MPI ranks required')
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_FSIZE, (6*1024**3, 6*1024**3))
    here = Path.cwd()
    root = here.parent
    contract = json.loads((here/'contract.json').read_text())
    if contract['authorization'] != '20261111_progress_single_registry_tauB':
        raise SystemExit('Missing scoped authorization')
    for name, expected in contract['source_hashes'].items():
        if sha(root/name) != expected:
            raise SystemExit('Source changed: '+name)
    for name, expected in contract['normalized_text_source_hashes'].items():
        actual = hashlib.sha256((root/name).read_text().strip().encode()).hexdigest()
        if actual != expected:
            raise SystemExit('PAW source changed: '+name)
    for name, expected in contract['runner_hashes'].items():
        if sha(here/name) != expected:
            raise SystemExit('Runner changed: '+name)
    sys.path.insert(0, str(root))
    sys.path.insert(0, str(here))
    import numpy as np
    from ase.io import read
    from tensor_response import TensorResponse
    from validate_7269 import serial
    from night_checks import response_acceptance
    data = here/'campaign'
    data.mkdir(exist_ok=False)
    start = time.monotonic()
    deadline = start+3.9*3600
    summary = {'job': os.environ['PBS_JOBID'], 'node': socket.gethostname(),
               'complete': False, 'status': 'input_verification', 'completed_states': [],
               'only_material': 'MoSSe S/Se E,tau_B,d (BiDB7266)',
               'coordinate_convention': 'x along a1, a2 at120degrees, +z lower-to-upper; Fz positive physical +z',
               'allocated_CPUs': 4, 'max_wall_hours': 4, 'automatic_retries': 0,
               'new_jobs_submitted_from_worker': 0}

    def report():
        summary['elapsed_seconds'] = time.monotonic()-start
        save(data/'summary.json', summary)

    def run(command, cwd, logname, seconds):
        available = deadline-time.monotonic()-180
        if available <= 0:
            raise RuntimeError('Global time budget exhausted')
        with (cwd/logname).open('w') as log:
            proc = subprocess.Popen(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = proc.wait(timeout=min(seconds, available))
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                raise RuntimeError('Command timed out; allocation stops, no retry: '+str(command[0]))
        files = [f for f in data.rglob('*') if f.is_file()]
        if any(f.stat().st_size > 6*1024**3 for f in files):
            raise RuntimeError('Individual output exceeded6GiB; stop without retry')
        if sum(f.stat().st_size for f in files) > 10*1024**3:
            raise RuntimeError('Campaign exceeded10GiB; stop without retry')
        if code:
            raise RuntimeError(f'Command failed ({code}); no retry: {command[0]}')

    def space_and_time_gate():
        if deadline-time.monotonic() < 85*60:
            raise RuntimeError('Less than85min reserved for next state; stop without submitting or retrying')
        usage = sum(f.stat().st_size for f in data.rglob('*') if f.is_file())
        if usage + 2*1024**3 > 10*1024**3 or shutil.disk_usage(data).free < 10*1024**3:
            raise RuntimeError('Campaign disk budget/free-space guard; no cleanup or retry')

    def state(mesh, label, suffix, field):
        space_and_time_gate()
        case = f'k{mesh}_{label}{suffix}'
        summary.update(status='running_state', current_state=case)
        report()
        base = data/case
        seed = base/'seed'
        seed.mkdir(parents=True, exist_ok=False)
        shutil.copyfile(here/f'POSCAR.{label}', seed/'POSCAR')
        template = root/'pilot7268'/('A0' if label == 'A' else 'lower0')/'seed'
        shutil.copyfile(template/'POTCAR', seed/'POTCAR')
        (seed/'KPOINTS').write_text(f'Full Gamma mesh for Layer NMEE\n0\nGamma\n{mesh} {mesh} 1\n0 0 0\n')
        restart = data/f'k{mesh}_{label}0'/'seed' if field else None
        if field and label == 'lower':
            restart = root/'pilot7268/lower0/seed'
        if restart:
            for name in ('WAVECAR', 'CHGCAR'):
                shutil.copyfile(restart/name, seed/name)
        text = replace_tags((template/'INCAR').read_text(), {
            'SYSTEM': f'7266 fixed BiDB geometry {case} Fz={field:+.3f} V/nm',
            'ISTART': 1 if restart else 0, 'ICHARG': 1 if restart else 2,
            'EFIELD': f'{-field/10:.8f}', 'NSW': 0, 'IBRION': -1,
            'NBANDS': 112 if label == 'A' else 64, 'KPAR': 1, 'NCORE': 1,
            'LDIPOL': 'T', 'IDIPOL': 3, 'LWANNIER90': 'T', 'LWRITE_MMN_AMN': 'T', 'DIPOL': '0.5 0.5 0.5'})
        for required in ('LSORBIT=T', 'BEXT=0 0 0', 'EDIFF=1e-9', 'ENCUT=450', 'IVDW=11', 'ISYM=-1', 'NELM=240', 'LDIPOL=T', 'IDIPOL=3'):
            if required not in text:
                raise RuntimeError('Electronic protocol changed: '+required)
        (seed/'INCAR').write_text(text)
        record = {'case': case, 'mesh': mesh, 'label': label, 'Fz_V_nm': field,
                  'fixed_ions': True, 'complete': False, 'stage': 'SCF',
                  'input_sha256': {name: sha(seed/name) for name in ('POSCAR', 'INCAR', 'KPOINTS', 'POTCAR')},
                  'job': os.environ['PBS_JOBID'], 'node': socket.gethostname(),
                  'restart_input_hashes': {name: sha(seed/name) for name in ('WAVECAR', 'CHGCAR')} if restart else {}}
        save(base/'result.json', record)
        before = time.monotonic()
        run(['mpirun', '-np', '4', '-ppn', '4', '-hostfile', nodefile,
             '/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'], seed, 'vasp.log', 55*60)
        record['scf_seconds'] = time.monotonic()-before
        output = (seed/'OUTCAR').read_text(errors='replace')
        log = (seed/'vasp.log').read_text(errors='replace')
        if ('General timing' not in output or 'aborting loop because EDIFF is reached' not in output
                or re.search(r'^\s*(DAV|RMM):\s*240\s', log, re.M)):
            raise RuntimeError('SCF convergence gate failed')
        record['electronic_iterations'] = len(re.findall(r'^\s*(DAV|RMM):', log, re.M))
        record['max_force_eV_A'] = float(np.max(np.linalg.norm(read(seed/'OUTCAR', index=-1).get_forces(apply_constraint=False), axis=1)))
        force_limit = .1
        if not np.isfinite(record['max_force_eV_A']) or record['max_force_eV_A'] > force_limit:
            record.update(stage='force_gate_failed', force_limit_eV_A=force_limit)
            save(base/'result.json', record)
            raise RuntimeError('Force gate failed before response or later fields')
        record['stage'] = 'export'
        save(base/'result.json', record)
        run([sys.executable, str(root/'export_tensor.py'), str(seed), str(base),
             '--layer-mode', 'bilayer' if label == 'A' else 'mono', '--mesh', str(mesh), str(mesh), '--cut', '.5'], base, 'export.log', 15*60)
        record['export'] = json.loads((base/'export.json').read_text())
        run([sys.executable, str(root/'build_stencil.py'), str(seed), str(base/'tensor_cache.npz'), str(base/'tensor_shell_cache.npz')], base, 'stencil.log', 5*60)
        run([sys.executable, str(root/'audit_gauge.py'), str(seed), str(base/'tensor_shell_cache.npz'), '--mesh', str(mesh), str(mesh), '--output', str(base/'gauge_audit.json')], base, 'gauge.log', 5*60)
        record['gauge_audit'] = json.loads((base/'gauge_audit.json').read_text())
        record['responses'] = {}
        with np.load(base/'tensor_shell_cache.npz') as cache:
            for nb in ((96, 112) if label == 'A' else (48, 56, 64)):
                response = TensorResponse(cache, nb)
                record['responses'][str(nb)] = {'standard': response.response(b=.002), 'half_b': response.response(b=.001)}
        # Convert NumPy values before atomic JSON writing.
        record = json.loads(json.dumps(record, default=serial, allow_nan=False))
        record['acceptance'] = response_acceptance(record, force_limit)
        record.update(complete=True, stage='accepted' if record['acceptance']['passed'] else 'response_gate_failed',
                      elapsed_seconds=time.monotonic()-before)
        save(base/'result.json', record)
        if not record['acceptance']['passed']:
            raise RuntimeError('Full response scientific gate failed: '+case)
        summary['completed_states'].append({'case': case, 'elapsed_seconds': record['elapsed_seconds'],
                                             'gap_eV': record['acceptance']['rows'][-1]['gap_eV'],
                                             'Cxx': record['acceptance']['rows'][-1]['Cxx']})
        report()
        return record

    report()
    try:
        if not json.loads((root/'final_anchor_acceptance.json').read_text())['passed']:
            raise RuntimeError('7269 method acceptance missing')
        mono = root/'pilot7268/lower0'
        record = json.loads((mono/'result.json').read_text())
        prior = json.dumps(record, sort_keys=True, separators=(',', ':')).encode()
        if hashlib.sha256(prior).hexdigest() != contract['prior_mono_canonical_json_sha256']:
            raise RuntimeError('Historical monolayer result changed')
        if sha(here/'POSCAR.lower') != sha(mono/'seed/POSCAR'):
            raise RuntimeError('Isolated geometry differs from archived zero reference')
        for name, expected in record['export']['input_sha256'].items():
            if sha(mono/'seed'/name) != expected:
                raise RuntimeError('Historical reference input differs: '+name)
        record['label'] = 'lower'
        record['gauge_audit'] = json.loads((mono/'gauge_audit.json').read_text())
        audited = response_acceptance(record, .1)
        if not record['complete'] or record['Fz_V_nm']!=0 or not audited['passed']:
            raise RuntimeError('Matched isolated zero reference failed')
        summary['reused_zero_reference'] = {'job':record['job'], 'acceptance':audited,
                                             'applies_to':['lower','upper']}
        report()
        double = {'0':state(6,'A','0',0.)}
        single = {'0':record}
        for label, table in [('A',double), ('lower',single)]:
            for suffix, field in [('m',-.1),('p',.1)]:
                table[suffix] = state(6,label,suffix,field)
        summary['Gxx'] = {}
        for label, table, bands in [('bilayer',double,(96,112)),('matched_isolate',single,(48,56,64))]:
            summary['Gxx'][label] = {}
            for nb in bands:
                plus = np.asarray(table['p']['responses'][str(nb)]['standard']['C'])
                minus = np.asarray(table['m']['responses'][str(nb)]['standard']['C'])
                scopes = 3 if label=='bilayer' else 1
                summary['Gxx'][label][str(nb)] = ((plus-minus)[:scopes,0,0]/.2).tolist()
        summary.update(status='six_by_six_complete_requires_review', complete=True)
        summary.pop('current_state', None)
        report()
    except Exception as error:
        summary.update(status='stopped_no_retry', stop_reason=str(error), complete=False)
        report()
        raise


if __name__ == '__main__':
    main()
