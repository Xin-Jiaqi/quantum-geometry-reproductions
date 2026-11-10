"""One explicitly authorized, bounded z-only structural check; no retry."""
from pathlib import Path
import json
import hashlib
import os
import re
import shutil
import socket
import subprocess
import time

import numpy as np
from ase.io import read, write
from ase.constraints import FixScaled


def main():
    nodefile = os.environ.get('PBS_NODEFILE', '')
    if not os.getenv('PBS_JOBID') or not nodefile or socket.gethostname().startswith(('super', 'cms-ssh')):
        raise SystemExit('PBS compute node required')
    if len(Path(nodefile).read_text().splitlines()) != 4:
        raise SystemExit('Exactly four ranks required')
    root = Path.cwd()
    contract = json.loads((root / 'structural_followup_proposal.json').read_text())
    if contract['status'] != 'user_authorized_20261110':
        raise SystemExit('Explicit structural authorization required')
    base = root / 'pilot7268' / 'A_zrelax'
    base.mkdir(exist_ok=False)
    seed = base / 'seed'
    seed.mkdir()
    source = root / 'pilot7268' / 'A0' / 'seed'
    initial = read(source / 'POSCAR')
    initial.set_constraint(FixScaled(range(len(initial)), mask=(True, True, False)))
    write(seed / 'POSCAR', initial, format='vasp', direct=True, sort=False)
    for name in ('POTCAR', 'KPOINTS', 'WAVECAR', 'CHGCAR'):
        shutil.copyfile(source / name, seed / name)
    incar = (source / 'INCAR').read_text()
    for key, value in {'SYSTEM': '7268 zero-field bounded z-only relaxation',
                       'ISTART': '1', 'ICHARG': '1', 'NSW': '12', 'IBRION': '2',
                       'ISIF': '2', 'POTIM': '0.3', 'EDIFFG': '-0.03',
                       'EFIELD': '0', 'LWANNIER90': '.FALSE.',
                       'LWRITE_MMN_AMN': '.FALSE.'}.items():
        pattern = rf'^\s*{key}\s*=.*$'
        incar = re.sub(pattern, f'{key}={value}', incar, flags=re.M) if re.search(pattern, incar, re.M) else incar + f'\n{key}={value}\n'
    for required in ('LSORBIT=T', 'BEXT=0 0 0', 'EDIFF=1e-9', 'ENCUT=450', 'IVDW=11'):
        if required not in incar:
            raise ValueError('Inherited electronic protocol mismatch: ' + required)
    (seed / 'INCAR').write_text(incar)
    record = {'job': os.environ['PBS_JOBID'], 'node': socket.gethostname(),
              'complete': False, 'passed': False, 'stage': 'z_relaxation',
              'max_ionic_steps': 12, 'force_target_eV_A': 0.03, 'gate_release': False}
    path = base / 'result.json'
    def save():
        path.write_text(json.dumps(record, indent=2) + '\n')
    save()
    start = time.monotonic()
    with (seed / 'vasp.log').open('w') as log:
        proc = subprocess.run(['mpirun', '-np', '4', '-ppn', '4', '-hostfile', nodefile,
                               '/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'],
                              cwd=seed, stdout=log, stderr=subprocess.STDOUT)
    record['elapsed_seconds'] = time.monotonic() - start
    record['returncode'] = proc.returncode
    save()
    if proc.returncode:
        raise SystemExit('VASP failure; no retry')
    output = (seed / 'OUTCAR').read_text(errors='replace')
    states = read(seed / 'OUTCAR', index=':')
    final = states[-1]
    forces = final.get_forces(apply_constraint=False)
    delta = final.positions - initial.positions
    delta[:, 2] -= delta[:, 2].mean()
    record.update(ionic_steps=len(states),
                  final_scaled_positions=final.get_scaled_positions().tolist(),
                  final_cell_A=final.cell.array.tolist(),
                  symbols=final.get_chemical_symbols(),
                  CONTCAR_sha256=hashlib.sha256((seed / 'CONTCAR').read_bytes()).hexdigest(),
                  max_force_eV_A=float(np.max(np.linalg.norm(forces, axis=1))),
                  forces_eV_A=forces.tolist(),
                  max_xy_displacement_A=float(np.max(np.abs(delta[:, :2]))),
                  max_z_displacement_after_translation_A=float(np.max(np.abs(delta[:, 2]))),
                  cell_max_change_A=float(np.max(np.abs(final.cell.array - initial.cell.array))),
                  ionic_converged='reached required accuracy - stopping structural energy minimisation' in output,
                  electronic_converged_steps=output.count('aborting loop because EDIFF is reached'),
                  force_history_eV_A=[float(np.max(np.linalg.norm(s.get_forces(apply_constraint=False), axis=1))) for s in states])
    record['passed'] = bool(record['ionic_converged'] and record['electronic_converged_steps'] >= len(states)
                            and record['max_force_eV_A'] < 0.03 and len(states) <= 12
                            and record['max_xy_displacement_A'] < 1e-6
                            and record['cell_max_change_A'] < 1e-8
                            and record['max_z_displacement_after_translation_A'] <= 0.15)
    record.update(complete=True, stage='structural_check_only_requires_zero_response_revalidation')
    save()
    if not record['passed']:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
