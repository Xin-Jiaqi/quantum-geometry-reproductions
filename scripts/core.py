#!/usr/bin/env python3
"""Run the independently verified paper recipes without merging their physical conventions."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / 'physics/capabilities.json'


def capabilities():
    data = json.loads(CATALOG.read_text())
    if data['schema_version'] != 1:
        raise ValueError('Unsupported capability catalog')
    rows = data['capabilities']
    if len({r['id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate capability ID')
    for row in rows:
        if row['status'] != 'MODEL_READY' or not row['limits'] or not row['units']:
            raise ValueError('Each capability needs a scope, units and limitations')
        directory = (ROOT / row['directory']).resolve()
        for path in (directory, directory / row['pythonpath'],
                     directory / row['example'], directory / row['documentation'], directory / 'tests'):
            if not path.resolve().is_relative_to(ROOT) or not path.exists():
                raise ValueError(f'Missing or outside-repository entry: {path}')
    return rows


def command(row, action):
    directory = ROOT / row['directory']
    env = os.environ.copy()
    # Each process imports its selected paper, without accidentally reusing another adapter.
    env['PYTHONPATH'] = str(directory / row['pythonpath'])
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        env[key] = '1'
    args = ([sys.executable, row['example']] if action == 'run' else
            [sys.executable, '-m', 'unittest', 'discover', '-s', 'tests', '-v'])
    return args, directory, env


def run_example(row):
    args, cwd, env = command(row, 'run')
    result = subprocess.run(args, cwd=cwd, env=env, text=True, capture_output=True, timeout=1200)
    if result.returncode:
        raise RuntimeError(f"{row['id']} recipe failed:\n{result.stderr}\n{result.stdout}")
    data = json.loads(result.stdout)
    if not isinstance(data, dict) or not data:
        raise ValueError(f"{row['id']} recipe must return a nonempty JSON object")
    commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True,
                            capture_output=True, check=True).stdout.strip()
    dirty = bool(subprocess.run(['git', 'status', '--porcelain'], cwd=ROOT, text=True,
                               capture_output=True, check=True).stdout)
    return dict(capability=row['id'], source_commit=commit, working_tree_dirty=dirty,
                units=row['units'], limits=row['limits'], result=data)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['list', 'check', 'run', 'test'])
    p.add_argument('capability', nargs='?', help='ID from list, or all')
    p.add_argument('--output', type=Path, help='Save run JSON; otherwise print it')
    a = p.parse_args(argv)
    rows = capabilities()
    if a.action in ('list', 'check'):
        if a.capability or a.output:
            p.error('list/check do not accept a capability or output')
        if a.action == 'list':
            for r in rows: print(f"{r['id']:18} {r['status']:12} {r['capability']}")
        else: print(f'OK: {len(rows)} capability entries resolve; this is not a physics test')
        return
    if a.capability not in [r['id'] for r in rows] + ['all']:
        p.error('Choose a capability ID or all')
    chosen = rows if a.capability == 'all' else [r for r in rows if r['id'] == a.capability]
    if a.action == 'test':
        if a.output: p.error('--output is only for run')
        for r in chosen:
            args, cwd, env = command(r, 'test')
            subprocess.run(args, cwd=cwd, env=env, check=True, timeout=1200)
    else:
        results = [run_example(r) for r in chosen]
        output = json.dumps(results if a.capability == 'all' else results[0],
                            ensure_ascii=False, indent=2, allow_nan=False) + '\n'
        if a.output:
            a.output.parent.mkdir(parents=True, exist_ok=True)
            a.output.write_text(output)
        else: print(output, end='')


if __name__ == '__main__': main()
