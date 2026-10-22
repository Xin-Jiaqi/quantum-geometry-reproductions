#!/usr/bin/env python3
"""Minimal synthetic unit test for parse_vasp_b1_cg.py (no VASP / POTCAR needed).

Runs the built-in selftest which fabricates a tiny 12-atom OUTCAR/OSZICAR/
vasprun.xml/stdout/POSCAR in a temp dir and asserts parsing + geometry metrics.
Run:
    python3 test_parse_vasp_b1_cg.py

Real-run usage (read-only, on the private WSL output directory):
    python3 parse_vasp_b1_cg.py OUTCAR OSZICAR vasprun.xml vasp.stdout POSCAR --outdir OUT
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import parse_vasp_b1_cg as p

def main():
    # synthetic fixtures + assertions are inside p.selftest()
    p.selftest()
    print("test_parse_vasp_b1_cg.py: ALL CHECKS PASSED")

if __name__ == "__main__":
    main()
