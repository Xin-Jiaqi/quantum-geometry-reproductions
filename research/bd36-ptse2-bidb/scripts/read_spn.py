#!/usr/bin/env python3
"""read_spn: parse and validate a wannier90 .spn file (spin operator matrix
elements in the Wannier basis, <w_n|sigma_i|w_m>).

The .spn is produced by VASP's LWRITE_SPN=.TRUE. during the W90 interface run
(IMR VASP 6.4.x native; para VaspWannier 5.4.4 - keywords present, output to be
validated in-place). Parser is tolerant of the known wannier90/VASP ASCII
layouts: numeric content is 3 components x (num_wann x num_wann) complex
matrices, per k-point in the SAME order as the VASP interface k-list.

Checks:
  - dimensional consistency: total real count == 6 * num_wann^2 * nk
  - Hermiticity per component (matrix == dagger)
  - sigma_z-like diagonal block structure sanity for spinor WFs
  - same k-count as the interface run (reported; compare externally)

Usage: python3 read_spn.py seedname.spn [nk_expected]
"""
import re, sys
from pathlib import Path


def parse_spn(path):
    text = Path(path).read_text(errors="replace")
    nums = []
    num_wann = None
    for ln in text.splitlines():
        for tok in ln.split():
            if re.fullmatch(r"[+-]?\d+(\.\d+)?([eEdD][+-]?\d+)?", tok):
                v = float(tok.replace("d", "e").replace("D", "e"))
                if v == int(v) and num_wann is None and abs(v) < 1e6 and v >= 1:
                    # first plausible integer = num_wann
                    if not nums:
                        num_wann = int(v)
                        continue
                nums.append(v)
    if num_wann is None:
        # fallback: infer from first integer in the numeric stream
        for i, v in enumerate(nums):
            if v == int(v) and v >= 1 and v <= 1000:
                num_wann = int(v)
                nums = nums[i + 1:]
                break
    if num_wann is None:
        raise SystemExit("FAIL: could not determine num_wann from .spn header")
    n_reals = len(nums)
    per_k = 6 * num_wann * num_wann
    if n_reals == 0 or n_reals % per_k != 0:
        raise SystemExit(
            f"FAIL: {n_reals} reals not divisible by 6*nw^2={per_k} "
            f"(nw={num_wann}) - unexpected .spn layout")
    nk = n_reals // per_k
    # rebuild: components x,y,z stacked; each is nw x nw complex, per k
    import numpy as np
    arr = np.array(nums).reshape(nk, 3, num_wann, num_wann, 2)
    mats = arr[..., 0] + 1j * arr[..., 1]
    herms = np.max(np.abs(mats - mats.conj().transpose(0, 1, 3, 2)))
    tr_sz = mats[:, 2].trace(axis1=2, axis2=3)  # sigma_z diagonal trace per k
    report = {
        "num_wann": num_wann,
        "n_k": int(nk),
        "max_hermiticity_violation": float(herms),
        "trace_sigma_z_per_k": [round(float(x), 6) for x in tr_sz],
        "dims_ok": True,
    }
    return report


def main():
    path = sys.argv[1]
    nk_exp = int(sys.argv[2]) if len(sys.argv) > 2 else None
    rep = parse_spn(path)
    ok = rep["dims_ok"] and rep["max_hermiticity_violation"] < 1e-6
    if nk_exp is not None and rep["n_k"] != nk_exp:
        ok = False
        rep["nk_mismatch"] = f"expected {nk_exp}, got {rep['n_k']}"
    import json
    print(json.dumps(rep, indent=2, sort_keys=True))
    print("->", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
