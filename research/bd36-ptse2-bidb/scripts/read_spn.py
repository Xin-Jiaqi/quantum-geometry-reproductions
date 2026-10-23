#!/usr/bin/env python3
"""read_spn (R1.2): parse a formatted wannier90/VASP .spn file and build the
Wannier-basis spin operator via the gauge (Wannier) transformation.

DEFINITIVE FORMAT (per wannierberri's production parser + VASP wiki LWRITE_SPN;
LWRITE_SPN is officially available for VASP > 6.4.2 only):
  line 1                : header string (ignored)
  line 2                : "nbnd nk"   (num_bands, num_kpts - TWO integers)
  then for each of nk k-points (in the SAME order as the interface k-list):
    3 * ntri lines, ntri = nbnd*(nbnd+1)/2, each "re im":
      order: for triangle index j = 0..ntri-1 (tril_indices, column-major):
             3 consecutive lines: sigma_x(re im), sigma_y(re im), sigma_z(re im)
      filling: lower triangle A[ipol, indn, indm] (row >= col) from the flat
               data reshaped (3, ntri, order='F'); upper triangle = conjugate.
      a REAL diagonal is required: sum |Im A[ipol, i, i]| ~ 0.
  BASIS: SPN.data[ik, m, n, ipol] = <u_{m,k}| S_ipol | u_{n,k}> - the spin
  matrix in the BLOCH EIGENSTATE basis (num_bands x num_bands), NOT directly
  the Wannier basis.

Wannier-basis spin operator (the observable of interest):
  S^W_ipol(k) = U(k)^dagger S_bloch_ipol(k) U(k),  U(k) = num_bands x num_wann
  (the Wannier matrix from the interface; num_wann <= num_bands). Layer-resolved
  spin then uses P_layer (Wannier centers) on S^W.

Unit test (run directly): python3 read_spn.py --selftest
"""
import json
import sys
from pathlib import Path

import numpy as np


def parse_spn(path, nk_expected=None):
    """Returns dict: num_bands, n_k, data[ik, m, n, ipol] (complex)."""
    text = Path(path).read_text(errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    # line 2 (index 1): "nbnd nk"
    try:
        nbnd, nk = (int(x) for x in lines[1].split()[:2])
    except (IndexError, ValueError):
        raise SystemExit(f"FAIL: cannot parse .spn header of {path} (want 'nbnd nk')")
    ntri = nbnd * (nbnd + 1) // 2
    n_lines = 3 * ntri
    indn, indm = np.tril_indices(nbnd)
    data = np.zeros((nk, nbnd, nbnd, 3), dtype=complex)
    idx = 2  # first k-block starts at line index 2
    for ik in range(nk):
        block = lines[idx:idx + n_lines]
        if len(block) != n_lines:
            raise SystemExit(
                f"FAIL: k-block {ik} has {len(block)} lines, expected {n_lines} "
                f"(nbnd={nbnd})")
        tmp = np.array([ln.split() for ln in block], dtype=float)
        if tmp.shape[1] != 2:
            raise SystemExit(f"FAIL: k-block {ik}: lines must be 're im' pairs")
        vals = tmp[:, 0] + 1j * tmp[:, 1]
        A = np.zeros((3, nbnd, nbnd), dtype=complex)
        A[:, indn, indm] = vals.reshape(3, ntri, order="F")
        A[:, indm, indn] = A[:, indn, indm].conj()
        diag_im = float(np.einsum("ijj->", np.abs(A.imag)))
        if diag_im > 1e-10:
            raise SystemExit(f"FAIL: k-block {ik}: non-real diagonal (sum|Im|={diag_im:.2e})")
        data[ik] = A.transpose(1, 2, 0)
        idx += n_lines
    if idx != len(lines):
        raise SystemExit(
            f"FAIL: {len(lines) - idx} trailing lines after {nk} k-blocks")
    if nk_expected is not None and nk != nk_expected:
        raise SystemExit(f"FAIL: nk={nk}, expected {nk_expected}")
    return {"num_bands": nbnd, "n_k": nk, "data": data}


def wannier_spin(spn, U):
    """S^W(k) = U(k)^dagger S(k) U(k); U: (n_k, num_bands, num_wann).
    Returns (n_k, num_wann, num_wann, 3)."""
    data = spn["data"]  # (nk, nbnd, nbnd, 3)
    U = np.asarray(U)
    if U.ndim == 2:
        U = np.broadcast_to(U, (data.shape[0],) + U.shape)
    if U.shape[0] != data.shape[0]:
        raise ValueError("U and spn k-count mismatch")
    nw = U.shape[2]
    out = np.zeros((data.shape[0], nw, nw, 3), dtype=complex)
    for ik in range(data.shape[0]):
        for ipol in range(3):
            out[ik, :, :, ipol] = U[ik].conj().T @ data[ik, :, :, ipol] @ U[ik]
    return out


def write_fixture(path, nbnd, nk, matrices):
    """Write a formatted .spn from an explicit (nk, nbnd, nbnd, 3) complex
    Hermitian set, EXACTLY in the documented layout (for tests)."""
    ntri = nbnd * (nbnd + 1) // 2
    indn, indm = np.tril_indices(nbnd)
    with open(path, "w") as f:
        f.write("Written on ... by BD36 test fixture\n")
        f.write(f" {nbnd} {nk}\n")
        for ik in range(nk):
            A = matrices[ik].transpose(2, 0, 1)  # (3, nbnd, nbnd)
            # wannier90 formatted .spn layout (conv_write_spn convention:
            # per triangle element m, the 3 spin components are written
            # CONSECUTIVELY): component-fastest order -> Fortran flatten.
            flat = A[:, indn, indm].flatten(order="F")
            for v in flat:
                f.write(f"{v.real: .14f} {v.imag: .14f}\n")
    return path


def selftest():
    nbnd, nk = 4, 3
    # known Hermitian Bloch spin matrices: S_z = diag(1,1,-1,-1),
    # S_x = sigma_x x I2, S_y = sigma_y x I2 (nbnd=4, two spin blocks)
    sz = np.diag([1.0, 1.0, -1.0, -1.0])
    sx = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
                  dtype=complex)
    sy = np.array([[0, -1j, 0, 0], [1j, 0, 0, 0],
                   [0, 0, 0, -1j], [0, 0, 1j, 0]], dtype=complex)
    mats = np.stack([sx, sy, sz])          # (3, nbnd, nbnd)
    mats = np.broadcast_to(mats, (nk, 3, nbnd, nbnd)).transpose(0, 2, 3, 1).copy()
    tmp = Path(__file__).parent / ".." / "LOGS" / "test_spn_fixture.spn"
    tmp.parent.mkdir(exist_ok=True)
    write_fixture(tmp, nbnd, nk, mats)
    spn = parse_spn(tmp, nk_expected=nk)
    assert spn["num_bands"] == nbnd
    err = np.max(np.abs(spn["data"] - mats))
    assert err < 1e-10, f"parse roundtrip err {err}"
    # Wannier transform with a known unitary U (nbnd x num_wann): permutation
    # U = [[0,0],[0,1],[1,0],[0,0]] (columns e_2, e_1 in 4-band space).
    U = np.array([[0, 0], [0, 1], [1, 0], [0, 0]], dtype=complex)
    SW = wannier_spin(spn, U)
    SW_ref = np.zeros_like(SW)
    for ik in range(nk):
        for ipol in range(3):
            SW_ref[ik, :, :, ipol] = U.conj().T @ mats[ik, :, :, ipol] @ U
    assert np.max(np.abs(SW - SW_ref)) < 1e-12, "Wannier transform mismatch"
    assert np.max(np.abs(SW - SW.conj().transpose(0, 2, 1, 3))) < 1e-12, \
        "S^W not Hermitian"
    print("read_spn selftest: PASS (roundtrip + Wannier transform + Hermiticity)")
    return 0


def main():
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    path = sys.argv[1]
    nk_exp = int(sys.argv[2]) if len(sys.argv) > 2 else None
    spn = parse_spn(path, nk_expected=nk_exp)
    rep = {
        "num_bands": spn["num_bands"],
        "n_k": spn["n_k"],
        "dims": list(spn["data"].shape),
        "max_abs": float(np.max(np.abs(spn["data"]))),
        "hermitian": bool(np.max(np.abs(
            spn["data"] - spn["data"].conj().transpose(0, 2, 1, 3))) < 1e-6),
        "basis": "Bloch eigenstates <u_mk|S_ipol|u_nk> (NOT Wannier basis)",
        "wannier_transform": "S_WF(k)=U(k)^dagger S(k) U(k); supply U via wannier_spin()",
    }
    print(json.dumps(rep, indent=2, sort_keys=True))
    print("-> PASS" if rep["hermitian"] else "-> FAIL")
    sys.exit(0 if rep["hermitian"] else 1)


if __name__ == "__main__":
    main()
