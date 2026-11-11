#!/usr/bin/env python3
"""read_spn (R1.3): parse a wannier90/VASP .spn file (formatted ASCII or
unformatted binary) and build the Wannier-basis spin operator.

DEFINITIVE FORMAT (sources: wannier90 conv_write_spn; VASP2WAN90_v2_fix patch
mlwf.F for VASP 5.4.4; VASP wiki LWRITE_SPN for VASP > 6.4.2; wannierberri
production parser):
  line 1 / record 1   : header string
  line 2 / record 2   : "nbnd nk" (num_bands, num_kpts)
  per k-point         : 3 * ntri values, ntri = nbnd*(nbnd+1)/2, in the pair
                        order  (m=1..nbnd, n=1..m)  [m outer, n inner],
                        three spin components CONSECUTIVE per pair:
                        sigma_x, sigma_y, sigma_z of the matrix element
                        S[n, m] = <u_nk|S_ipol|u_mk>  (n <= m, i.e. the
                        UPPER-triangle element is stored, row=n, col=m).
  The full matrix is Hermitian-completed: S[m, n] = conj(S[n, m]) for n < m,
  and the diagonal is real (checked).
  BASIS: .spn is in the BLOCH EIGENSTATE basis (num_bands x num_bands),
  NOT the Wannier basis.

Wannier-basis operator (observable): S^W_ipol(k) = U(k)^dagger S(k) U(k),
U(k) = num_bands x num_wann Wannier matrix from the interface.

FORMATTED vs UNFORMATTED:
  - VASP 5.4.4 VaspWannier patch: default UNFORMATTED (binary); set
    INCAR  LSPN_FMTED = .TRUE.  for ASCII. WRITE_SPN is SERIAL-ONLY in that
    patch (WDES%COMM%NCPU==1 required; silently skipped otherwise).
  - VASP > 6.4.2 native LWRITE_SPN: always FORMATTED, auto-adds
    'spn_formatted = .true.' to the win.
  - wannier90 conv_write_spn: UNFORMATTED binary.
parse_spn auto-detects by peeking the file.

Unit tests (run: python3 read_spn.py --selftest) use INDEPENDENTLY constructed
official-format fixtures (hand-written text and packed binary), including a
sigma_y sign test that distinguishes the triangle orientation (a self-loop
writer/parser test alone is NOT accepted).
"""
import json
import struct
import sys
from pathlib import Path

import numpy as np

# known 4-band Pauli matrices (spin blocks) for the independent sign tests
SX4 = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
               dtype=complex)
SY4 = np.array([[0, -1j, 0, 0], [1j, 0, 0, 0],
                [0, 0, 0, -1j], [0, 0, 1j, 0]], dtype=complex)
SZ4 = np.diag([1.0, 1.0, -1.0, -1.0])


def _is_formatted(raw):
    """Formatted .spn always has: line1 = header string, line2 = 'nbnd nk'
    (two integers). A binary (unformatted) file has no newline-separated
    second line that parses as two ints (its header record is raw chars)."""
    head = raw.split(b"\n", 2)
    if len(head) < 2:
        return False
    try:
        tok = head[1].strip().split()
        if len(tok) >= 2:
            int(tok[0]); int(tok[1])
            return True
    except ValueError:
        pass
    return False


def parse_spn(path, nk_expected=None, fmt="auto"):
    """Returns dict: num_bands, n_k, data[ik, m, n, ipol] (complex).
    fmt: 'auto' | 'formatted' | 'unformatted'."""
    raw = Path(path).read_bytes()
    if fmt == "auto":
        fmt = "formatted" if _is_formatted(raw) else "unformatted"
    if fmt == "formatted":
        return _parse_formatted(raw, nk_expected, path)
    return _parse_unformatted(raw, nk_expected, path)


def _parse_formatted(raw, nk_expected, path):
    lines = [ln for ln in raw.decode("utf-8", errors="replace").splitlines()
             if ln.strip()]
    try:
        nbnd, nk = (int(x) for x in lines[1].split()[:2])
    except (IndexError, ValueError):
        raise SystemExit(f"FAIL: cannot parse .spn header of {path} (want 'nbnd nk')")
    return _assemble(lines[2:], nbnd, nk, nk_expected, path)


def _parse_unformatted(raw, nk_expected, path):
    """Fortran sequential unformatted records (ifort/gfortran default):
    each record = [int32 len][data][int32 len] (little-endian; big-endian
    fallback). Record 1: header chars; record 2: (nbnd, nk) int32;
    per k: complex128[3*ntri]."""
    pos, records = 0, []
    n = len(raw)
    while pos < n:
        if pos + 8 > n:
            raise SystemExit(f"FAIL: truncated record marker in {path}")
        (l1,) = struct.unpack_from("<i", raw, pos)
        big = False
        if not (0 <= l1 <= n - pos - 8):
            (l1,) = struct.unpack_from(">i", raw, pos)
            big = True
        if not (0 <= l1 <= n - pos - 8):
            raise SystemExit(f"FAIL: bad record length {l1} in {path}")
        data = raw[pos + 4: pos + 4 + l1]
        (l2,) = struct.unpack_from(("<i" if not big else ">i"), raw, pos + 4 + l1)
        if l2 != l1:
            raise SystemExit(f"FAIL: record marker mismatch in {path}")
        records.append(data)
        pos += 8 + l1
    if len(records) < 2:
        raise SystemExit(f"FAIL: too few records in {path}")
    nbnd, nk = struct.unpack("<ii", records[1])
    ntri = nbnd * (nbnd + 1) // 2
    nvals = 3 * ntri
    lines = []
    for ik in range(nk):
        rec = records[2 + ik]
        if len(rec) != 16 * nvals:
            raise SystemExit(
                f"FAIL: k-block {ik}: {len(rec)} bytes, expected {16 * nvals} "
                f"(nbnd={nbnd})")
        vals = np.frombuffer(rec, dtype="<c16")
        lines.extend(f"{v.real: .14f} {v.imag: .14f}" for v in vals)
    return _assemble(lines, nbnd, nk, nk_expected, path)


def _assemble(lines, nbnd, nk, nk_expected, path):
    """lines: per-k blocks of 're im' lines, component-fastest per pair,
    pair order (m outer, n=1..m inner), value = upper element S[n, m]."""
    ntri = nbnd * (nbnd + 1) // 2
    n_lines = 3 * ntri
    rows, cols = np.tril_indices(nbnd)   # pair j: (m_j, n_j) = (rows, cols)
    data = np.zeros((nk, nbnd, nbnd, 3), dtype=complex)
    idx = 0
    for ik in range(nk):
        block = lines[idx:idx + n_lines]
        if len(block) != n_lines:
            raise SystemExit(
                f"FAIL: k-block {ik} has {len(block)} lines, expected {n_lines} "
                f"(nbnd={nbnd}) in {path}")
        tmp = np.array([ln.split() for ln in block], dtype=float)
        if tmp.shape[1] != 2:
            raise SystemExit(f"FAIL: k-block {ik}: lines must be 're im' pairs")
        vals = tmp[:, 0] + 1j * tmp[:, 1]           # (3*ntri,) component-fastest
        A = np.zeros((3, nbnd, nbnd), dtype=complex)
        # value v_j = S[n_j, m_j] -> upper slot (row=n_j, col=m_j)
        A[:, cols, rows] = vals.reshape(3, ntri, order="F")
        A[:, rows, cols] = A[:, cols, rows].conj()  # Hermitian completion
        diag_im = float(np.einsum("ijj->", np.abs(A.imag)))
        if diag_im > 1e-10:
            raise SystemExit(
                f"FAIL: k-block {ik}: non-real diagonal (sum|Im|={diag_im:.2e})")
        data[ik] = A.transpose(1, 2, 0)
        idx += n_lines
    if idx != len(lines):
        raise SystemExit(f"FAIL: {len(lines) - idx} trailing lines in {path}")
    if nk_expected is not None and nk != nk_expected:
        raise SystemExit(f"FAIL: nk={nk}, expected {nk_expected}")
    return {"num_bands": nbnd, "n_k": nk, "data": data, "fmt": "parsed"}


def wannier_spin(spn, U):
    """S^W(k) = U(k)^dagger S(k) U(k); U: (n_k, num_bands, num_wann)."""
    data = spn["data"]
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


# ----------------------------------------------------------------------
# INDEPENDENT official-format fixtures (no writer involved)
# ----------------------------------------------------------------------
def official_text_spn(matrices):
    """matrices: (nk, nbnd, nbnd, 3) complex. Emit the OFFICIAL formatted
    layout by hand: header, 'nbnd nk', then per k the pairs (m=1..nbnd,
    n=1..m) with 3 components consecutive of S[n, m] = M[n, m] (n <= m)."""
    nbnd = matrices.shape[1]
    lines = ["Written by BD36 independent fixture", f"{nbnd} {matrices.shape[0]}"]
    rows, cols = np.tril_indices(nbnd)
    for ik in range(matrices.shape[0]):
        for j in range(len(rows)):
            n, m = cols[j], rows[j]          # pair (m, n); upper element (n, m)
            for ipol in range(3):
                v = matrices[ik, n, m, ipol]  # S[n, m] with n <= m
                lines.append(f"{v.real: .14f} {v.imag: .14f}")
    return "\n".join(lines) + "\n"


def official_binary_spn(matrices):
    """Packed Fortran unformatted records: header chars, (nbnd, nk) int32,
    per k complex128[3*ntri] (values same official pair order)."""
    nbnd, nk = matrices.shape[1], matrices.shape[0]
    rows, cols = np.tril_indices(nbnd)
    header = b"Written by BD36 independent binary fixture" + b"\x00" * 60
    recs = [header, struct.pack("<ii", nbnd, nk)]
    for ik in range(nk):
        vals = []
        for j in range(len(rows)):
            n, m = cols[j], rows[j]
            for ipol in range(3):
                v = matrices[ik, n, m, ipol]
                vals.append(complex(v.real, v.imag))
        recs.append(b"".join(struct.pack("<dd", v.real, v.imag) for v in vals))
    out = b""
    for r in recs:
        out += struct.pack("<i", len(r)) + r + struct.pack("<i", len(r))
    return out


def selftest():
    # --- independent sigma_y SIGN test (the convention discriminator) ---
    nbnd, nk = 2, 1
    mats2 = np.zeros((nk, nbnd, nbnd, 3), dtype=complex)
    sy2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sx2 = np.array([[0, 1], [1, 0]], dtype=complex)
    sz2 = np.diag([1.0, -1.0]).astype(complex)
    mats2[0, :, :, 0] = sx2
    mats2[0, :, :, 1] = sy2
    mats2[0, :, :, 2] = sz2
    # hand-written official text:
    # pairs (m,n): (1,1): S[0,0]; (2,1): S[0,1]. per pair x,y,z consecutive.
    txt = ("official sigma_y sign fixture\n"
           "2 1\n"
           # pair (1,1): S[0,0] = (sx[0,0]=0, sy[0,0]=0, sz[0,0]=1)
           " 0.0 0.0\n 0.0 0.0\n 1.0 0.0\n"
           # pair (2,1): S[0,1] = (sx[0,1]=1, sy[0,1]=-i, sz[0,1]=0)
           #   sigma_y[0,1] = <e0|sigma_y|e1> = -i  (upper element of Mphys)
           " 1.0 0.0\n 0.0 -1.0\n 0.0 0.0\n"
           # pair (2,2): S[1,1] = (sx[1,1]=0, sy[1,1]=0, sz[1,1]=-1)
           " 0.0 0.0\n 0.0 0.0\n -1.0 0.0\n")
    fa = Path(__file__).parent / ".." / "LOGS" / "r13_spn_sy2_official.spn"
    fa.write_text(txt)
    spn = parse_spn(fa, nk_expected=1)
    err = np.max(np.abs(spn["data"] - mats2))
    assert err < 1e-12, f"official-format sigma_y sign FAILED (err {err})"
    assert abs(spn["data"][0, 0, 1, 1] - (-1j)) < 1e-12, \
        "sigma_y upper element M[0,1] must be -i (sigma_y[0,1]=<e0|sy|e1>=-i); " \
        "got {spn['data'][0,0,1,1]} - triangle orientation / conjugation wrong"
    # --- 4-band sigma_y ⊗ I2 official text (10 pairs) ---
    nbnd4, nk4 = 4, 1
    mats4 = np.zeros((nk4, nbnd4, nbnd4, 3), dtype=complex)
    mats4[0, :, :, 0] = SX4
    mats4[0, :, :, 1] = SY4
    mats4[0, :, :, 2] = SZ4
    fb = Path(__file__).parent / ".." / "LOGS" / "r13_spn_4band_official.spn"
    fb.write_text(official_text_spn(mats4))
    spn4 = parse_spn(fb, nk_expected=1)
    assert np.max(np.abs(spn4["data"] - mats4)) < 1e-12, "4-band official FAIL"
    # --- binary (unformatted) fixture ---
    fbb = Path(__file__).parent / ".." / "LOGS" / "r13_spn_4band_binary.spn"
    fbb.write_bytes(official_binary_spn(mats4))
    spnb = parse_spn(fbb, nk_expected=1)
    assert np.max(np.abs(spnb["data"] - mats4)) < 1e-12, "binary parse FAIL"
    # --- Wannier gauge transform (permutation U) against independent numpy ---
    U = np.array([[0, 0], [0, 1], [1, 0], [0, 0]], dtype=complex)
    SW = wannier_spin(spn4, U)
    SW_ref = np.zeros_like(SW)
    for ipol in range(3):
        SW_ref[0, :, :, ipol] = U.conj().T @ mats4[0, :, :, ipol] @ U
    assert np.max(np.abs(SW - SW_ref)) < 1e-12, "Wannier transform mismatch"
    assert np.max(np.abs(SW - SW.conj().transpose(0, 2, 1, 3))) < 1e-12, \
        "S^W not Hermitian"
    print("read_spn selftest: PASS (official sigma_y sign, 4-band, binary, "
          "Wannier transform, Hermiticity)")
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
        "basis": "Bloch eigenstates <u_mk|S_ipol|u_nk>; upper-triangle S[n,m] "
                 "(n<=m) stored, Hermitian-completed",
        "wannier_transform": "S_WF(k)=U(k)^dagger S(k) U(k)",
    }
    print(json.dumps(rep, indent=2, sort_keys=True))
    # stdout is intentionally strict JSON: downstream Phase-0 uses json.load().
    # Status is carried by exit code; mixing a status line would corrupt JSON.
    sys.exit(0 if rep["hermitian"] else 1)


if __name__ == "__main__":
    main()
