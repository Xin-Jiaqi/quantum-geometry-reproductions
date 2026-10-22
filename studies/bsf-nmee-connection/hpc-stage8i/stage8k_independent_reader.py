#!/usr/bin/env python3
"""Stage8k INDEPENDENT WAVECAR reader (no pymatgen), VASP 6.4.3 ncl.

Layout (verified against file size 4,899,393,312 B and EIGENVAL):
  REC1  (114528 B): recl, spin, rtag (3 f64) + pad
  REC2  (114528 B): nk, nb, encut, lat(9), efermi + pad
  k-records (97*114528 B each):
     header recl8 (114528 B): nplane(1) + k(3) + enocc[nb][e,0,occ] + pad
     96 band records (114528 B each): nplane complex64 coeffs + (recl8-nplane) f64 pad
  ncl spinor: 14294 complex64 per band = 7147 G x [up,dn] interleaved (4 float32 per G)
Provides: header, kpoint finder, enocc energies vs EIGENVAL, and spinor spin
expectations under BOTH interleaved and grouped channel layouts.
"""
import numpy as np

RECL = 114528.0
RECL8 = 14316
HDR_BYTES = RECL8 * 8
KP_BYTES = 97 * RECL8 * 8


def read_header(path):
    with open(path, "rb") as f:
        head = f.read(24)
        recl, spin, rtag = np.frombuffer(head, dtype=np.float64)
        recl8 = int(recl / 8)
        f.seek(recl8 * 8)
        rec2 = np.frombuffer(f.read(recl8 * 8), dtype=np.float64)
        nk, nb = int(rec2[0]), int(rec2[1])
        encut = float(rec2[2])
        a = rec2[3:12].reshape(3, 3).tolist()
        efermi = float(rec2[12])
        lat = [a[0], a[1], a[2]]
    hdr = {"recl": float(recl), "spin": int(spin), "rtag": int(rtag), "nk": nk,
           "nb": nb, "encut": encut, "lat": lat, "efermi": efermi,
           "recl8": recl8, "hdr_bytes": recl8 * 8,
           "kp_bytes": (1 + nb) * recl8 * 8, "k0_off": 2 * recl8 * 8}
    return hdr


def kpoint_meta(path, header, kidx):
    """nplane, kvec, enocc(nb,3) for kpoint index kidx (seek only, fast)."""
    off = header["k0_off"] + kidx * header["kp_bytes"]
    nb = header["nb"]
    with open(path, "rb") as f:
        f.seek(off)
        nplane = int(np.frombuffer(f.read(8), dtype=np.float64)[0])
        k = np.frombuffer(f.read(24), dtype=np.float64)
        enocc = np.frombuffer(f.read(3 * nb * 8), dtype=np.float64).reshape(nb, 3)
    return nplane, k, enocc


def band_coeffs(path, header, kidx, band):
    """Read nplane complex64 coefficients for one band (seek only)."""
    off = header["k0_off"] + kidx * header["kp_bytes"] + header["hdr_bytes"] + band * header["hdr_bytes"]
    nplane, _, _ = kpoint_meta(path, header, kidx)
    with open(path, "rb") as f:
        f.seek(off)
        coeffs = np.frombuffer(f.read(nplane * 8), dtype=np.complex64)
    return coeffs, nplane


def split_spinor(coeffs, nplane, layout):
    """Split ncl coefficient stream into up/dn per-G arrays (nG each).
    layout 'interleaved': coeff[2i]=up(G_i), coeff[2i+1]=dn(G_i)
    layout 'grouped': coeff[i]=up(G_i), coeff[i+nG]=dn(G_i)
    """
    nG = nplane // 2
    if layout == "interleaved":
        up = coeffs[0::2]
        dn = coeffs[1::2]
    elif layout == "grouped":
        up = coeffs[:nG]
        dn = coeffs[nG:]
    else:
        raise ValueError(layout)
    return up.astype(np.complex128), dn.astype(np.complex128), nG


def spin_expected(up, dn):
    """Pseudo-wavefunction spin (spin-1/2 units, <S>=sigma/2).
    <Sx>=Re(<u|d>), <Sy>=Im(<u|d>)?? sign conventions:
      sigma_x = 2 Re(u*d) ; sigma_y = 2 Im(u*d) ; sigma_z = |u|^2-|d|^2
    so Sx = Re(u*d), Sy = Im(u*d), Sz = (|u|^2-|d|^2)/2.
    """
    nu = np.vdot(up, up).real
    nd = np.vdot(dn, dn).real
    tot = nu + nd
    ud = np.vdot(up, dn)
    if tot < 1e-12:
        return None, None, None, tot
    return (ud.real / tot, ud.imag / tot, (nu - nd) / tot / 2, tot)


def find_kindex(header, target, kpts):
    """Match fractional k (mod 1, nearest) -> index."""
    d = np.max(np.abs(((kpts - np.asarray(target) + 0.5) % 1) - 0.5), axis=1)
    return int(np.argmin(d)), float(d.min())
