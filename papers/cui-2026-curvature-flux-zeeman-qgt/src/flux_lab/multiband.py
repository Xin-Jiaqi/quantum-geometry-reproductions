"""Audited NONDEGENERATE Zeeman-QGT primitives.

Bloch Hamiltonian, its k derivatives and the three dimensionless Pauli spin
operators must be in the SAME fixed orthonormal orbital basis. No .spn or
wannier90_hr.dat conversion is inferred; degenerate inputs fail closed.

T[n,m,a,b] = i <n|D_a H|m> <m|sigma_b|n> / (E_m-E_n).
The published four sectors are determined by a/b index exchange, not merely
by taking real/imaginary parts. All returned numbers are natural units.
"""
from __future__ import annotations
import numpy as np


def covariant_h_derivative(h, dh, connection):
    """D_a H = partial_a H - i[A_a,H] (A_a Hermitian position connection)."""
    h = np.asarray(h, dtype=complex)
    dh = np.asarray(dh, dtype=complex)
    connection = np.asarray(connection, dtype=complex)
    if h.ndim != 2 or h.shape[0] != h.shape[1] or dh.shape != (3,*h.shape) or connection.shape != dh.shape:
        raise ValueError('H: (N,N), dH and A: (3,N,N) required')
    return dh - 1j*(connection@h - h@connection)


def _validate(h, velocity, spin, tol):
    h = np.asarray(h, dtype=complex)
    v = np.asarray(velocity, dtype=complex)
    s = np.asarray(spin, dtype=complex)
    if h.ndim != 2 or h.shape[0] != h.shape[1] or h.shape[0] < 2:
        raise ValueError('H must be square, N>=2')
    if v.shape != (3,*h.shape) or s.shape != v.shape:
        raise ValueError('all three velocity and physical spin operators required')
    if not np.all(np.isfinite(h)) or not np.all(np.isfinite(v)) or not np.all(np.isfinite(s)):
        raise ValueError('nonfinite operator')
    if not np.allclose(h,h.conj().T,atol=tol,rtol=0):
        raise ValueError('non-Hermitian H')
    for op in (*v,*s):
        if not np.allclose(op,op.conj().T,atol=tol,rtol=0):
            raise ValueError('non-Hermitian input operator')
    e, u = np.linalg.eigh(h)
    if np.any(np.diff(e) <= tol):
        raise ValueError('degenerate or unresolved bands: use non-Abelian response')
    return e,u,v,s


def zeeman_qgt(h, velocity, spin, *, tol=1e-9, method='spectral'):
    """Return per-(n,m) and per-band fourfold decomposition; NO transport claim.

    'velocity' means D_k H in energy*length (NOT group velocity / hbar).
    'spin' is dimensionless Pauli operator; do not supply hbar*sigma/2.
    Supports 'spectral' or independent 'projector' oracle.
    """
    es, u, v, s = _validate(h,velocity,spin,tol)
    nband = es.size
    T = np.zeros((nband,nband,3,3),complex)
    for n in range(nband):
        pn=np.outer(u[:,n],u[:,n].conj())
        for m in range(nband):
            if m==n: continue
            den=es[m]-es[n]
            if method=='spectral':
                r=np.array([1j*np.vdot(u[:,n],op@u[:,m])/den for op in v])
                si=np.array([np.vdot(u[:,m],op@u[:,n]) for op in s])
                T[n,m]=np.outer(r,si)
            elif method=='projector':
                pm=np.outer(u[:,m],u[:,m].conj())
                T[n,m]=np.array([[1j*np.trace(pn@v[a]@pm@s[b])/den for b in range(3)] for a in range(3)])
            else:
                raise ValueError('method must be spectral or projector')
    total=np.sum(T,axis=1)
    sym=(total+np.swapaxes(total,1,2))/2
    anti=total-np.swapaxes(total,1,2)
    normal_curvature=-anti.imag
    anomalous_curvature=-anti.real
    vector=lambda x: np.stack((x[:,1,2],x[:,2,0],x[:,0,1]),axis=-1)
    return dict(energies=es, pairs=T, tensor=total, gN=sym.real,
                gA=sym.imag, omegaN=vector(normal_curvature),
                omegaA=vector(anomalous_curvature), units='natural',
                scope='nondegenerate Zeeman geometry only')
