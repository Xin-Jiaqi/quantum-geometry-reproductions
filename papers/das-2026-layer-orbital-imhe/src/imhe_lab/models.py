"""Independent periodic layer/orbital benchmark; NOT the authors' R5G model."""
from __future__ import annotations
import numpy as np
from .core import Fields, orbital_moment


def mesh(n=48):
    if type(n) is not int or n < 8:
        raise ValueError('mesh n must be an integer >=8')
    axis = -np.pi+2*np.pi*(np.arange(n)+.5)/n
    x, y = np.meshgrid(axis, axis, indexing='ij')
    points = np.column_stack((x.ravel(), y.ravel()))
    # Full torus d^2 k/(2pi)^2, with dimensionless reciprocal coordinates.
    return points, np.full(n*n, 1/(n*n))


def layer_hamiltonian(points, gate=.23, magnetic_order=0.):
    """Three spinless orbitals at z=(-1,0,+1), periodic Hermitian hoppings.

    H(k,magnetic_order)^*=H(-k,-magnetic_order). Gate is a BACKGROUND energy.
    The later u probe independently couples P=diag(-1,0,+1).
    """
    points = np.asarray(points, float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError('points must be finite (Nk,2)')
    if not np.isfinite([gate, magnetic_order]).all():
        raise ValueError('finite background parameters required')
    h = np.zeros((len(points), 3, 3), complex)
    dh = np.zeros((2,)+h.shape, complex)
    ddh = np.zeros((2, 2)+h.shape, complex)
    p = np.diag([-1., 0., 1.])
    h[:] = np.diag([-1.1, .15, 1.35])+gate*p
    h[:, 0, 1] = .42+1j*magnetic_order
    h[:, 1, 0] = .42-1j*magnetic_order
    h[:, 0, 2] = h[:, 2, 0] = .19
    h[:, 1, 2] = h[:, 2, 1] = .37
    def add(i, j, amplitude, r, sine=False):
        r = np.asarray(r)
        phase = points@r
        value = np.sin(phase) if sine else np.cos(phase)
        first = np.cos(phase) if sine else -np.sin(phase)
        second = -value
        z = amplitude*value
        h[:, i, j] += z
        for a in range(2):
            dh[a, :, i, j] += amplitude*r[a]*first
            for b in range(2):
                ddh[a, b, :, i, j] += amplitude*r[a]*r[b]*second
        if i != j:
            h[:, j, i] += z.conj()
            for a in range(2):
                dh[a, :, j, i] += (amplitude*r[a]*first).conj()
                for b in range(2):
                    ddh[a, b, :, j, i] += (amplitude*r[a]*r[b]*second).conj()
    add(0, 0, .24, (1, 0)); add(1, 1, .21, (0, 1)); add(2, 2, .19, (1, 1))
    add(0, 1, .17, (1, 0)); add(0, 1, .28j, (0, 1), True)
    add(0, 2, .13, (1, -1)); add(0, 2, .15j, (1, 0), True)
    add(1, 2, .11, (0, 1)); add(1, 2, .20j, (1, 1), True)
    return h, dh, ddh, np.broadcast_to(p, h.shape).copy()


def layer_fields(points, gate=.23, magnetic_order=0.):
    h, dh, ddh, p = layer_hamiltonian(points, gate, magnetic_order)
    m, dm = orbital_moment(h, dh, ddh)
    return Fields(h, dh, p, np.zeros_like(dh), m, dm)
