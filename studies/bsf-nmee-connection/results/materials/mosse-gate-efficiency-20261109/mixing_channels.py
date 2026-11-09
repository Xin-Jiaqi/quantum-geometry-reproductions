"""Fixed zero-field layer subspaces for causal hopping interventions.

VV and CC couple occupied or empty states of the decoupled layer blocks;
VC contains both directions between occupied and empty sectors. These are
Hamiltonian controls, not additive pieces of the nonlinear response.
"""
import numpy as np


def adj(a):
    return a.conj().swapaxes(-1, -2)


def split_channels(h, layers, occupied_per_layer, velocity=None):
    """Return decoupled H, VV/CC/VC hoppings and optional exact derivatives.

    Projectors are constructed at B=0, independently in each layer, and must
    stay fixed when the Zeeman probe is applied. Within-sector degeneracies
    cause no ambiguity. A closed layer gap is rejected.
    """
    h = np.asarray(h, dtype=complex)
    n = h.shape[-1]
    assert sorted(np.concatenate(layers).tolist()) == list(range(n))
    p = np.zeros_like(h)
    intra = np.zeros_like(h)
    dp = np.zeros_like(h) if velocity is not None else None
    vintra = np.zeros_like(h) if velocity is not None else None
    gaps = []
    for ids in layers:
        ids = np.asarray(ids)
        block = h[..., ids[:, None], ids]
        e, u = np.linalg.eigh(block)
        no = occupied_per_layer
        gaps.append(float(np.min(e[..., no] - e[..., no - 1])))
        if gaps[-1] < 1e-8:
            raise ValueError('A decoupled layer has no separating valence gap')
        occ = u[..., :no]
        p[..., ids[:, None], ids] = occ @ adj(occ)
        intra[..., ids[:, None], ids] = block
        if velocity is not None:
            vb = velocity[..., ids[:, None], ids]
            vintra[..., ids[:, None], ids] = vb
            f = (np.arange(len(ids)) < no).astype(float)
            df = f[:, None] - f[None, :]
            de = e[..., :, None] - e[..., None, :]
            ratio = np.divide(df, de, out=np.zeros_like(de), where=df != 0)
            dp[..., ids[:, None], ids] = u @ (ratio * (adj(u) @ vb @ u)) @ adj(u)
    q = np.eye(n) - p
    t = h - intra
    vv, cc = p @ t @ p, q @ t @ q
    channels = {'VV': vv, 'CC': cc, 'VC': t - vv - cc}
    diagnostics = {'layer_gaps_eV': gaps,
                   'projector_idempotency': float(np.max(abs(p @ p - p))),
                   'reconstruction_eV': float(np.max(abs(intra + sum(channels.values()) - h)))}
    if velocity is None:
        return intra, channels, p, diagnostics
    dt = velocity - vintra
    dvv = dp @ t @ p + p @ dt @ p + p @ t @ dp
    dcc = -dp @ t @ q + q @ dt @ q - q @ t @ dp
    derivatives = {'VV': dvv, 'CC': dcc, 'VC': dt - dvv - dcc}
    return intra, channels, p, diagnostics, vintra, derivatives
