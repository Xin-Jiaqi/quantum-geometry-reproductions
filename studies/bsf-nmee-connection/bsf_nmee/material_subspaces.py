"""Response-independent selection of valence and empty material subspaces."""
import numpy as np


def adj(a):
    return a.conj().swapaxes(-1, -2)


def active_layer_basis(h, layers, occupied=14, weighting='amplitude'):
    """Two top valence, two low empty, and two coupled high-empty states/layer.

    Select high-empty combinations from the SVD of their coupling to the
    OTHER layer's top valence states. No response enters this selection.
    'mixing' weights the coupling by the high-empty minus mean-valence
    energy denominator. Returns a full-basis isometry and diagnostics.
    """
    if weighting not in ('amplitude', 'mixing'):
        raise ValueError('Unknown basis selection')
    eig = [np.linalg.eigh(h[:, ids[:, None], ids]) for ids in layers]
    basis = np.zeros((*h.shape[:1], h.shape[-1], 12), complex)
    diagnostics = []
    for i, ids in enumerate(layers):
        e, u = eig[i]; other = layers[1-i]
        eo, uo = eig[1-i]
        top = u[:, :, occupied-2:occupied]
        low = u[:, :, occupied:occupied+2]
        high = u[:, :, occupied+2:]
        vother = uo[:, :, occupied-2:occupied]
        inter = h[:, ids[:, None], other]
        coupling = adj(high)@inter@vother
        candidate = coupling.copy()
        if weighting == 'mixing':
            delta = e[:, occupied+2:] - eo[:, occupied-2:occupied].mean(axis=1)[:, None]
            if delta.min() <= 0:
                raise ValueError('High-empty minus valence separation is not positive')
            candidate /= delta[:, :, None]
        left, singular, _ = np.linalg.svd(candidate, full_matrices=False)
        if singular[:, -1].min() < 1e-10:
            raise ValueError('High-empty coupled sector loses rank; no unique smooth rank-two choice')
        remote = high@left
        selected = np.concatenate([top, low, remote], axis=2)
        basis[:, ids[:, None], np.arange(6*i, 6*i+6)] = selected
        residue = coupling-left@(adj(left)@coupling)
        diagnostics.append(dict(minimum_singular=float(singular[:, -1].min()),
                                discarded_coupling_fraction=float(np.linalg.norm(residue)/np.linalg.norm(coupling)),
                                valence_boundary_gap=float(np.min(e[:, occupied-2]-e[:, occupied-3])),
                                empty_boundary_gap=float(np.min(e[:, occupied+2]-e[:, occupied+1]))))
    return basis, diagnostics


def orthonormalize(q):
    metric = adj(q)@q
    e, u = np.linalg.eigh(metric)
    if e.min() <= 1e-8:
        raise ValueError('Projected active basis loses rank')
    return q@((u/np.sqrt(e)[:, None, :])@adj(u)), float(np.min(e))
