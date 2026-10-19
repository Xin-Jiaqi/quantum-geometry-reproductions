"""Unitary real-space layer-permutation bookkeeping for NMEE E*B tensors.

Only crystallographic unitary spatial transformations. Magnetic/antiunitary
operations require extra spin/time-reversal transformation rules.
"""
import numpy as np
from .clean_core.symmetry import transform_nmee_tensor


def transform_layer_tensors(tensors, rotation, layer_permutation):
    """Connection operator between stacking states.

    tensors[source_layer] = Gamma[spin_alpha, B_gamma, E_beta] (3,3,3).
    layer_permutation[source_layer] = destination_layer.
    For spin and B axial, E polar: two det(R) factors cancel.
    Return tensor dict in the *destination* stacking/layer labeling.
    """
    names=set(tensors)
    if set(layer_permutation)!=names or set(layer_permutation.values())!=names:
        raise ValueError('layer permutation must be bijection of exactly all input layers')
    return {layer_permutation[src]:transform_nmee_tensor(tensor,rotation)
            for src,tensor in tensors.items()}


def connection_residual(tensors,rotation,layer_permutation):
    """Check tensor covariance of a *single* symmetry-invariant stacking state.
    Do NOT use as symmetry test between two distinct ferroic switching partners;
    for those states compare mapped tensors to independently computed tensors.
    """
    images=transform_layer_tensors(tensors,rotation,layer_permutation)
    return max(float(np.max(abs(np.asarray(tensors[k])-images[k]))) for k in tensors)
