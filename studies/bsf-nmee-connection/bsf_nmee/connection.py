"""Use group-theory-operations-toolkit as the ONLY interface-orbit backend.

The toolkit's fractional rotation acts on stacking shifts; its Cartesian
rotation acts on physical polar/axial response tensors. Never interchange them.
Only unitary, z-normal-preserving point operations and two identically labeled
layers are supported. This is a symmetry *relationship*, not a proof of ferroic
switchability, interface energy degeneracy or a microscopic response value.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping
import numpy as np


@dataclass(frozen=True)
class Connection:
    source_shift: tuple[float, float]
    target_shift: tuple[float, float]
    operation_name: str
    layer_exchanging: bool
    rotation_cartesian: np.ndarray
    rotation_fractional: np.ndarray
    layer_permutation: dict[str, str]
    centering: str

    def as_dict(self):
        return {
            'source_shift': list(self.source_shift),
            'target_shift': list(self.target_shift),
            'operation_name': self.operation_name,
            'layer_exchanging': self.layer_exchanging,
            'rotation_cartesian': self.rotation_cartesian.tolist(),
            'rotation_fractional': self.rotation_fractional.tolist(),
            'layer_permutation': self.layer_permutation,
            'centering': self.centering,
            'scope': 'unitary spatial operation; symmetry-related states only',
        }


def _validate_operation(op):
    """Require toolkit OperationRecord-like objects, with both bases explicit."""
    if not all(hasattr(op, attr) for attr in ('matrix_fractional','matrix_cartesian','name')):
        raise TypeError('Use a group-theory-operations-toolkit OperationRecord with both coordinate matrices')
    rf = np.asarray(op.matrix_fractional, float)
    rc = np.asarray(op.matrix_cartesian, float)
    for label, mat in (('fractional', rf), ('cartesian', rc)):
        if mat.shape != (3, 3) or not np.isfinite(mat).all():
            raise ValueError(f'invalid {label} operation')
        if np.max(np.abs(mat[:2, 2])) > 1e-9 or np.max(np.abs(mat[2, :2])) > 1e-9:
            raise ValueError(f'{label} operation mixes the plane and layer normal')
        if abs(abs(mat[2,2]) - 1) > 1e-9:
            raise ValueError('operation must preserve/reverse layer normal')
    if not np.allclose(rc.T @ rc, np.eye(3), atol=1e-8):
        raise ValueError('physical tensor rotation must be Cartesian and orthogonal')
    if abs(rf[2,2]-rc[2,2]) > 1e-9:
        raise ValueError('fractional and Cartesian layer exchange disagree')
    return rf, rc


def connect_interface(shift, operation, *, centering='P', expected_target=None):
    """Map an input interface shift via toolkit's published BSF interface API.

    `expected_target` checks, but never forces, the computed target. If a
    coordinate setting differs, do NOT relabel by hand: provide the correct
    OperationRecord and Bravais centering. `operation` must be applicable to
    the *monolayer atomic motif* for the connection to be a physical symmetry.
    """
    from group_theory_operations.stacking import equivalent_interface_translation
    rf, rc = _validate_operation(operation)
    tau = np.asarray(shift, float)
    if tau.shape != (2,) or not np.isfinite(tau).all():
        raise ValueError('shift must be a finite 2-component fractional translation')
    mapped = equivalent_interface_translation(tau, operation, centering=centering)
    if expected_target is not None:
        # Canonicalize using the toolkit itself, avoiding duplicate centering logic.
        check = equivalent_interface_translation(expected_target, np.eye(3), centering=centering)
        if not np.allclose(mapped, check, atol=1e-9, rtol=0):
            raise ValueError(f'connection maps {tuple(tau)} to {mapped}, not to {check}')
    exchange = bool(rc[2,2] < 0)
    permutation = ({'top':'bottom', 'bottom':'top'} if exchange else
                   {'top':'top', 'bottom':'bottom'})
    return Connection(tuple(float(x) for x in tau), tuple(map(float,mapped)),
                      str(operation.name), exchange, rc, rf, permutation, centering)


def _transform_one(tensor, r):
    g = np.asarray(tensor)
    if g.shape != (3, 3, 3):
        raise ValueError('NMEE tensor order must be [spin_alpha, B_gamma, E_beta] (3,3,3)')
    if not np.issubdtype(g.dtype, np.number) or not np.isfinite(g).all():
        raise ValueError('tensor must be finite numerical values')
    # Both spin and B are axial: det(R)^2 = 1; E is polar.
    return np.einsum('ia,jb,kc,abc->ijk',r,r,r,g,optimize=True)


def transform_layer_responses(layer_tensors: Mapping[str, np.ndarray], connection: Connection):
    """Predict response at the partner stacking for a unitary connection.

    This predicts the other state; it must NOT be treated as the symmetry
    residual of one structure unless mapped shift equals the original shift.
    """
    if set(layer_tensors) != {'top','bottom'}:
        raise ValueError('expected top and bottom layer tensors')
    return {connection.layer_permutation[source]:
            _transform_one(value, connection.rotation_cartesian)
            for source, value in layer_tensors.items()}


def response_summary(layer_tensors: Mapping[str, np.ndarray]):
    top = np.asarray(layer_tensors['top'])
    bottom = np.asarray(layer_tensors['bottom'])
    return {'total':top+bottom,'layer_difference':top-bottom}
