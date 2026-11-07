"""Stage 8: explicitly varied bilayers; frozen stage 2--7 sources are untouched.

The default reproduces SlidingBilayer. `direct` means identical complete
monolayer Hamiltonians, not merely AA orbital coordinates. Anisotropy changes
one intralayer nearest-neighbour hopping, preserving the layer mirror.
"""
from dataclasses import dataclass
import numpy as np
from nmee_lab.clean_core.models import kane_mele, I2, SX, SY, SZ
from .sliding_tight_binding import SlidingBilayer, SlidingSpec, M_SPIN, M_GLOBAL
from .valley_quadrature import valley_quadrature


@dataclass(frozen=True)
class AuditSpec(SlidingSpec):
    staggered_eV: float = .010
    stacking: str = 'mirror'
    bond_anisotropy_eV: float = 0.
    same_scale: float = 1.
    cross_scale: float = 1.


class AuditBilayer(SlidingBilayer):
    def __init__(self, spec=AuditSpec()):
        if spec.stacking not in ('mirror', 'direct'):
            raise ValueError('stacking must be mirror or direct')
        if not np.isfinite([spec.staggered_eV, spec.bond_anisotropy_eV,
                            spec.same_scale, spec.cross_scale]).all():
            raise ValueError('nonfinite model parameter')
        super().__init__(spec)
        self.mono = kane_mele(bond_nm=spec.bond_nm,
                             rashba_soc=spec.rashba_soc_eV,
                             intrinsic_soc=spec.intrinsic_soc_eV,
                             staggered_potential=spec.staggered_eV)
        self.hops = [(a,b,d,t*(spec.same_scale if a == b else spec.cross_scale))
                     for a,b,d,t in self.hops]

    def matrices(self, k, axis=None):
        k = np.atleast_2d(np.asarray(k, float))
        h = self.mono.hamiltonian(k) if axis is None else self.mono.velocity(k, axis)
        # First A--B bond is (0,bond_nm) in the full orbital embedding.
        d = np.array([0., self.spec.bond_nm])
        correction = -self.spec.bond_anisotropy_eV*np.exp(1j*k@d)
        if axis is not None:
            correction *= 1j*d[axis]
        h[:, :2, 2:] += correction[:, None, None]*I2
        h[:, 2:, :2] += correction.conj()[:, None, None]*I2
        lower = M_SPIN@h@M_SPIN if self.spec.stacking == 'mirror' else h
        t = self.tunneling(k, axis)
        result = np.zeros((len(k),8,8), complex)
        result[:, :4, :4], result[:, 4:, 4:] = h, lower
        result[:, :4, 4:], result[:, 4:, :4] = t, t.conj().swapaxes(1,2)
        return result


def integration_grid(model, exterior=61, radial=60, angular=96):
    """C3-closed quadrature of the BZ, NOT a symmetry projection of the result.

Every rotated k is evaluated independently, including for anisotropic models.
This prevents the exterior grid from leaking linear errors into cubic scalars.
It remains necessary to refine the radial/exterior quadrature.
"""
    k, w, meta = valley_quadrature(model.model().reciprocal,
                                  exterior, 1.5, radial, angular)
    theta = 2*np.pi/3
    r = np.array([[np.cos(theta),-np.sin(theta)],
                  [np.sin(theta),np.cos(theta)]])
    return np.concatenate([k,k@r.T,k@(r@r).T]), np.tile(w/3,3), meta


def symmetry_residuals(spec):
    """Off-grid Hamiltonian checks, independent of integrated response fitting."""
    from dataclasses import replace
    k = np.array([[.13,.29],[8.1,-3.2],[-5.8,12.2]])
    obj = AuditBilayer(spec)
    h = obj.matrices(k)
    theta = 2*np.pi/3
    r = np.array([[np.cos(theta),-np.sin(theta)],
                  [np.sin(theta),np.cos(theta)]])
    spin_rotation = np.kron(np.eye(4), np.diag(np.exp(-.5j*theta*np.array([1,-1]))))
    d_rot = np.array(spec.shift)@obj.a@r.T@np.linalg.inv(obj.a)
    rotated = AuditBilayer(replace(spec, shift=tuple(d_rot))).matrices(k@r.T)
    mirrored = AuditBilayer(replace(spec, shift=tuple(-np.array(spec.shift)))).matrices(k)
    time_reversal = np.kron(np.eye(4),1j*SY)
    inversion = np.kron(SX,np.kron(SX,I2))
    # Single-layer generators: C3 and vertical m_x survive, horizontal m_z does not.
    mono = obj.mono.hamiltonian(k)
    mx = np.kron(I2,SX)
    mono_mx = obj.mono.hamiltonian(k*np.array([-1,1]))
    return {
        'C3_eV': float(np.max(abs(rotated-spin_rotation@h@spin_rotation.conj().T))),
        'mirror_pair_eV': float(np.max(abs(mirrored-M_GLOBAL@h@M_GLOBAL))),
        'time_reversal_eV': float(np.max(abs(obj.matrices(-k)-time_reversal@h.conj()@time_reversal.conj().T))),
        'graphene_inversion_eV': float(np.max(abs(obj.matrices(-k)-inversion@h@inversion))),
        'single_layer_vertical_mirror_eV': float(np.max(abs(mono_mx-mx@mono@mx))),
        'single_layer_horizontal_mirror_eV': float(np.max(abs(mono-M_SPIN@mono@M_SPIN))),
    }
