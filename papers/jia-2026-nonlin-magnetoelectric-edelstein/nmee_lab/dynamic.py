"""Full mixed E/B response from SM (2), (6), (7), before the DC split.

Complex Fourier amplitudes: s_alpha/mu_B = Gamma * E_beta * B_gamma.
H_B=-g*mu_B*B*S, H_E=E*(i partial_k + A); H in eV, k in nm^-1.
The output frequency is z_E+z_B. Each cosine input has half-amplitude
Fourier components: multiply by 1/4 for each chosen pair of signed drives.
Finite positive Im(z) is adiabatic switching, not a collision lifetime.
"""
from dataclasses import asdict, dataclass
import numpy as np
from scipy.special import expit

from .clean_core.core import KB, MU_B
from .covariant import CovariantBloch2D, _herm, covariant_derivative
from .ptbackend.liouville_general import R


@dataclass(frozen=True)
class DynamicSpec:
    electric_energy_eV: complex
    magnetic_energy_eV: complex
    mu_eV: float = 0.
    temperature_K: float = 100.
    observable: str = 'x'
    magnetic: str = 'z'
    beta: int = 0

    def validate(self):
        for z in (self.electric_energy_eV, self.magnetic_energy_eV):
            if not np.isfinite(z) or complex(z).imag <= 0:
                raise ValueError('Both drive energies require finite, positive imaginary parts')
        if not np.isfinite(self.temperature_K) or self.temperature_K <= 0:
            raise ValueError('Use T>0 for the Fermi derivative; T=0 metallic delta functions are not sampled')
        if not np.isfinite(self.mu_eV) or self.beta not in (0, 1):
            raise ValueError('Finite chemical potential and beta=0/1 required')


def _equilibrium(H, velocity, mu, temperature):
    """f(H) and its Frechet derivative, including exactly degenerate blocks.

    Use the analytic derivative inside a thermal-scale window, rather than
    an absolute energy tolerance which can erase a low-temperature surface.
    """
    e, u = np.linalg.eigh(H)
    ud = u.conj().swapaxes(-1, -2)
    thermal = KB * temperature
    f = expit((mu - e) / thermal)
    delta = e[:, :, None] - e[:, None, :]
    close = abs(delta) < 1e-5 * thermal
    divided = np.zeros_like(delta)
    np.divide(f[:, :, None] - f[:, None, :], delta, out=divided, where=~close)
    x = (mu - (e[:, :, None] + e[:, None, :]) / 2) / thermal
    np.copyto(divided, -expit(x) * expit(-x) / thermal, where=close)
    rho = (u * f[:, None, :]) @ ud
    derivative = u @ (divided * (ud @ velocity @ u)) @ ud
    return e, u, ud, rho, derivative


def evaluate_dynamic(model: CovariantBloch2D, points, spec: DynamicSpec, *, chunk_size=2048):
    """Equal-weight 2D BZ quadrature; exact degeneracies are retained.

    Accepts the same physical H/S/position-connection callbacks as the static
    interface. Returns both time orderings and their sum, in (nm V T)^-1
    for s/mu_B. Mesh and frequency/broadening convergence remain explicit
    caller choices. No intrinsic/extrinsic DC extraction is performed.
    """
    spec.validate()
    pts = np.asarray(points, float)
    if pts.ndim != 2 or pts.shape[1] != 2 or not len(pts) or not np.isfinite(pts).all():
        raise ValueError('Finite points with shape (Nk,2) required')
    if not isinstance(chunk_size, int) or isinstance(chunk_size, bool) or chunk_size < 1:
        raise ValueError('chunk_size must be a positive integer')
    parts = {'B_then_E': 0j, 'E_then_B': 0j}
    for start in range(0, len(pts), chunk_size):
        block = pts[start:start+chunk_size]
        result = _evaluate_block(model, block, spec)
        for name, value in result['parts'].items():
            parts[name] += value * len(block) / len(pts)
    return dict(result, response=sum(parts.values()), parts=parts, mesh_n=len(pts))


def _evaluate_block(model, pts, spec):
    rec = np.asarray(model.reciprocal_nm, float)
    if rec.shape != (2, 2) or not np.isfinite(rec).all() or abs(np.linalg.det(rec)) < 1e-12:
        raise ValueError('Finite nonsingular 2D reciprocal lattice required')
    if not np.isfinite(model.g_factor):
        raise ValueError('Finite g factor required')
    if any(s not in model.operators for s in (spec.observable, spec.magnetic)):
        raise ValueError('Observable and magnetic operators required in the same gauge')
    if spec.magnetic not in model.derivative_operators:
        raise ValueError('Explicit magnetic-operator derivative required (including physical zero)')
    H = _herm(model.hamiltonian(pts), 'H')
    dH = _herm(model.derivative_h(pts, spec.beta), 'dH')
    A = _herm(model.connection(pts, spec.beta), 'position connection')
    S = _herm(model.operators[spec.magnetic](pts), 'magnetic spin')
    O = _herm(model.operators[spec.observable](pts), 'observable')
    dS = _herm(model.derivative_operators[spec.magnetic](pts, spec.beta), 'dS')
    if H.shape[0] != len(pts) or any(x.shape != H.shape for x in (dH, A, S, O, dS)):
        raise ValueError('All matrices must have matching (Nk,N,N) shapes')
    V = covariant_derivative(H, dH, A)
    DS = covariant_derivative(S, dS, A)
    e, u, ud, rho, drho = _equilibrium(H, V, spec.mu_eV, spec.temperature_K)
    def comm(a, b):
        return a @ b - b @ a
    def resolve(x, z):
        return R(e, u, ud, x, z)
    ze, zb = spec.electric_energy_eV, spec.magnetic_energy_eV
    rhoE = 1j * resolve(drho, ze)
    rhoB = resolve(comm(S, rho), zb)
    # D R_B X = R_B ([D H, R_B X] + D X).
    derivativeB = resolve(comm(V, rhoB) + comm(DS, rho) + comm(S, drho), zb)
    densities = {'B_then_E': resolve(1j * derivativeB, ze + zb),
                 'E_then_B': resolve(comm(S, rhoE), ze + zb)}
    measure = abs(float(np.linalg.det(rec))) / (2 * np.pi)**2
    prefactor = -model.g_factor * MU_B * measure
    parts = {name: complex(prefactor * np.einsum('nij,nji->n', O, value).mean())
             for name, value in densities.items()}
    return dict(protocol='finite_frequency_full_EB_adiabatic_covariant',
                response=sum(parts.values()), parts=parts,
                settings=asdict(spec), model=model.label,
                output_energy_eV=complex(ze + zb), mesh_n=len(pts),
                normalization='s/mu_B per E[V/nm] B[T], 2D (nm V T)^-1',
                warnings=['Full mixed response, not isolated intrinsic DC',
                          'Positive imaginary energies denote adiabatic switching, not fitted lifetimes',
                          'Equal-weight BZ sampling; converge the mesh at the stated T and drive energies'])


def cumnas_dynamic_model(theta_rad=0., params=None):
    """Paper Eq(5), physical k in nm^-1; fixed-orbital, A=0 convention.

    Operators x/y/z are total Pauli spin; A_x/... and B_x/... are sublattice
    spin. The paper does not fix a physical lattice length: a_nm is an input.
    """
    from .ptbackend import cumnas as cm
    p = cm.Model() if params is None else params
    if not isinstance(p, cm.Model) or not np.isfinite(p.a_nm) or p.a_nm <= 0 or not np.isfinite(theta_rad):
        raise ValueError('CuMnAs Model, positive finite a_nm, and finite theta required')
    ops = dict(cm.S)
    for label, projector in (('A', cm.PA), ('B', cm.PB)):
        for axis, spin in cm.S.items():
            ops[f'{label}_{axis}'] = (projector @ spin + spin @ projector) / 2
    def zero(k, beta):
        return np.zeros((len(k), 4, 4), complex)
    return CovariantBloch2D(
        hamiltonian=lambda k: cm.hamiltonian(k * p.a_nm, theta_rad, p),
        derivative_h=lambda k, b: cm.velocity(k * p.a_nm, theta_rad, b, p),
        connection=zero,
        operators={s: (lambda k, op=op: np.broadcast_to(op, (len(k), 4, 4))) for s, op in ops.items()},
        derivative_operators={s: zero for s in ops},
        reciprocal_nm=2 * np.pi * np.eye(2) / p.a_nm,
        g_factor=p.g, label='CuMnAs_Eq5_fixed_orbital')
