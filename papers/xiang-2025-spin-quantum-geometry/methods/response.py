"""Full linear/quadratic electric response, before a tau-power decomposition.

SM(1),(2): rho1=L_z^-1 D f(H), rho2=L_(z1+z2)^-1 D rho1,
L_z=-i z+i[H,.]. Natural e=hbar=1; thermal_energy is k_B T.
All matrices use one orbital gauge; D=partial-i[A,.].
"""
from dataclasses import dataclass
import numpy as np
from .finite_eta import equilibrium, resolvent, fd4, hermitian


@dataclass(frozen=True)
class SpinModel:
    hamiltonian: object  # H(kx,ky)
    connection: object  # A(kx,ky,axis), explicit zero allowed for a fixed basis
    spin: object  # physical Pauli spin(kx,ky,component)
    reciprocal: object  # 2x2 rows, reciprocal units matching k


def _derivative(model, fn, k, axis, step):
    def at(x):
        q=np.array(k,float);q[axis]=x
        return fn(q)
    A=hermitian(model.connection(*k,axis),'position connection')
    value=fn(k)
    return fd4(at,k[axis],step)-1j*(A@value-value@A)


def spin_observable(model,k,*,component='z',current_axis=None,step=3e-4):
    """s^c=sigma^c/2 or conventional J^{a;c}={D_a H,sigma^c}/4.

    This evaluates the full velocity operator (including scalar curvature).
    It does not replace a nonlinear H by its linear-Dirac vertex.
    """
    S=hermitian(model.spin(*k,component),'physical Pauli spin')
    if current_axis is None:return S/2
    if current_axis not in (0,1):raise ValueError('current_axis=0/1 or None')
    V=_derivative(model,lambda q:hermitian(model.hamiltonian(*q)),k,current_axis,step)
    return (V@S+S@V)/4


def response_at_k(model,k,*,energy1,axis1=0,energy2=None,axis2=None,
                  mu=0.,thermal_energy=.1,component='z',current_axis=None,
                  step=3e-4,backend='band'):
    """Linear or symmetric second-order coefficient per Fourier amplitude.

    For two labeled drives the quadratic cross term is 2*K*E1*E2;
    for identical drives the second harmonic is K*E1**2. Cosine Fourier
    amplitudes are E/2. No scattering or intrinsic-DC coefficient is inferred.
    Exact degeneracies are permitted by the full density-matrix resolvent.
    """
    k=np.asarray(k,float)
    if k.shape!=(2,) or not np.isfinite(k).all():raise ValueError('finite 2D k required')
    if axis1 not in (0,1) or (energy2 is not None and axis2 not in (0,1)):
        raise ValueError('supply an axis=0/1 for each drive')
    if not np.isfinite([mu,thermal_energy,step]).all() or min(thermal_energy,step)<=0:
        raise ValueError('finite mu, positive thermal_energy and derivative step required')
    for z in (energy1,) if energy2 is None else (energy1,energy2):
        if not np.isfinite(z) or complex(z).imag<=0:raise ValueError('retarded complex energies require Im(z)>0')
    def H(q):return hermitian(model.hamiltonian(*q),'H')
    def rho0(q):return equilibrium(H(q),mu,thermal_energy)
    def solve(q,X,z):
        return resolvent(H(q),X,eta=complex(z).imag,omega=complex(z).real,backend=backend)
    def first(q,z,axis):
        return solve(q,_derivative(model,rho0,q,axis,step),z)
    if energy2 is None:
        rho=first(k,energy1,axis1);out_energy=energy1
    else:
        a=_derivative(model,lambda q:first(q,energy1,axis1),k,axis2,step)
        b=_derivative(model,lambda q:first(q,energy2,axis2),k,axis1,step)
        out_energy=energy1+energy2
        rho=solve(k,(a+b)/2,out_energy)
    O=spin_observable(model,k,component=component,current_axis=current_axis,step=step)
    if O.shape!=rho.shape:raise ValueError('observable and H have different sizes')
    return dict(response=complex(np.trace(O@rho)),density=rho,
                output_energy=complex(out_energy),
                protocol='full_adiabatic_electric_linear' if energy2 is None else 'full_adiabatic_electric_quadratic_symmetric',
                units='natural e=hbar=1; spin=sigma/2, current={D H,sigma}/4',
                scope='Finite-frequency full response, NOT intrinsic tau^0')


def integrate_response(model,points,**settings):
    """Equal-weight 2D BZ integral. Caller converges mesh and derivative step."""
    pts=np.asarray(points,float);rec=np.asarray(model.reciprocal,float)
    if pts.ndim!=2 or pts.shape[1]!=2 or not len(pts):raise ValueError('nonempty 2D grid required')
    if rec.shape!=(2,2) or not np.isfinite(rec).all() or abs(np.linalg.det(rec))<1e-12:
        raise ValueError('nonsingular 2D reciprocal lattice required')
    values=[response_at_k(model,k,**settings)['response'] for k in pts]
    return dict(response=complex(np.mean(values)*abs(np.linalg.det(rec))/(2*np.pi)**2),
                settings=settings,mesh_n=len(pts),normalization='integral d2k/(2pi)^2',
                scope='Full response at specified finite energies; not isolated intrinsic DC')
