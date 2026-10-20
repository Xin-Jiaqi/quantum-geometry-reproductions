"""Finite-eta length-gauge response in a fixed orbital Bloch basis.

e=hbar=1; k dimensionless and Fourier field prefactors factored out.
The eta regulator is adiabatic, NOT a scattering model.
Computes the full rho2 at finite eta (including Drude poles), NOT
the intrinsic tau^0 coefficient of the Xiang et al. PRL.
"""
from __future__ import annotations
import numpy as np

def hermitian(A, name='matrix'):
    A=np.asarray(A,dtype=complex)
    if A.ndim!=2 or A.shape[0]!=A.shape[1] or not np.isfinite(A).all():
        raise ValueError(name+' must be finite and square')
    if np.max(np.abs(A-A.conj().T))>1.e-9:
        raise ValueError(name+' must be Hermitian')
    return A

def equilibrium(H, mu=0., temperature=.2):
    """Finite-T Fermi density in a fixed orbital basis; supports doublets."""
    H=hermitian(H,'H')
    if temperature<=0:raise ValueError('temperature must be positive')
    en,U=np.linalg.eigh(H)
    f=1/(1+np.exp(np.clip((en-mu)/temperature,-700,700)))
    return (U*f[None,:])@U.conj().T

def resolvent(H,X,eta=.1,order=1,omega=0.,backend='band'):
    """[(p eta-i omega)+i[H,.]]^-1 X; two algebraically different backends."""
    H=hermitian(H,'H');X=np.asarray(X,dtype=complex)
    N=H.shape[0]
    if X.shape!=H.shape or not np.isfinite(X).all():raise ValueError('X shape/finite')
    if eta<=0 or order not in (1,2):raise ValueError('eta>0 and order=1 or 2 required')
    z=order*eta-1j*omega
    if backend=='band':
        en,U=np.linalg.eigh(H)
        B=U.conj().T@X@U
        return U@(B/(z+1j*(en[:,None]-en[None,:])))@U.conj().T
    if backend=='orbital':
        eye=np.eye(N,dtype=complex)
        L=z*np.eye(N*N,dtype=complex)+1j*(np.kron(eye,H)-np.kron(H.T,eye))
        return np.linalg.solve(L,X.reshape(N*N,order='F')).reshape((N,N),order='F')
    raise ValueError('backend must be band or orbital')

def fd4(fn,x,step=3.e-4):
    """Fourth-order central derivative in a fixed basis."""
    if step<=0:raise ValueError('step positive')
    return (fn(x-2*step)-8*fn(x-step)+8*fn(x+step)-fn(x+2*step))/(12*step)

def second_order(Hfun,kx,ky,observable,mu=0.,temperature=.2,eta=.1,
                 omega=0.,axis=0,step=3.e-4,backend='band'):
    """Symmetrized (+omega,-omega) DC rho2 matrix and scalar observable.

    rho1(omega)=-L1(omega)^-1 d_k rho0;
    rho2(0)=-L2(0)^-1 d_k [rho1(+omega)+rho1(-omega)]/2.
    The response remains finite-eta, with a user-specified preparation.
    A real E(t)=E cos(omega t) needs Fourier prefactors outside this function.
    """
    if axis not in (0,1):raise ValueError('axis must be 0 or 1')
    H0=hermitian(Hfun(kx,ky),'H(k)')
    O=hermitian(observable,'observable')
    if O.shape!=H0.shape:raise ValueError('observable has wrong size')
    k=np.array([kx,ky],float)
    def at(x):
        q=k.copy();q[axis]=x
        return hermitian(Hfun(*q),'H(k)')
    def r1(x,w):
        drho=fd4(lambda y:equilibrium(at(y),mu,temperature),x,step)
        return -resolvent(at(x),drho,eta,1,w,backend)
    ds=fd4(lambda x:(r1(x,omega)+r1(x,-omega))/2,k[axis],step)
    rho2=-resolvent(H0,ds,eta,2,0.,backend)
    val=np.trace(O@rho2)
    return {'value_real':float(val.real),'value_imag':float(val.imag),
            'rho2':rho2,'eta':float(eta),'omega':float(omega),
            'protocol':'adiabatic_finite_eta_full_density_matrix_DC',
            'scientific_scope':'NOT intrinsic tau^0; NOT a collision model',
            'backend':backend}
