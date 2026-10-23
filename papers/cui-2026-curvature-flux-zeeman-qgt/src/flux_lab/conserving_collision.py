"""Round28: gauge-covariant local-density-conserving linearized dephasing.

Finite Nx x Ny SOC tight-binding torus from Round27; one-body f(H) at T>0.
Collision C=-gamma Q(delta_rho - i[G,rho0]), QX=X-sum_i P_i X P_i.
Q conserves every site's particle density, and is covariant under local U(1)
phase U=exp(iL), because all P_i commute with L. The gauge-covariant
reference term is indispensable: naive dephasing of delta_rho violates the
pure-gauge Ward despite having zero local collision charge source.

The affine linearization relative to a prescribed f(H) is NOT an autonomous
CP thermal Lindblad or a microscopic detailed-balance self-energy/vertex.
It has no phonons, impurity ladders, metallic Drude or thermodynamic limit.

Convention e^{-izt}, Im z>0; hbar=1; all response units are synthetic.
i d_t(delta_rho)=[H,delta_rho]+[V,rho0] + i C.
(z-ad_H+i gamma Q) delta_rho = [V,rho0]+i gamma Q i[G,rho0].
"""
from __future__ import annotations
import numpy as np
from .finite_q_gauge_ward import fermi_density,bond_currents,continuity_residual

def _validate(torus,V,G,gamma,z,mu,temperature):
    torus.check()
    h=torus.H;n=len(h)
    V=np.asarray(V,complex)
    if V.shape!=h.shape or not np.all(np.isfinite(V)):
        raise ValueError("matched finite complex drive V required")
    if G is None:G=np.zeros_like(h)
    else:G=np.asarray(G,complex)
    if (G.shape!=h.shape or not np.all(np.isfinite(G)) or
        not np.allclose(G,np.diag(np.diag(G)),atol=1e-12,rtol=0) or
        not np.allclose(G,G.conj().T,atol=1e-12,rtol=0)):
        raise ValueError("real diagonal Peierls generator G required")
    if not np.isfinite(gamma) or gamma<0:
        raise ValueError("gamma must be nonnegative finite")
    z=complex(z)
    if not np.isfinite(z) or z.imag<=0:
        raise ValueError("retarded Im z >0 required")
    if not np.isfinite(mu) or not np.isfinite(temperature) or temperature<=0:
        raise ValueError("finite mu and T>0 required")
    return h,V,G,z

def local_block(X,projectors):
    """Orthogonal Hilbert-Schmidt projector onto local onsite spin blocks."""
    X=np.asarray(X,complex)
    out=np.zeros_like(X)
    for P in projectors:out+=P@X@P
    return out

def collision_projector(X,projectors):
    return np.asarray(X,complex)-local_block(X,projectors)

def pure_gauge_drive(torus,L,z):
    h=torus.H;L=np.asarray(L,complex)
    return 1j*(L@h-h@L),1j*complex(z)*L

def zeeman_vertex(torus,site_values,*,spin_axis='y',kappa=.7):
    """Physical spin Pauli, not orbital pseudospin; arbitrary local B_i."""
    torus.check()
    spin={'x':np.array([[0,1],[1,0]],complex),
          'y':np.array([[0,-1j],[1j,0]],complex),
          'z':np.diag([1.,-1.]).astype(complex)}
    if spin_axis not in spin:raise ValueError("physical spin axis x/y/z required")
    a=np.asarray(site_values,float)
    if a.shape!=(len(torus.projectors),) or not np.all(np.isfinite(a)) or not np.isfinite(kappa):
        raise ValueError("one finite Zeeman field per site and finite kappa required")
    return -kappa*np.kron(np.diag(a),spin[spin_axis])

def _full_q_matrix(projectors):
    n=projectors.shape[-1]
    I=np.eye(n*n,dtype=complex)
    local=np.zeros((n*n,n*n),complex)
    for P in projectors:
        local+=np.kron(P.T,P)
    return I-local

def collision_response(torus,V,*,G=None,gamma=.12,z=.42+.13j,
                       mu=.2,temperature=.3,method='spectral',
                       conserve_local_charge=True,gauge_reference=True):
    """Return density-matrix response, local densities, currents, Ward error.

    spectral: diagonal unperturbed resolvent with convergent iteration of
    site-block feedback. liouville: independent full orbital n^2 inversion.
    conserve_local_charge=False is an unphysical global RTA negative control.
    gauge_reference=False is a bad-gauge negative control.
    """
    h,V,G,z=_validate(torus,V,G,gamma,z,mu,temperature)
    rho=fermi_density(h,mu,temperature)
    dr_ref=1j*(G@rho-rho@G) if gauge_reference else np.zeros_like(h)
    P=torus.projectors
    qop=(lambda X:collision_projector(X,P)) if conserve_local_charge else (lambda X:X)
    source=V@rho-rho@V
    b=source+1j*gamma*qop(dr_ref)
    n=len(h)
    if method=='spectral':
        e,U=np.linalg.eigh(h)
        den=z+1j*gamma-(e[:,None]-e[None,:])
        def R0(X):
            Xb=U.conj().T@X@U
            return U@(Xb/den)@U.conj().T
        y=R0(b)
        if not conserve_local_charge:
            dr=y
            iterations=1
        else:
            dr=y
            for iterations in range(1,501):
                # (z-ad_H+i gamma I)dr=b+i gamma LOCAL(dr).
                nxt=y+1j*gamma*R0(local_block(dr,P))
                if np.max(abs(nxt-dr))<2e-13:
                    dr=nxt
                    break
                dr=nxt
            else:
                raise RuntimeError("conserving collision fixed-point did not converge")
    elif method=='liouville':
        Lmat=np.kron(np.eye(n),h)-np.kron(h.T,np.eye(n))
        Q=_full_q_matrix(P) if conserve_local_charge else np.eye(n*n)
        mat=z*np.eye(n*n)-Lmat+1j*gamma*Q
        dr=np.linalg.solve(mat,b.reshape(-1,order='F')).reshape((n,n),order='F')
        iterations=1
    else:
        raise ValueError("method must be spectral or liouville")
    diss=-gamma*qop(dr-dr_ref)
    J=bond_currents(torus)
    current_contact=1j*(G[None]@J-J@G[None])
    current=np.einsum('bij,ji->b',J,dr,optimize=True)
    current+=np.einsum('bij,ji->b',current_contact,rho,optimize=True)
    dens=np.einsum('aij,ji->a',P,dr,optimize=True)
    collision_site=np.einsum('aij,ji->a',P,diss,optimize=True)
    ward=continuity_residual(torus,dr,V,rho,z,current_contact)
    return dict(dr=dr,rho=rho,density=dens,current=current,
                collision=diss,collision_site=collision_site,ward=ward,
                source=source,gauge_reference=dr_ref,
                current_contact=current_contact,
                iterations=iterations,collision_rate=gamma,z=z,
                scope="finite synthetic linearized U1 density-conserving toy; not CP bath or metallic KME")
