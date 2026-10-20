"""Independent metal Fermi-domain quadrature for Xiang et al. PRL Eq.(17).

At T=0 the occupation factor f_- - f_+ = Theta(E(k)-|mu|).
For fixed ky the Fermi boundary in u=cos(a*kx) solves a *quadratic*,
which gives *all* crossings even when E(kx,ky) is nonmonotonic.
The outer ky integral is split at every Fermi topology change where
E(kx=0,ky), E(kx=pi/a,ky), or the x-direction extremum equals |mu|.

The second calculation implements the paper's Fermi-surface identity
C_xy = -sum_n integral (df/d epsilon_n) v^x_n s^y_n
     = integral_delta(E-|mu|) (partial_x E) n_y
using one-dimensional contour roots, without the occupied-sea C kernel.

All wavevector values inside are the dimensionless angles x=a*kx,y=a*ky;
spatial integration d^2 k/(2pi)^2 has a factor 1/a^2. No broadening or
fitted normalization, no eigensolver, no SciPy dependency.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np


@lru_cache(maxsize=None)
def _nodes(order):
    if not isinstance(order,int) or order<4:raise ValueError('Gauss order >= 4 required')
    return np.polynomial.legendre.leggauss(order)


def _energy_coeff(y, *, a=1., beta=.25, v=.825, M=.1):
    """E^2(kx,ky) = c2 * cos(x)^2 + c1 * cos(x) + c0."""
    y=np.asarray(y,dtype=float)
    w=v/a
    B=2*M/a**2
    A=beta+B*(2-np.cos(y))
    c2=B*B-w*w
    c1=-2*A*B
    c0=w*w*(1+np.sin(y)**2)+A*A
    return c2,c1,c0,A,B,w


def _energy_at_branches(y, **params):
    c2,c1,c0,_,_,_=_energy_coeff(y,**params)
    # E at cos(x)=+1, -1 and at interior extremum of E(x) in [0,pi].
    up=np.sqrt(np.maximum(c2+c1+c0,0.))
    down=np.sqrt(np.maximum(c2-c1+c0,0.))
    critical=np.clip(-c1/(2*c2),-1.,1.)
    stationary=np.sqrt(np.maximum(c2*critical**2+c1*critical+c0,0.))
    return (up,down,stationary)


def _bisect(func,lo,hi,*,iterations=62):
    fl=func(lo)
    if fl*func(hi)>0:raise ValueError('unbracketed contour feature')
    for _ in range(iterations):
        mid=.5*(lo+hi)
        fm=func(mid)
        if fm==0:return mid
        if fl*fm<=0:hi=mid
        else:lo=mid;fl=fm
    return .5*(lo+hi)


def y_topology_splits(mu, *, scan=400, **params):
    """Locate all interior ky values where a Fermi contour touches/endpoints.

    A dense bracketing scan is used only to isolate topology events, never to
    evaluate the response. Checks are frozen against contour-root convergence.
    """
    threshold=abs(float(mu))
    if not np.isfinite(threshold):raise ValueError('finite chemical potential required')
    grid=np.linspace(0,np.pi,scan+1)
    vals=np.stack(_energy_at_branches(grid,**params))
    roots=[]
    for branch in range(3):
        diff=vals[branch]-threshold
        for j in range(scan):
            if diff[j]*diff[j+1]<0:
                r=_bisect(lambda yy: float(_energy_at_branches(yy,**params)[branch])-threshold,
                          float(grid[j]),float(grid[j+1]))
                if 1.e-9<r<np.pi-1.e-9:roots.append(float(r))
    return [0.,*sorted({round(x,13) for x in roots}),float(np.pi)]


def contour_roots_x(y,mu,**params):
    """Return all distinct kx roots on [0,pi], handles two crossings/ky."""
    c2,c1,c0,_,_,_=_energy_coeff(y,**params)
    D=c1*c1-4*c2*(c0-abs(float(mu))**2)
    if D<0:return []
    sq=np.sqrt(max(0.,float(D)))
    roots=[]
    for u in ((-c1-sq)/(2*c2),(-c1+sq)/(2*c2)):
        if -1+1e-12<float(u)<1-1e-12:
            x=float(np.arccos(np.clip(u,-1.,1.)))
            if all(abs(x-z)>1e-11 for z in roots):roots.append(x)
    return sorted(roots)


def _kernel_at(x,y,*,alpha='x',spin='x',mass_sign=1,**params):
    """Real d-vector formula evaluated without diagonalization."""
    from .d_vector_independent import dvector_dirac,pauli_trace_geometry
    a=params.get('a',1.)
    d,grad=dvector_dirac(np.asarray(x)/a,np.asarray(y)/a,mass_sign=mass_sign,**params)
    return pauli_trace_geometry(d,grad,alpha,spin)


def integrate_contour(mu, *, outer_order=96, inner_order=96, mass_sign=1,
                      photon_eV=.0041356676943898555,
                      a=1.,beta=.25,v=.825,M=.1):
    """Exact Fermi-boundary partition + tensor Gauss integral on smooth cells.

    Only C_xy and D_xx are returned; finite-frequency complex K also computed.
    Preserves full nonresonant denominator, not an arbitrary eta.
    """
    params=dict(a=a,beta=beta,v=v,M=M)
    bounds=y_topology_splits(mu,**params)
    n_outer,w_outer=_nodes(outer_order)
    n_inner,w_inner=_nodes(inner_order)
    c_total=0.;d_total=0.;K_real=0.;K_imag=0.
    max_roots=0; n_subdomains=0
    for ylo,yhi in zip(bounds[:-1],bounds[1:]):
        ys=.5*(ylo+yhi)+.5*(yhi-ylo)*n_outer
        for y,wy in zip(ys,w_outer):
            roots=contour_roots_x(float(y),mu,**params)
            max_roots=max(max_roots,len(roots))
            edges=[0.,*roots,np.pi]
            for xlo,xhi in zip(edges[:-1],edges[1:]):
                xmid=.5*(xlo+xhi)
                c2,c1,c0,*_=_energy_coeff(y,**params)
                E=np.sqrt(c2*np.cos(xmid)**2+c1*np.cos(xmid)+c0)
                if E<=abs(mu):continue
                xs=.5*(xlo+xhi)+.5*(xhi-xlo)*n_inner
                X,Y=np.broadcast_arrays(xs,np.array(y))
                # geometric Pauli formulas are independent from the spectral Kubo engine
                qC=_kernel_at(X,Y,alpha='x',spin='y',mass_sign=mass_sign,**params)
                qD=_kernel_at(X,Y,alpha='x',spin='x',mass_sign=mass_sign,**params)
                gap=qD['gap']
                if photon_eV>=float(np.min(gap)):
                    raise ValueError('resonant frequency not supported')
                wt=.25*(yhi-ylo)*(xhi-xlo)*wy*w_inner/(a*a*np.pi*np.pi)
                # 4 quarters of BZ because both response kernels are kx/ky even
                C=qC['Z01'];D=qD['D_per_omega01']
                den=1-(photon_eV/gap)**2
                c_total+=float(np.sum(wt*C))
                d_total+=float(np.sum(wt*D))
                K_real+=float(np.sum(wt*C/den))
                K_imag+=float(np.sum(wt*photon_eV*D/den))
                n_subdomains+=1
    return dict(mu_eV=mu,C_xy=c_total,D_peromega_xx=d_total,
                K_real_xy=K_real,K_imag_xx=K_imag,
                topology_splits_ky=bounds,max_x_crossings=max_roots,
                occupied_subdomains=n_subdomains,outer_order=outer_order,
                inner_order=inner_order)


def integrate_conduction_fermi_line(mu, *, order=160,a=1.,beta=.25,v=.825,M=.1):
    """Independent Eq.(13) Fermi line, no Fermi sea integrand evaluation.

    For epsilon_+ = E and epsilon_-=-E, each Fermi crossing contributes
       sign(partial_{kx} E) * sigma^y_+ = sign(partial_x E)*(-v/a sin x /E).
    This equals the occupied-sea result via integration by parts on torus.
    """
    threshold=abs(float(mu));params=dict(a=a,beta=beta,v=v,M=M)
    bounds=y_topology_splits(mu,**params)
    nodes,weights=_nodes(order)
    accum=0.;root_count=0
    for ylo,yhi in zip(bounds[:-1],bounds[1:]):
        ys=.5*(ylo+yhi)+.5*(yhi-ylo)*nodes
        for y,wy in zip(ys,weights):
            c2,c1,c0,_,_,w=_energy_coeff(y,**params)
            for x in contour_roots_x(float(y),mu,**params):
                u=np.cos(x)
                dE_dx=(2*c2*u+c1)*(-np.sin(x))/(2*threshold)
                sy=-w*np.sin(x)/threshold
                accum+=float(.5*(yhi-ylo)*wy*np.sign(dE_dx)*sy)
                root_count+=1
    return dict(mu_eV=mu,C_xy_eq13=float(accum/(a*a*np.pi*np.pi)),
                gauss_order=order,contour_roots_evaluated=root_count,
                topology_splits_ky=bounds)
