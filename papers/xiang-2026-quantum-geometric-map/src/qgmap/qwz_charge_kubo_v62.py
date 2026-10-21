"""v62 independent periodic charge Kubo anomalous Hall at specified mu and T.

Two exact-intrinsic-band oracles: analytic d-vector scalar Berry curvature,
and Kubo interband velocity eigenbasis (v54). These address the *same* bulk
charge sigma_xy observable at finite Fermi filling, not finite-device Rxy.
Metallic dissipative sigma_xx and generic rho_xy require a separate bath model.
"""
import numpy as np
from scipy.special import expit
from .periodic_kubo_v59 import qwz_h,qwz_derivatives
from .orbital_geometry_v54 import band_geometry


def qwz_charge_kubo(mass=-1.,mu=0.,temperature=0.,N=45):
    if not isinstance(N,int) or N<9 or not np.isfinite([mass,mu,temperature]).all() or temperature<0:
        raise ValueError('finite inputs; T>=0 and N>=9')
    if min(abs(mass),abs(mass-2),abs(mass+2))<1e-8:
        raise ValueError('mass at gapless critical point')
    ka=np.linspace(-np.pi,np.pi,N,endpoint=False)
    via_pauli=[];via_kubo=[];minimal=float('inf')
    for x in ka:
        for y in ka:
            dx=np.array([np.sin(x),np.sin(y),mass+np.cos(x)+np.cos(y)])
            u=np.array([np.cos(x),0.,-np.sin(x)])
            v=np.array([0.,np.cos(y),-np.sin(y)])
            r=np.linalg.norm(dx)
            e=np.array([-r,r]);minimal=min(minimal,2*r)
            f=expit(-(e-mu)/temperature) if temperature else (e<mu).astype(float)+.5*(e==mu)
            weight=float(f[0]-f[1])
            omega=+float(np.dot(dx,np.cross(u,v)))/(2*r**3)
            via_pauli.append(omega*weight)
            geometry=band_geometry(qwz_h(x,y,mass),qwz_derivatives(x,y))
            via_kubo.append(float(np.dot(f,geometry.curvature)))
    return {'mass':float(mass),'mu':float(mu),'temperature':float(temperature),'N':N,
            'sigma_pauli_q2_over_h':float(2*np.pi*np.mean(via_pauli)),
            'sigma_interband_q2_over_h':float(2*np.pi*np.mean(via_kubo)),
            'min_direct_gap':float(minimal),
            'status':'PERIODIC_INTRINSIC_CHARGE_KUBO_NO_SIGMA_XX_NO_ORBITAL_EB'}
