"""Reproducible *synthetic* 4-band spin/orbital BZ example, no DFT inputs."""
from __future__ import annotations
import argparse,json
import numpy as np
from flux_lab.core import PAULI
from flux_lab.causal_response import reciprocal_residual

SPIN=np.stack([np.kron(np.eye(2),s) for s in PAULI])
TAU=np.stack([np.kron(s,np.eye(2)) for s in PAULI])

def synthetic_spin_orbital(kx,ky,tilt=.23):
    m=.28+.65*(2-np.cos(kx)-np.cos(ky))
    H=(np.sin(kx)*SPIN[0]+np.sin(ky)*SPIN[1]+m*SPIN[2]
       +2.5*TAU[2]+.7*(TAU[0]@SPIN[2])+tilt*np.sin(kx)*np.eye(4))
    V=np.stack([np.cos(kx)*SPIN[0]+.65*np.sin(kx)*SPIN[2]+tilt*np.cos(kx)*np.eye(4),
                np.cos(ky)*SPIN[1]+.65*np.sin(ky)*SPIN[2]])
    return H,V,SPIN

def bz_example(n=16,mu=2.3,temperature=.025,omega_drive=.035):
    if n<6:raise ValueError("kmesh must be >=6")
    k=-np.pi+2*np.pi*(np.arange(n)+.5)/n
    average=np.zeros((2,3),complex);max_res=0.
    for x in k:
        for y in k:
            H,V,S=synthetic_spin_orbital(x,y)
            q=reciprocal_residual(H,V,S,mu=mu,temperature=temperature,
                                  omega_drive=omega_drive)
            average+=q["delta_chi_JB"]
            max_res=max(max_res,q["max_residual"])
    average/=n*n
    return dict(kmesh=n,mu=mu,temperature=temperature,omega_drive=omega_drive,
                chi_yz_im=float(average[1,2].imag),
                chi_yz_re=float(average[1,2].real),
                max_local_reciprocal_residual=max_res,
                normalization="BZ d2k/(2pi)^2; natural units; interband only")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--kmesh",type=int,default=16)
    a=p.parse_args()
    print(json.dumps(bz_example(a.kmesh),indent=2))
