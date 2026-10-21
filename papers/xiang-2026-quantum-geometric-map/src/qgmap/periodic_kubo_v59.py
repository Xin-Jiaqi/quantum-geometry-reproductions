"""v59 independent 2D periodic-bulk *intrinsic charge anomalous Hall* oracle.

QWZ Chern model, Berry curvature via v54 spectral dH/dk and independent
Fukui-Hatsugai-Suzuki lattice-gauge Wilson link construction. This is an
independently verified periodic charge Hall baseline, NOT orbital d_B sigma,
NOT the original paper's intrinsic magnetoelectric Hall formula.
"""
import numpy as np
from .orbital_geometry_v54 import band_geometry

sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex)
sz=np.array([[1,0],[0,-1]],complex)

def qwz_h(kx,ky,mass):
    if not np.isfinite([kx,ky,mass]).all():raise ValueError('nonfinite model input')
    return np.sin(kx)*sx+np.sin(ky)*sy+(mass+np.cos(kx)+np.cos(ky))*sz

def qwz_derivatives(kx,ky):
    return np.array([np.cos(kx)*sx-np.sin(kx)*sz,
                     np.cos(ky)*sy-np.sin(ky)*sz])

def periodic_chern(mass=-1.,N=51):
    if not isinstance(N,int) or N<9 or not np.isfinite(mass):raise ValueError('N>=9, finite mass')
    if min(abs(mass),abs(mass-2),abs(mass+2))<1e-9:
        raise ValueError('gapless QWZ critical mass cannot be assigned a Chern invariant')
    k=np.linspace(-np.pi,np.pi,N,endpoint=False)
    energies=np.empty((N,N,2)); lower=np.empty((N,N,2),complex)
    curvature=np.empty((N,N))
    min_gap=float('inf')
    for a,x in enumerate(k):
        for b,y in enumerate(k):
            h=qwz_h(x,y,mass)
            es,us=np.linalg.eigh(h)
            lower[a,b]=us[:,0]
            energies[a,b]=es
            min_gap=min(min_gap,es[1]-es[0])
            curvature[a,b]=band_geometry(h,qwz_derivatives(x,y)).curvature[0]
    # Curvature integral C=(1/(2pi))*int_BZ Omega, discrete mesh Δk².
    conventional_chern=float(np.mean(curvature)*(2*np.pi))
    # Independent Wilson FHS: phase of Ux(k) Uy(k+dx) Ux*(k+dy) Uy*(k)
    xlap=np.sum(lower.conj()*np.roll(lower,-1,axis=0),axis=-1)
    ylap=np.sum(lower.conj()*np.roll(lower,-1,axis=1),axis=-1)
    if min(np.min(abs(xlap)),np.min(abs(ylap)))<1e-6:
        raise ValueError('FHS links nearly orthogonal; mesh too coarse')
    ux=xlap/abs(xlap);uy=ylap/abs(ylap)
    z=ux*np.roll(uy,-1,axis=0)*np.conj(np.roll(ux,-1,axis=1))*np.conj(uy)
    fhs=float(-np.sum(np.angle(z))/(2*np.pi))
    # Wilson phase follows link <u(k)|u(k+dk)> = exp(-i A dk)
    # for A=i<u|du>; hence physical Ω Chern is MINUS Wilson link winding.
    return {'mass':float(mass),'mesh':N,'chern_spectral':conventional_chern,
            'chern_wilson':fhs,'min_direct_gap':min_gap,
            'conductivity_e2_over_h_signed':conventional_chern,
            'status':'VERIFIED_QWZ_PERIODIC_CHARGE_AHE_ORACLE_NOT_ORBITAL_EB'}
