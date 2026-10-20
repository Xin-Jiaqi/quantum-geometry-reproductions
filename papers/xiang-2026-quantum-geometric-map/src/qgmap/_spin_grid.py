"""Numerical helpers preserved from the v30 offline multiband method, e=ℏ=1."""
import numpy as np
I=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1,-1]).astype(complex)

def mesh(N=61,offset=(.31,.23)):
    if N<4: raise ValueError("N>=4 required")
    x=(np.arange(N)+offset[0])*2*np.pi/N-np.pi
    y=(np.arange(N)+offset[1])*2*np.pi/N-np.pi
    return np.meshgrid(x,y,indexing="ij")

def fermi(e,mu,T):
    if T<=0: raise ValueError("Finite energy-temperature T>0 required")
    z=np.clip((e-mu)/T,-500,500)
    f=1/(1+np.exp(z))
    return f,-f*(1-f)/T
