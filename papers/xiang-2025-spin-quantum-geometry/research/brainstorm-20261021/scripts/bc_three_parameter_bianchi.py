"""Exploratory exact-geometry consistency check: Bianchi identity on (kx, sliding delta, spin-link twist theta).
This is ordinary Berry curvature mathematics; it is NOT a novel material prediction.
"""
import numpy as np
I=np.eye(2,dtype=complex);sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]],complex);sz=np.diag([1,-1]).astype(complex)
def X(a,b):return np.kron(a,b)
S=X(sz,sx)/2
Se,Ue=np.linalg.eigh(S)
def rot(theta):return (Ue*np.exp(-1j*theta*Se))@Ue.conj().T
def H(kx,ky,delta,theta):
    hp=(1+0.24*np.sin(delta))*X(sx,I)+0.21*np.cos(kx)*X(sy,I)+(0.33+0.16*np.cos(delta))*X(sy,sx)
    hl=0.49*np.sin(kx)*X(sz,I)+1.3*X(sz,sz)+0.24*np.sin(ky)*X(sz,sx)+0.11*np.cos(delta)*X(sz,sy)
    U=rot(theta)
    return U@hp@U.conj().T+hl

def proj(kx,ky,delta,theta):
    e,U=np.linalg.eigh(H(kx,ky,delta,theta))
    return U[:,:2]@U[:,:2].conj().T

def F(kx,ky,delta,theta,var,dh=1e-3):
    P=proj(kx,ky,delta,theta)
    Pk=(proj(kx+dh,ky,delta,theta)-proj(kx-dh,ky,delta,theta))/(2*dh)
    if var=='theta':
      Pvar=(proj(kx,ky,delta,theta+dh)-proj(kx,ky,delta,theta-dh))/(2*dh)
    elif var=='delta':
      Pvar=(proj(kx,ky,delta+dh,theta)-proj(kx,ky,delta-dh,theta))/(2*dh)
    return (1j*np.trace(P@(Pk@Pvar-Pvar@Pk))).real

def integral(delta,theta,var,N=29):
    ks=np.linspace(-np.pi,np.pi,N,endpoint=False)
    return np.mean([F(kx,ky,delta,theta,var) for kx in ks for ky in ks])

if __name__=='__main__':
    delta,theta=0.47,0.31
    h=0.01
    left=(integral(delta+h,theta,'theta')-integral(delta-h,theta,'theta'))/(2*h)
    right=(integral(delta,theta+h,'delta')-integral(delta,theta-h,'delta'))/(2*h)
    print('∂delta ∫F_k-theta:',left)
    print('∂theta ∫F_k-delta:',right)
    print('Bianchi residual:',left-right)
