"""Exact global PT^2=-1 synthetic six-band three-doublet Hamiltonian.

NumPy-only version of the offline V19 model. NOT MnBi2Te4.
"""
import numpy as np

def _doublet_J(ncluster):
    return np.kron(np.eye(ncluster), np.array([[0,1],[-1,0]],complex))

class PTThreeDoublets:
    def __init__(self,seed=202519):
        rng=np.random.default_rng(seed)
        N=6
        A=rng.normal(size=(N,N))+1j*rng.normal(size=(N,N))
        G=(A+A.conj().T)/2
        self.J=_doublet_J(3)
        self.G=(G-self.J@G.conj()@self.J.conj().T)/2
        self.G=self.G/np.linalg.norm(self.G,ord=2)*.62
        self.g,self.V=np.linalg.eigh(self.G)
        self.centers=np.linspace(-2.25,2.45,3)
        self.slopes=np.linspace(-.20,.23,3)
        self.curvature=np.linspace(.08,.13,3)

    def __call__(self,kx,ky=0.):
        k=float(kx)
        U=(self.V*np.exp(-1j*np.sin(k)*self.g)[None,:])@self.V.conj().T
        e=self.centers+self.slopes*np.sin(k)+self.curvature*np.cos(k)
        D=np.diag(np.repeat(e,2))
        return U@D@U.conj().T

    def dH(self,kx):
        k=float(kx);H=self(k)
        dp=self.slopes*np.cos(k)-self.curvature*np.sin(k)
        U=(self.V*np.exp(-1j*np.sin(k)*self.g)[None,:])@self.V.conj().T
        Hprime=U@np.diag(np.repeat(dp,2))@U.conj().T
        return Hprime-1j*np.cos(k)*(self.G@H-H@self.G)

    def check_PT(self,kx):
        H=self(kx)
        return float(np.max(abs(self.J@H.conj()@self.J.conj().T-H)))
