"""v60: finite-open QWZ local Chern marker versus the SAME periodic Hamiltonian.

This validates occupied-projector topology, NOT equality of finite-device
Hall voltage and periodic bulk Hall conductivity, nor bulk orbital mixed EB.
"""
import numpy as np
from .periodic_kubo_v59 import sx,sy,sz

def qwz_open(Lx,Ly,mass):
    if not isinstance(Lx,int) or not isinstance(Ly,int) or min(Lx,Ly)<3 or not np.isfinite(mass):
        raise ValueError('Lx,Ly>=3; finite mass')
    n=2*Lx*Ly;H=np.zeros((n,n),complex);xy=[]
    Tx=(sz-1j*sx)/2;Ty=(sz-1j*sy)/2
    for y in range(Ly):
        for x in range(Lx):
            i=2*(y*Lx+x)
            H[i:i+2,i:i+2]=mass*sz
            xy.extend(((x,y),(x,y)))
            for xx,yy,hop in ((x+1,y,Tx),(x,y+1,Ty)):
                if xx<Lx and yy<Ly:
                    j=2*(yy*Lx+xx)
                    H[i:i+2,j:j+2]=hop
                    H[j:j+2,i:i+2]=hop.conj().T
    return H,np.asarray(xy,float)

def local_chern_marker(H,coords,mu=0):
    """Bianco-Resta -4π Im diag(P X (1-P) Y P), sum orbital markers by site."""
    H=np.asarray(H,complex);coords=np.asarray(coords,float);n=len(H)
    if H.shape!=(n,n) or coords.shape!=(n,2) or not np.allclose(H,H.conj().T,atol=1e-12):
        raise ValueError('matched Hermitian H and orbital coordinates')
    e,u=np.linalg.eigh(H)
    occ=e<mu
    P=u[:,occ]@u[:,occ].conj().T
    Q=np.eye(n)-P
    x=coords[:,0];y=coords[:,1]
    pxq=P@(x[:,None]*Q)
    yP=y[:,None]*P
    diagonal=np.einsum('ij,ji->i',pxq,yP)
    marker=-4*np.pi*np.imag(diagonal)
    return {'marker_orbital':marker,'marker_site':marker.reshape(-1,2).sum(axis=1),
            'sum_over_all_orbitals':float(marker.sum()),
            'occupation':int(occ.sum()),'min_distance_to_mu':float(np.min(abs(e-mu)))}

def compare_finite_periodic(L=11,mass=-1,margin=3,mu=0,kmesh=31):
    from .periodic_kubo_v59 import periodic_chern
    H,xy=qwz_open(L,L,mass)
    local=local_chern_marker(H,xy,mu)
    c=local['marker_site'].reshape(L,L)
    if margin<1 or 2*margin>=L:raise ValueError('nonempty central region required')
    bulk=periodic_chern(mass,kmesh)
    return {'L':L,'mass':mass,'margin':margin,
            'interior_marker_mean':float(np.mean(c[margin:L-margin,margin:L-margin])),
            'center_marker':float(c[L//2,L//2]),
            'full_finite_trace':local['sum_over_all_orbitals'],
            'periodic_chern':bulk['chern_wilson'],
            'periodic_spectral_chern':bulk['chern_spectral'],
            'finite_min_distance_mu':local['min_distance_to_mu'],
            'status':'SAME_QWZ_HAMILTONIAN_FINITE_PROJECTOR_VS_BULK_CHERN_ONLY_NOT_TRANSPORT'}
