"""Projective orbital TB. Energies in eV, coordinates in Å, Pauli spin."""
from pathlib import Path
import numpy as np
from scipy.linalg import null_space

def polar_cols(a):
    u,s,vh=np.linalg.svd(a,full_matrices=False)
    return u@vh,s

def localized_model(data, frozen_empty=8, tail_scale=4.0):
    """Freeze the 28 d+p occupied and first empty states. Downfold the remaining empty space."""
    mesh=int(round(len(data['k'])**.5))
    hgrid=np.zeros((mesh,mesh,44,44),complex)
    singular=[]
    for k,coeff,energy in zip(data['k'],data['coeff'],data['energy']):
        nf=28+frozen_empty
        frozen,s=polar_cols(coeff[:,:nf].astype(complex))
        rest=null_space(frozen.conj().T)
        overlaps=rest.conj().T@coeff[:,nf:]
        weights=1/(1+(energy[nf:]-energy[nf])**2/tail_scale**2)
        lift,_=polar_cols((overlaps*weights).conj().T)
        effective=rest@lift.conj().T
        h=(frozen*energy[:nf])@frozen.conj().T+(effective*energy[nf:])@effective.conj().T
        i,j=np.rint(k[:2]*mesh).astype(int)%mesh
        hgrid[i,j]=h
        singular.append(s.min())
    # Remove the TR-odd part of the small B-prepared seed before applying ±b.
    tr=np.kron(np.eye(22),np.array([[0,1],[-1,0]]))
    neg=(-np.arange(mesh))%mesh
    partner=tr@hgrid[neg][:,neg].conj()@tr.T
    return (hgrid+partner)/2, float(min(singular))

def interpolate(data,hgrid,n):
    """Periodic H(R) in the shortest image of each orbital pair. v includes orbital embedding."""
    cell=data['cell'][:2,:2]
    tau=np.repeat(data['pos'][:,:2],[10,10,6,6,6,6],axis=0)
    diff=tau[None,:,:]-tau[:,None,:]
    mesh=hgrid.shape[0]
    hr=np.fft.fft2(hgrid,axes=(0,1))/mesh**2
    grid=np.zeros((n,n,44,44),complex)
    vx=np.zeros_like(grid)
    shifts=np.array([(i,j) for i in (-mesh,0,mesh) for j in (-mesh,0,mesh)])
    for i,j in np.ndindex(mesh,mesh):
        images=shifts+[i,j]
        dist=(images[:,None,None,:]+diff)@cell
        norm=(dist**2).sum(-1)
        closest=np.isclose(norm,norm.min(0),atol=1e-7,rtol=0)
        count=closest.sum(0)
        for (a,b),mask,dr in zip(images,closest,dist):
            if mask.any():
                weight=hr[i,j]*mask/count
                grid[a%n,b%n]+=weight
                vx[a%n,b%n]+=1j*weight*dr[:,:,0]
    h=np.fft.ifft2(grid,axes=(0,1))*n*n
    v=np.fft.ifft2(vx,axes=(0,1))*n*n
    assert np.max(abs(h-h.conj().swapaxes(-1,-2)))<1e-5
    return h.reshape(-1,44,44),v.reshape(-1,44,44)

NAMES=['dxy','dyz','dz2','dxz','dx2-y2']*2+['py','pz','px']*4
BOTTOM=np.repeat([True,False,True,False,True,False],[10,10,6,6,6,6])
def response(h,v,hop,offset,b=.002,keep=None,filled=28,ops=None):
    """Static length-gauge electric susceptibility differentiated in +b sigma_z convention."""
    inds=np.arange(44) if keep is None else np.repeat(np.isin(NAMES,keep),2).nonzero()[0]
    lower=BOTTOM[inds]
    total=np.zeros(2)
    edges=[-np.inf,np.inf]
    scale=np.where(lower[:,None]==lower[None,:],1.,hop)
    sx=np.kron(np.diag(lower[::2]),[[0,1],[1,0]])
    sz=np.tile([1.,-1.],len(inds)//2)
    for sl in range(0,len(h),64):
        mat=h[sl:sl+64][:,inds][:,:,inds]*scale
        vel=v[sl:sl+64][:,inds][:,:,inds]*scale
        pos=np.diag((~lower).astype(float)-.5) if ops is None else .5*np.eye(len(inds))-ops[0][sl:sl+64][:,inds][:,:,inds]
        spin_x=sx if ops is None else ops[1][sl:sl+64][:,inds][:,:,inds]
        spin_z=np.diag(sz) if ops is None else ops[2][sl:sl+64][:,inds][:,:,inds]
        for sign in (1,-1):
            energy,u=np.linalg.eigh(mat+offset*pos+sign*b*spin_z)
            uv=u[:,:,:filled];uc=u[:,:,filled:]
            vcv=uc.conj().transpose(0,2,1)@vel@uv
            svc=uv.conj().transpose(0,2,1)@spin_x@uc
            delta=energy[:,:filled,None]-energy[:,None,filled:]
            chi=-2*(svc*vcv.swapaxes(1,2)).imag/delta**2
            total[0 if sign==1 else 1]+=chi.sum()
            edges=[max(edges[0],energy[:,filled-1].max()),min(edges[1],energy[:,filled].min())]
    gap=edges[1]-edges[0]
    return {'hopping':float(hop),'detuning_eV':float(offset),'response':float((total[0]-total[1])/(2*b*len(h))*6.331927853335898e4) if gap>0 else None,'gap_eV':float(gap)}
