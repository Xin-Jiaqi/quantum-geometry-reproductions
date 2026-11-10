"""Exact endpoint product split in the same N/Delta^2 convention as KM.

Each factor is phase invariant for resolved bands. Band matching is a
controlled ambiguity near crossings; different matchings retain the exact sum.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment
adj=lambda a:a.conj().swapaxes(-1,-2)

def match(reference,field,overlap,nocc):
    ov=adj(field)@overlap@reference
    perms=[];powers=[]
    for o in ov:
        perm=np.zeros(o.shape[0],dtype=int)
        for sl in (slice(0,nocc),slice(nocc,len(o))):
            inds=np.arange(len(o))[sl]
            rows,cols=linear_sum_assignment(-abs(o[sl,sl])**2)
            perm[inds[cols]]=inds[rows]
        perms.append(perm);powers.append(abs(o[perm,np.arange(len(o))])**2)
    return np.array(perms),np.array(powers)

def reorder(e,u,perm):
    return np.take_along_axis(e,perm,axis=1),np.take_along_axis(u,perm[:,None,:],axis=2)

def split(np_,dp,nm,dm,field_interval):
    energy=(np_+nm)*(dp-dm)/(2*field_interval)
    matrix=(np_-nm)*(dp+dm)/(2*field_interval)
    total=(np_*dp-nm*dm)/field_interval
    return energy,matrix,total

def eig(d,nb,b):
    h=np.eye(nb)*d['energies'][:,:nb,None]+b*d['Z'][:,:nb,:nb]
    return np.linalg.eigh(h)

def factors(d,nb,e,u):
    nk=len(d['k']);mesh=round(nk**.5);nv=int(d['nocc'])
    occ=u[:,:,:nv];empty=u[:,:,nv:];deriv=np.zeros_like(occ,dtype=complex)
    links=d['links'][:,:,:nb,:nb];smallest=1.
    for a in range(4):
        v=links[:,a]@occ[d['neighbors'][:,a]]
        overlap=adj(occ)@v
        smallest=min(smallest,float(np.linalg.svd(overlap,compute_uv=False).min()))
        deriv+=(-1)**a*mesh*d['cell'][a//2,0]/(4*np.pi)*v@np.linalg.inv(overlap)
    rcv=(adj(empty)@(1j*deriv)).swapaxes(1,2)
    delta=e[:,:nv,None]-e[:,None,nv:]
    velocity=-1j*delta*rcv
    ops=[d['spin_total'][:,:nb,:nb],d['spin'][:,:nb,:nb]]
    ops.append(ops[0]-ops[1])
    factor=1e10*(2.00231930436256*5.7883818060e-5/2)/np.linalg.norm(np.cross(d['cell'][0],d['cell'][1]))
    n=np.array([-2*np.imag((adj(occ)@s@empty)*velocity)*factor for s in ops])
    den=delta**-2
    # Order-independent common occupied/empty gap after tracked permutations.
    gap=float(e[:,nv:].min()-e[:,:nv].max())
    assert gap>.05 and smallest>1e-5,(gap,smallest)
    spacing=min(float(np.diff(np.sort(e[:,:nv],axis=1),axis=1).min()),float(np.diff(np.sort(e[:,nv:],axis=1),axis=1).min()))
    return n,den,{'gap_eV':gap,'min_neighbor_singular':smallest,'min_same_sector_spacing_eV':spacing}
