"""Construct a complete in-plane MMN finite-difference shell, without tensor averaging."""
from pathlib import Path
import argparse
import json
import os
import socket
import numpy as np


def complete_stencil(seed, cache_file, destination):
    with np.load(cache_file) as data:
        cache = {k:data[k] for k in data.files}
    k, mesh, cell = cache['k'], cache['mesh'], cache['cell']
    nk, nb = cache['energies'].shape
    entries = {}
    with (seed/'wannier90.mmn').open() as f:
        next(f)
        bands, npoints, nneighbors = map(int,next(f).split())
        if (bands,npoints)!=(nb,nk):raise ValueError('MMN dimensions')
        for _ in range(nk*nneighbors):
            i,j,*g=map(int,next(f).split());i-=1;j-=1
            raw=np.fromstring(' '.join(next(f) for _ in range(nb*nb)),sep=' ').reshape(-1,2)
            v=k[j]+g-k[i]
            if abs(v[2])>1e-7:continue
            integer=v[:2]*mesh
            if np.max(abs(integer-np.rint(integer)))>1e-7:raise ValueError('Non-grid link')
            step=tuple(map(int,np.rint(integer)))
            if step==(0,0):continue
            if step not in entries:entries[step]={}
            if i in entries[step]:raise ValueError('Duplicate MMN step')
            entries[step][i]=(j,(raw[:,0]+1j*raw[:,1]).reshape(nb,nb).T)
    steps=sorted(entries)
    if any(len(entries[s])!=nk for s in steps):raise ValueError('Incomplete directional shell')
    if any(tuple(-np.array(s)) not in entries for s in steps):raise ValueError('Missing reverse shell')
    b_reduced=np.column_stack([np.array(steps)/mesh,np.zeros(len(steps))])
    b_cart=(b_reduced @ (2*np.pi*np.linalg.inv(cell).T))[:,:2]
    norms=np.linalg.norm(b_cart,axis=1)
    groups=[]
    for norm in norms:
        if not any(np.isclose(norm,x,rtol=1e-7,atol=1e-10) for x in groups):groups.append(norm)
    if len(groups)>3:raise ValueError('More than3 radial shells needs explicit convergence/stencil design')
    group_indices=[np.flatnonzero(np.isclose(norms,x,rtol=1e-7,atol=1e-10)) for x in groups]
    moments=[]
    for ids in group_indices:
        b=b_cart[ids];moments.append([np.sum(b[:,0]**2),np.sum(b[:,1]**2),np.sum(b[:,0]*b[:,1])])
    matrix=np.array(moments).T
    shell_weights=np.linalg.lstsq(matrix,np.array([1.,1.,0.]),rcond=None)[0]
    residual=float(np.max(abs(matrix@shell_weights-[1.,1.,0.])))
    if residual>1e-8 or np.min(shell_weights)<=0:raise ValueError('No positive complete first-derivative shell')
    weights=np.zeros(len(steps))
    for ids,w in zip(group_indices,shell_weights):weights[ids]=w
    coeff=weights[:,None]*b_cart
    links=np.stack([np.stack([entries[s][i][1] for s in steps]) for i in range(nk)])
    neighbors=np.array([[entries[s][i][0] for s in steps] for i in range(nk)])
    reverse_error=0.
    for q,step in enumerate(steps):
        rev=steps.index(tuple(-np.array(step)))
        if not np.array_equal(neighbors[neighbors[:,q],rev],np.arange(nk)):raise ValueError('Shell neighbor reversal')
        reverse_error=max(reverse_error,float(np.max(abs(links[:,q]-links[neighbors[:,q],rev].conj().swapaxes(-1,-2)))))
    if reverse_error>1e-6:raise ValueError('Shell link adjoint mismatch')
    cache.update(links=links,neighbors=neighbors,stencil_steps=np.array(steps),derivative_weights=coeff)
    np.savez_compressed(destination,**cache)
    report={'stencil_steps':steps,'b_cart_Ainv':b_cart.tolist(),'shell_weights_A2':shell_weights.tolist(),'derivative_weights_A':coeff.tolist(),'first_moment_error':residual,'reverse_link_error':reverse_error,'source_cache':str(cache_file),'output':str(destination),'interpretation':'geometric derivative stencil; no projection/zeroing/averaging of response tensor'}
    destination.with_suffix('.json').write_text(json.dumps(report,indent=2))
    return report


if __name__=='__main__':
    if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute only')
    p=argparse.ArgumentParser();p.add_argument('seed',type=Path);p.add_argument('cache',type=Path);p.add_argument('destination',type=Path);a=p.parse_args()
    complete_stencil(a.seed,a.cache,a.destination)
