"""SM Table I: in-plane PHE rank-three tensor, axial magnetic index, 3D point ops."""
import numpy as np
from scipy.linalg import null_space

LABELS=("xx,x","yy,x","xy,x","xx,y","yy,y","xy,y")
PAIRS=((0,0,0),(1,1,0),(0,1,0),(0,0,1),(1,1,1),(0,1,1))
def rotation_z(theta):
    r=np.eye(3);r[:2,:2]=[[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]
    return r
def _tensor(c):
    t=np.zeros((2,2,2))
    for v,(i,j,k) in zip(c,PAIRS):t[i,j,k]=t[j,i,k]=v
    return t
def representation(r,antiunitary=False):
    r=np.asarray(r,float)
    if r.shape!=(3,3) or not np.allclose(r@r.T,np.eye(3),atol=1e-12):
        raise ValueError("3D orthogonal matrix needed")
    if not np.allclose(r[2,:2],0) or not np.allclose(r[:2,2],0):
        raise ValueError("2D plane must be preserved")
    xy=r[:2,:2]; factor=np.linalg.det(r)*(-1 if antiunitary else 1)
    def apply(c):
        transformed=factor*np.einsum("ia,jb,kc,abc->ijk",xy,xy,xy,_tensor(c))
        return [transformed[i,j,k] for i,j,k in PAIRS]
    return np.column_stack([apply(v) for v in np.eye(6)])
def invariant_basis(ops):
    mats=[representation(r,anti)-np.eye(6) for r,anti in ops]
    return null_space(np.concatenate(mats) if mats else np.zeros((0,6)),rcond=1e-12)
def _table():
    xz=np.diag([1.,-1.,1.]); yz=np.diag([-1.,1.,1.])
    c2y=np.diag([-1.,1.,-1.]);c3=rotation_z(2*np.pi/3);c6=rotation_z(np.pi/3)
    return {
      "1":[],"2":[(c2y,False)],"m":[(xz,False)],
      "m'm2'":[(xz,False),(yz,True)],"2'":[(c2y,True)],"m'":[(xz,True)],
      "3":[(c3,False)],"6'":[(c6,True)],"3m":[(c3,False),(yz,False)],
      "6'mm'":[(c6,True),(yz,False)],"3m'":[(c3,False),(yz,True)]
    }
def table_i_acceptance():
    A=np.array([1.,-1.,0,0,0,-1.]);B=np.array([0,0,1.,1.,-1.,0])
    eye=np.eye(6)
    shapes={
      "1":eye,"2":eye[:,[2,3,4]],"m":eye[:,[2,3,4]],
      "m'm2'":eye[:,[2,3,4]],"2'":eye[:,[0,1,5]],"m'":eye[:,[0,1,5]],
      "3":np.column_stack([A,B]),"6'":np.column_stack([A,B]),
      "3m":A[:,None],"6'mm'":A[:,None],"3m'":B[:,None]
    }
    results={}
    for name,ops in _table().items():
        basis=invariant_basis(ops);reference=shapes[name]
        err=float(np.max(abs(basis@basis.T-reference@np.linalg.pinv(reference))))
        results[name]=dict(dim=int(basis.shape[1]),expected=int(reference.shape[1]),
                           error=err,passed=bool(basis.shape[1]==reference.shape[1] and err<5e-12))
    return results
