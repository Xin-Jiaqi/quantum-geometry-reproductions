"""Symmetry preflight for Gamma_{i,a;j}=d_(B_a)d_(E_j)<S_i>.

S_i and B_a are axial, E_j is polar. The two determinant factors cancel,
so the three spatial indices transform like one rank-3 polar tensor.

Time reversal alone is NOT a no-go for this *B-derivative* response:
spin and B are both TR-odd, electric field is TR-even.

Use operations of the SAME fixed stacking. Do not include a domain-mapping
mirror as a little-group symmetry of a displaced configuration.
"""
import numpy as np


def spatial_transform(T,R):
    R=np.asarray(R,float); T=np.asarray(T,float)
    if R.shape!=(3,3) or not np.allclose(R@R.T,np.eye(3),atol=1e-10):
        raise ValueError("R must be an orthogonal Cartesian 3x3 operation")
    if T.shape!=(3,3,3):
        raise ValueError("Tensor Gamma must be 3x3x3")
    return np.einsum('ia,jb,kc,abc->ijk',R,R,R,T,optimize=True)


def component_allowed(operations,indices=(0,2,0)):
    """Reynolds project a source tensor over a complete spatial point group."""
    if not operations:
        raise ValueError("Nonempty complete group required")
    t=np.zeros((3,3,3));t[indices]=1.
    projection=sum(spatial_transform(t,R) for R in operations)/len(operations)
    return bool(np.linalg.norm(projection)>1e-10)


def examples():
    I=np.eye(3);mz=np.diag([1,1,-1]);inv=-I
    theta=2*np.pi/3
    c=np.array([[np.cos(theta),-np.sin(theta),0],
                [np.sin(theta),np.cos(theta),0],[0,0,1]])
    return {
       'C1':component_allowed([I]),
       'C3':component_allowed([I,c,c@c]),
       'horizontal_mirror_mz':component_allowed([I,mz]),
       'inversion':component_allowed([I,inv]),
       'C3_plus_mz':component_allowed([I,c,c@c,mz,c@mz,c@c@mz]),
       'TR_alone_forbids':False,
       'note':'A mirror mapping two displaced domains is NOT a symmetry of one fixed displaced domain.'
    }
