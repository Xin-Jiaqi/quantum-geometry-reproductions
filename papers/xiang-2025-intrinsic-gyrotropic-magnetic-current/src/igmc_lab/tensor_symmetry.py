"""Independent small magnetic-point-group Neumann checks for IGMC rank-2 tensors.

C transforms as a time-even polar-axial cross tensor (J/B),
D as time-odd polar-axial cross tensor (J/dB/dt) in clean nondissipative AC.
Printed PRL Eq.(14) is retained as an explicit *alternative* protocol.
"""
import numpy as np

IDENTITY = np.eye(3)
C4Z = np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
M_XZ = np.diag([1.,-1.,1.])
INVERSION = -np.eye(3)

def transform(tensor, rotation, *, antiunitary=False, channel='C', protocol='channel_resolved_symmetry'):
    """Return eta*det(R)*R*T*R^t for a polar-axial rank-2 tensor."""
    t=np.asarray(tensor,dtype=float);r=np.asarray(rotation,dtype=float)
    if t.shape!=(3,3) or r.shape!=(3,3) or not np.allclose(r@r.T,IDENTITY,atol=1e-12):
        raise ValueError('tensor 3x3 and orthogonal spatial rotation required')
    if channel not in ('C','D'):raise ValueError('channel must be C or D')
    if protocol=='channel_resolved_symmetry':
        eta=(-1 if channel=='D' else +1) if antiunitary else 1
    elif protocol=='paper_printed_eq14':
        eta=-1 if antiunitary else 1
    else:raise ValueError('explicit symmetry protocol required')
    return eta*round(np.linalg.det(r))*r@t@r.T

def allowed_basis_dim(generators, *, channel='C', protocol='channel_resolved_symmetry'):
    """Number of independent real tensor components, computed by linear nullspace."""
    constraints=[]
    eye=np.eye(9)
    for r,anti in generators:
        mapping=np.column_stack([transform(eye[:,j].reshape(3,3),r,
            antiunitary=anti,channel=channel,protocol=protocol).reshape(9)
                                 for j in range(9)])
        constraints.append(eye-mapping)
    if not constraints:return 9
    s=np.linalg.svd(np.concatenate(constraints),compute_uv=False)
    return int(9-np.count_nonzero(s>1e-9))

def example_generators():
    return {
      '4mm1_prime':[(C4Z,False),(M_XZ,False),(IDENTITY,True)],
      '4m_prime_m_prime':[(C4Z,False),(M_XZ,True)],
      'inversion':[(INVERSION,False)],
    }
