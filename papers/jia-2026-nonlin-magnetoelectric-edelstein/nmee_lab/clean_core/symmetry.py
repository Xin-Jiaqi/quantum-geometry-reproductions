"""Spatial symmetry bookkeeping for intrinsic NMEE, not a full magnetic-group engine."""
import numpy as np

def transform_nmee_tensor(gamma:np.ndarray,R:np.ndarray):
    """Transform Gamma[spin_alpha, Zeeman_gamma, E_beta] under spatial R.

    Spin and B are both axial (det R) R, E is polar R, so two det factors cancel.
    Does not apply an antiunitary/time-reversal rule or interchange layers.
    """
    tensor=np.asarray(gamma,float);rot=np.asarray(R,float)
    if tensor.shape!=(3,3,3) or rot.shape!=(3,3):raise ValueError('3x3x3 tensor and 3x3 rotation required')
    if not np.allclose(rot@rot.T,np.eye(3),atol=1e-10):raise ValueError('R not orthogonal')
    return np.einsum('ia,jb,kc,abc->ijk',rot,rot,rot,tensor,optimize=True)

def symmetry_residual(tensor,R):
    transformed=transform_nmee_tensor(tensor,R)
    return float(np.max(abs(np.asarray(tensor)-transformed)))
