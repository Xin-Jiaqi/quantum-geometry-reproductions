import numpy as np
from bsf_nmee.mixing_channels import split_channels, adj


def example():
    rng = np.random.default_rng(918)
    a = rng.normal(size=(3,8,8)) + 1j*rng.normal(size=(3,8,8))
    h = np.diag([-3,-2,2,3]*2) + .06*(a+adj(a))
    v = .03*(a+adj(a))
    return h,v


def test_subspace_channels_reconstruct_and_are_gauge_covariant():
    h,v=example();layers=[np.arange(4),np.arange(4,8)]
    intra,parts,p,diag=split_channels(h,layers,2)
    assert np.max(abs(h-intra-sum(parts.values())))<1e-12
    assert np.max(abs(p@parts['VC']@p))<1e-12
    rng=np.random.default_rng(514);w=np.zeros((8,8),complex)
    for ix in layers:
        q,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
        w[np.ix_(ix,ix)]=q
    _,rot,_,_=split_channels(adj(w)@h@w,layers,2)
    for key in parts:assert np.max(abs(rot[key]-adj(w)@parts[key]@w))<1e-12


def test_modified_velocity_includes_projector_derivative():
    h,v=example();layers=[np.arange(4),np.arange(4,8)]
    _,_,_,_,_,derivs=split_channels(h,layers,2,velocity=v)
    plus=split_channels(h+1e-5*v,layers,2)[1]
    minus=split_channels(h-1e-5*v,layers,2)[1]
    for key in derivs:assert np.max(abs(derivs[key]-(plus[key]-minus[key])/2e-5))<1e-9
