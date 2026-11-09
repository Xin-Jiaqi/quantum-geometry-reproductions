import numpy as np
from bsf_nmee.material_subspaces import active_layer_basis, orthonormalize, adj


def test_selected_high_sector_preserves_valence_coupling_without_response_fit():
    rng=np.random.default_rng(81)
    layers=[np.arange(10),np.arange(10,20)]
    e=np.r_[np.arange(6)-8,np.arange(4)+1.]
    h=np.tile(np.diag(np.r_[e,e]),(3,1,1)).astype(complex)
    t=.05*(rng.normal(size=(3,10,10))+1j*rng.normal(size=(3,10,10)))
    h[:,:10,10:]=t;h[:,10:,:10]=adj(t)
    b,diag=active_layer_basis(h,layers,occupied=6)
    np.testing.assert_allclose(adj(b)@b,np.broadcast_to(np.eye(12),(3,12,12)),atol=1e-12)
    assert max(d['discarded_coupling_fraction'] for d in diag)<1e-12
    # Independent rotations of the original layer basis change no chosen subspace.
    rotation=np.zeros((20,20),complex)
    for ids in layers:
        u,_=np.linalg.qr(rng.normal(size=(10,10))+1j*rng.normal(size=(10,10)))
        rotation[np.ix_(ids,ids)]=u
    br,_=active_layer_basis(adj(rotation)@h@rotation,layers,occupied=6)
    np.testing.assert_allclose(br@adj(br),adj(rotation)@(b@adj(b))@rotation,atol=1e-11)


def test_truncated_basis_metric_is_restored():
    rng=np.random.default_rng(11)
    q=rng.normal(size=(2,20,12))+1j*rng.normal(size=(2,20,12))
    q,min_metric=orthonormalize(q)
    assert min_metric>0
    np.testing.assert_allclose(adj(q)@q,np.broadcast_to(np.eye(12),(2,12,12)),atol=1e-12)
