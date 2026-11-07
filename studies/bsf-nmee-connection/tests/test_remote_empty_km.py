import numpy as np

from bsf_nmee.remote_empty_km import ingredients, assemble, selective_extra, flatten
from bsf_nmee.mixing_channels import adj


def test_remote_orbital_preserves_time_reversal_and_velocity():
    k = np.array([[.13, .29], [8.1, -3.2], [-5.8, 12.2]])
    direction = np.array([.5, np.sqrt(3)/2])
    tr = np.kron(np.eye(6), [[0., 1.], [-1., 0.]])
    for site in [0, 1]:
        for channel in ['all', 'VC', {'VC': 1., 'CC': .2}]:
            h, v = assemble(ingredients(k, site=site), .25, .01, .12, channel)
            hm, _ = assemble(ingredients(-k, site=site), .25, .01, .12, channel)
            np.testing.assert_allclose(h, tr@hm.conj()@adj(tr), atol=1e-12)
            hp, _ = assemble(ingredients(k+1e-5*direction, site=site), .25, .01, .12, channel)
            hn, _ = assemble(ingredients(k-1e-5*direction, site=site), .25, .01, .12, channel)
            np.testing.assert_allclose(v, (hp-hn)/2e-5, atol=1e-9)


def test_extra_orbital_zero_coupling_recovers_km_block():
    from bsf_nmee.bilayer_audit import AuditBilayer, AuditSpec
    from bsf_nmee.bilayer_protocol import ab_ba_shift
    k = np.array([[.13, .29], [8.1, -3.2]])
    model = AuditBilayer(AuditSpec(shift=ab_ba_shift(.25), staggered_eV=.045,
                                  tunneling_eV=.01))
    h, _ = assemble(ingredients(k), .25, .01, 0.)
    ix = np.r_[0:4, 6:10]
    np.testing.assert_allclose(h[:, ix[:, None], ix], model.matrices(k), atol=1e-14)
    assert np.max(abs(h[:, [4, 5, 10, 11]][:, :, ix])) == 0


def test_selective_coupling_and_flattened_velocities_and_mirror():
    from bsf_nmee.bilayer_audit import AuditBilayer
    k = np.array([[.13, .29], [8.1, -3.2], [-5.8, 12.2]])
    shift = np.array([.23, .41]); direction = np.array([.5, np.sqrt(3)/2])
    mirror = np.kron([[0., 1.], [1., 0.]], np.kron(np.eye(3), np.diag([1., -1.])))
    def at(points, displacement=shift):
        return selective_extra(ingredients(points, shift=displacement), .12, .01, .12, .1)
    h, v = at(k); hm, _ = at(k, -shift)
    np.testing.assert_allclose(hm, mirror@h@mirror, atol=1e-12)
    theta=2*np.pi/3
    rot=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    lattice=AuditBilayer().a
    spinrot=np.kron(np.eye(6),np.diag(np.exp(-.5j*theta*np.array([1.,-1.]))))
    hr,_=at(k@rot.T,shift@lattice@rot.T@np.linalg.inv(lattice))
    np.testing.assert_allclose(hr,spinrot@h@adj(spinrot),atol=1e-12)
    hf, vf = flatten(h, v)
    np.testing.assert_allclose(np.linalg.eigvalsh(hf), np.tile([-.015]*4+[.015]*8,(len(k),1)), atol=1e-12)
    hp, vp = at(k+1e-5*direction); hn, vn = at(k-1e-5*direction)
    np.testing.assert_allclose(v, (hp-hn)/2e-5, atol=1e-9)
    np.testing.assert_allclose(vf, (flatten(hp,vp)[0]-flatten(hn,vn)[0])/2e-5, atol=1e-9)
