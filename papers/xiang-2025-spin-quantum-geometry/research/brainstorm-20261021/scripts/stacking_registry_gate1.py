"""Gate 1: explicit interlayer registry, PT/C2z sewing and spin-transport falsification.

Independent synthetic model, NOT Xiang et al. original code or a ferroelectric material.
Run: python scripts/stacking_registry_gate1.py
Requires: numpy only. All output is a model-dependent RTA proxy, not a full Kubo result.
"""
import numpy as np

I = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]], complex)
sz = np.diag([1, -1]).astype(complex)
K = np.kron
I4 = np.eye(4, dtype=complex)
tlz, spin_y = K(sz, I), K(I, sy)
PT, C2 = K(sx, 1j * sy), K(I, 1j * sz)
tau_x, tau_y = K(sx, I), K(sy, I)

def hopping(kx, ky, dx, dy, amp=0.85, ell=0.8, Rmax=3):
    """All real-space bonds R+delta with Gaussian distance-dependent amplitudes."""
    value = deriv = 0j
    for nx in range(-Rmax, Rmax + 1):
        for ny in range(-Rmax, Rmax + 1):
            rx, ry = nx + dx, ny + dy
            weight = amp * np.exp(-(rx * rx + ry * ry) / (2 * ell * ell))
            term = weight * np.exp(1j * (kx * rx + ky * ry))
            value += term
            deriv += 1j * rx * term
    return value, deriv

def bloch(kx, ky, dx, dy, chiral=False, tunnel_amp=0.85, gamma_tied_to_t=False):
    """Rashba + staggered exchange, with one *optional assumed* odd-dispersion channel.

    gamma(delta) sin(kx) represents a separately postulated imaginary,
    layer-even intralayer hopping. NOT derived from the interlayer Gaussian hops.
    Setting gamma_tied_to_t=True imposes proportionality to tunnel_amp**2 by hand,
    which is a sensitivity study, NOT a microscopic perturbation derivation.
    """
    t, tx = hopping(kx, ky, dx, dy, amp=tunnel_amp)
    pref = 0.65 * (tunnel_amp / 0.85)**2 if gamma_tied_to_t else 0.65
    gamma = pref * np.sin(2 * np.pi * dx) if chiral else 0.0
    eps = -0.42 * (np.cos(kx) + 0.8 * np.cos(ky)) + gamma * np.sin(kx)
    veps = 0.42 * np.sin(kx) + gamma * np.cos(kx)
    alpha, exchange = 0.45, 0.8
    local = (alpha * (np.sin(ky)*K(sz, sx) - np.sin(kx)*K(sz, sy))
             + exchange * K(sz, sz))
    dlocal = -alpha * np.cos(kx) * K(sz, sy)
    inter = t.real * tau_x - t.imag * tau_y
    dinter = tx.real * tau_x - tx.imag * tau_y
    H = eps*I4 + inter + local
    v = veps*I4 + dinter + dlocal
    j = (v @ spin_y + spin_y @ v) / 4
    jstag = (j @ tlz + tlz @ j) / 2
    return H, v, j, jstag, inter, local

def fermi_surface(dx, dy, mu=1.3, T=0.13, N=41, chiral=False,
                  tunnel_amp=0.85, gamma_tied_to_t=False):
    """Diagonal occupation (constant-tau) proxy, omitting overall e*tau factors.

    The PT Kramers doublets are traced as projectors. This is NOT a full
    AC/intrinsic spin current, and is zero in the T=0 band insulator.
    """
    coef = np.zeros(2)
    vmax, cmin = -np.inf, np.inf
    for kx in np.linspace(-np.pi, np.pi, N, endpoint=False):
        for ky in np.linspace(-np.pi, np.pi, N, endpoint=False):
            H,v,j,js,*_ = bloch(kx,ky,dx,dy,chiral,tunnel_amp,gamma_tied_to_t)
            eigenvalues,U = np.linalg.eigh(H)
            for b in (0,2):
                P = U[:,b:b+2] @ U[:,b:b+2].conj().T
                energy = eigenvalues[b:b+2].mean()
                velocity = np.trace(P@v).real / 2
                if b == 0: vmax = max(vmax, energy)
                else: cmin = min(cmin, energy)
                fprime = 1 / (4*T*np.cosh(np.clip((energy-mu)/(2*T),-80,80))**2)
                coef += fprime * velocity * np.array(
                    [np.trace(P@j).real, np.trace(P@js).real]) / (N*N)
    return coef, (vmax, cmin)

def checks():
    assert np.linalg.norm(PT@PT.conj()+I4) < 1e-12
    for dx,dy in ((0,0),(.16,.23),(-.16,-.23),(.32,-.11)):
      for kx,ky in ((.5,.2),(1.2,-.7),(-1.1,.6)):
        for chiral in (False,True):
          H,v,j,js,inter,local = bloch(kx,ky,dx,dy,chiral)
          def anti(A): return PT @ A.conj() @ PT.conj().T
          assert np.linalg.norm(anti(H)-H) < 1e-12
          assert np.linalg.norm(anti(v)-v) < 1e-12
          assert np.linalg.norm(anti(j)+j) < 1e-12
          assert np.linalg.norm(anti(js)-js) < 1e-12
          hf,vf,jf,jsf,*_ = bloch(-kx,-ky,-dx,-dy,chiral)
          def turn(A): return C2 @ A @ C2.conj().T
          assert np.linalg.norm(turn(H)-hf) < 1e-11
          assert np.linalg.norm(turn(v)+vf) < 1e-11
          assert np.linalg.norm(turn(j)-jf) < 1e-11
          assert np.linalg.norm(turn(js)-jsf) < 1e-11
          # Bloch embedding: k_x -> k_x+2pi is a unitary sewing, not literal periodicity.
          W = K(np.diag([np.exp(1j*np.pi*dx), np.exp(-1j*np.pi*dx)]), I)
          hp, vp, *_ = bloch(kx+2*np.pi,ky,dx,dy,chiral)
          assert np.linalg.norm(W@H@W.conj().T-hp) < 1e-11
          assert np.linalg.norm(W@v@W.conj().T-vp) < 1e-11
          S = K(sz, sx)/2
          Tinter = 1j*(inter@S-S@inter)
          Tlocal = 1j*(local@S-S@local)
          assert np.linalg.norm(Tinter+Tlocal-1j*(H@S-S@H)) < 1e-12
    print('PASS: PT, C2z partner, embedded k sewing, exact torque continuity')
    for chiral in (False,True):
      for dx in (0,.12,.24,.33):
        plus,_ = fermi_surface(dx,.18,chiral=chiral,N=23)
        minus,_ = fermi_surface(-dx,-.18,chiral=chiral,N=23)
        assert np.linalg.norm(plus+minus)<3e-11
    print('PASS: C2z switching response parity and PT total cancellation')

    dx,dy=.24,.18
    print('metallic mu=1.3, T=.13, independent odd-dispersion hopping:')
    for amp in (0.,.4,.85,1.5):
       c41, _ = fermi_surface(dx,dy,chiral=True,tunnel_amp=amp,N=41)
       c61, _ = fermi_surface(dx,dy,chiral=True,tunnel_amp=amp,N=61)
       print('  t_amp=%4.2f stag(N41)=%+.9f stag(N61)=%+.9f' %
             (amp,c41[1],c61[1]))
       assert abs(c61[0])<1e-12 and abs(c41[1]-c61[1])<5e-6
    coupled, _=fermi_surface(dx,dy,chiral=True,tunnel_amp=0.,
                            gamma_tied_to_t=True,N=23)
    assert abs(coupled[1])<1e-12
    print('PASS: nonzero at t=0 in the independent gamma channel; '
          'forcing gamma proportional to t^2 removes the t=0 signal by assumption')

    c0,(vmax,cmin)=fermi_surface(dx,dy,mu=.346,T=.01,chiral=True,N=41)
    assert vmax<.346<cmin and abs(c0[1])<1e-7
    print('PASS: true gap, chemical potential in gap, no zero-T Drude signal')
    print('  valence_max %.7f, conduction_min %.7f, C_stag(T=.01) %.4g' %
          (vmax,cmin,c0[1]))

if __name__ == '__main__':
    checks()
