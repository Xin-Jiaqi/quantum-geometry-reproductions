"""Exact weak-hopping static NMEE coefficients from gapped projector algebra.

No hopping-amplitude, spin-source, k, or Zeeman finite differences.
Uses P^2=P and [H,P]=0 in the decoupled occupied/empty eigenbasis.
Only occupied-to-empty spectral denominators appear, so degeneracies
within either separated subspace do not require nondegenerate formulas.

Variables: x and y are independent hopping amplitudes, b is Bz in Tesla,
lambda couples the physical Sx operator, q is a local kx derivative.
Source H = H0 + x X + y Y - g mu_B b Sz + lambda Sx
           + q v0 + x*q d_k X + y*q d_k Y.
The coefficient of x*y*b in i Tr P[P_lambda,P_q] is the mixed
static NMEE derivative. For x=y denoting two formal copies of the SAME
hopping, divide by 2 to obtain the quadratic coefficient.
Gauge: the original physical Wannier orbital embedding; do not use
a k-dependent gauge transformation without its connection terms.
"""
from dataclasses import replace
import numpy as np
from nmee_lab.clean_core.core import MU_B
from .sliding_tight_binding import SlidingBilayer, SlidingSpec


def _subsets(mask):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            return
        sub = (sub-1) & mask


def projector_coefficient(H0, X, Y, v0, vx, vy, spinx, spinz, mu,
                          *, g_factor=10., part="all"):
    """d_x d_y d_B mixed spin-momentum Berry curvature (array of N)."""
    H0 = np.asarray(H0)
    n = len(H0)
    if H0.ndim != 3 or H0.shape[1] != 8 or H0.shape[2] != 8:
        raise ValueError("Expected Nx8x8 8-band model")
    e, U = np.linalg.eigh(H0)
    if not np.all((e[:, 3] < mu) & (mu < e[:, 4])):
        raise ValueError("Decoupled state must be insulating at fixed mu")
    UD = U.conj().swapaxes(1, 2)
    frame = lambda a: UD @ a @ U
    vertices = {
        1: frame(X), 2: frame(Y),
        4: frame(-g_factor*MU_B*spinz), 8: frame(spinx),
        16: frame(v0),
    }
    if part in ("all", "x"):
        vertices[17] = frame(vx)
    if part in ("all", "y"):
        vertices[18] = frame(vy)
    if part not in ("all", "none", "x", "y"):
        raise ValueError("part must be all/none/x/y")
    P = {0: np.broadcast_to(
        np.diag([1]*4+[0]*4).astype(complex), (n,8,8)).copy()}
    eigdiff = e[:, :, None] - e[:, None, :]
    cross = np.zeros((8,8),bool)
    cross[:4,4:] = True
    cross[4:,:4] = True
    for mask in range(1,32):
        prod = np.zeros_like(H0)
        for a in _subsets(mask):
            if a and a != mask:
                prod += P[a] @ P[mask ^ a]
        source = np.zeros_like(H0)
        for a, vertex in vertices.items():
            if mask & a == a:
                remaining = P[mask ^ a]
                source -= vertex @ remaining - remaining @ vertex
        pm = np.zeros_like(H0)
        pm[:,:4,:4] = -prod[:,:4,:4]
        pm[:,4:,4:] = prod[:,4:,4:]
        pm[:,cross] = source[:,cross]/eigdiff[:,cross]
        P[mask] = pm
    answer = np.zeros(n,dtype=complex)
    for a in _subsets(7):  # 7 = x|y|B
        left = 7 ^ a
        for b in _subsets(left):
            c = left ^ b
            pl, pk = P[b|8], P[c|16]
            answer += 1j*np.einsum("nij,nji->n",
                                    P[a], pl@pk-pk@pl, optimize=True)
    if np.max(np.abs(answer.imag)) > 1e-6*max(1., np.max(abs(answer.real))):
        raise ArithmeticError("Unexpected complex projector response")
    return np.real(answer)


def _polar_nodes(radius, nr, nt):
    if not radius > 0 or nr < 2 or nt < 4:
        raise ValueError("Bad valley disk quadrature")
    z, gw = np.polynomial.legendre.leggauss(nr)
    r = radius*(z+1)/2
    rw = radius*gw/2
    phi = 2*np.pi*(np.arange(nt)+.5)/nt
    q = (r[:,None,None] * np.stack(
        (np.cos(phi),np.sin(phi)),axis=-1)[None,:,:]).reshape(-1,2)
    weight = np.repeat(r*rw,nt)/(2*np.pi*nt)
    return q, weight


def analytic_valley_ABC(spec=SlidingSpec(), *, radius=.75, nr=26, nt=72):
    """Return A,B,C in Gamma=A s²+B c²+C sc+O(hopping^4).

    k disks and declared A/B orbital decomposition only; not full BZ.
    The spec shift is in fractional primitive coordinates; reference
    t_perp in spec does not enter derivatives at zero hopping.
    """
    no_hop=SlidingBilayer(replace(spec,tunneling_eV=0.))
    same=SlidingBilayer(replace(spec,tunneling_eV=1.,cross_sublattice=False))
    all_hop=SlidingBilayer(replace(spec,tunneling_eV=1.,cross_sublattice=True))
    reciprocal=no_hop.model().reciprocal
    centers=np.array([[1/3,2/3],[2/3,1/3]])@reciprocal
    q,w=_polar_nodes(radius,nr,nt)
    answer={name:{term:0. for term in ("v0","vx","vy","total")}
            for name in ("A_same2","B_cross2","C_mixed")}
    valleys=[]
    for valley,center in enumerate(centers):
        k=q+center
        h0=no_hop.matrices(k); v0=no_hop.matrices(k,axis=0)
        hs=same.matrices(k); vs=same.matrices(k,axis=0)
        ha=all_hop.matrices(k); va=all_hop.matrices(k,axis=0)
        Xs=hs-h0; Xc=ha-hs
        Wx=vs-v0; Wc=va-vs
        per={}
        for name,x,y,wx,wy,factor in (
            ("A_same2",Xs,Xs,Wx,Wx,.5),
            ("B_cross2",Xc,Xc,Wc,Wc,.5),
            ("C_mixed",Xs,Xc,Wx,Wc,1.)):
            vals={p:float(w@projector_coefficient(
                h0,x,y,v0,wx,wy,
                no_hop.spins["total_x"],no_hop.spins["total_z"],
                spec.mu_eV,g_factor=no_hop.mono.g_factor,part=p))*factor
                for p in ("none","x","y","all")}
            decomposition={
                "v0":vals["none"], "vx":vals["x"]-vals["none"],
                "vy":vals["y"]-vals["none"],"total":vals["all"]}
            per[name]=decomposition
            for part,num in decomposition.items():
                answer[name][part] += num
        valleys.append(per)
    return {
        "protocol":"0 K gapped equilibrium/static d_B Berry curvature; no driven DC",
        "scope":"two valley disks, not entire BZ",
        "radius_inv_nm":radius,"nr":nr,"nt":nt,"shift":list(spec.shift),
        "coefficients":answer,"valley_decomposition":valleys,
        "caveat":"v0/vx/vy split requires declared orbital embedding; total A/B/C do not depend on this bookkeeping at fixed physical Hamiltonian",
    }
