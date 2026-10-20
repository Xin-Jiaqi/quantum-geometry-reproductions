"""Exact 4x4 Schur/Feshbach resolvent with REMOTE-sector contact vertices.

This is an algebraically exact rewriting of the 8-band STATIC response,
not a controlled low-energy 4-band-only approximation or faster algorithm.
The decoupled middle-four subspace defines the P/Q basis at each fixed k.
Its k dependence causes NO error because v=dH/dk is transformed as an
operator *before* applying the resolvent derivative.

All tests are for synthetic 8-band, gapped T=0, equilibrium static NMEE;
no driven-DC claim or ferroelectric energetics.
"""
from dataclasses import replace
import numpy as np
from nmee_lab.clean_core.core import MU_B
from .sliding_tight_binding import SlidingBilayer, SlidingSpec


def reference_basis(model: SlidingBilayer, k):
    """8x8 unitary with decoupled middle 4 followed by 4 outer bands."""
    ref=SlidingBilayer(replace(model.spec,tunneling_eV=0.))
    _,eigvec=np.linalg.eigh(ref.matrices(np.asarray(k)[None,:])[0])
    return np.column_stack((eigvec[:,2:6],eigvec[:,:2],eigvec[:,6:]))


def feshbach(z, H, frame):
    """G=(z-H)^(-1)=U g V + R,  g=[z-A-B†rB]^-1."""
    h=frame.conj().T@H@frame
    A=h[:4,:4];B=h[4:,:4];D=h[4:,4:]
    r=np.linalg.inv(z*np.eye(4)-D)
    sigma=B.conj().T@r@B
    g=np.linalg.inv(z*np.eye(4)-A-sigma)
    U=np.vstack([np.eye(4),r@B])
    V=np.hstack([np.eye(4),B.conj().T@r])
    R=np.zeros((8,8),complex);R[4:,4:]=r
    G=frame@(U@g@V+R)@frame.conj().T
    return dict(frame=frame,g=g,U=U,V=V,R=R,G=G,sigma=sigma)


def inserted_source(components,O):
    """G O G: dressed 4D vertex and all three contact pieces.

    O_eff(z)=V(z) O U(z), corresponding to d/dsource[A+Sigma(z)].
    GOG=U*g*O_eff*g*V
        +U*g*V*O*R+R*O*U*g*V+R*O*R.
    All terms below are embedded back into the declared orbital basis.
    """
    F=components['frame'];g=components['g']
    U=components['U'];V=components['V'];R=components['R']
    o=F.conj().T@O@F
    oe=V@o@U
    pieces=dict(
        dressed=U@g@oe@g@V,
        left_contact=U@g@V@o@R,
        right_contact=R@o@U@g@V,
        remote_contact=R@o@R)
    return {key:F@M@F.conj().T for key,M in pieces.items()},oe


def _occupied_doublet_contours(e,n):
    """Numerical validation only: exact 8-band energies set contour locations.

    Separately integrate occupied clusters 0/1 and 2/3. The Feshbach
    *identity* requires no eigensystem and holds for any contour off poles.
    """
    if n<12 or e[2]<=e[1] or e[4]<=e[3]:
        raise ValueError("isolated occupied doublets and contour n>=12 required")
    circles=[]
    for low,high,left,right in ((0,1,None,2),(2,3,1,4)):
        z0=(e[low]+e[high])/2
        width=(e[high]-e[low])/2
        gap=min((e[low]-e[left]) if left is not None else np.inf,
                e[right]-e[high])
        radius=width+.3*gap
        ph=np.exp(2j*np.pi*np.arange(n)/n)
        circles.append((z0+radius*ph,radius*ph/n))
    return circles


def eigen_static_chi(H,spin,velocity):
    e,U=np.linalg.eigh(H)
    s=U.conj().T@spin@U;v=U.conj().T@velocity@U
    return float((-2*np.imag(s[:4,4:]*v[4:,:4].T
                   /(e[:4,None]-e[None,4:])**2)).sum())


def feshbach_static_at_k(spec:SlidingSpec,k,*,step_T=.25,contour_n=64):
    """Return exact static Gamma=partial_B chi and dressed-only diagnostic.

    Dressed-only omits the mixed U*g/R and Q-contact insertions, hence
    is NOT a physical reduced NMEE or a gauge-invariant observable.
    """
    if step_T<=0:raise ValueError("positive magnetic step needed")
    obj=SlidingBilayer(spec)
    k=np.asarray(k,float)
    if k.shape!=(2,):raise ValueError("k shape must be (2,)")
    frame=reference_basis(obj,k)
    H=obj.matrices(k[None,:])[0]
    vel=obj.matrices(k[None,:],axis=0)[0]
    sx=obj.spins['total_x'];sz=obj.spins['total_z']
    rows=[]
    for sign in (+1,-1):
        Hb=H-sign*obj.mono.g_factor*MU_B*step_T*sz
        e,U=np.linalg.eigh(Hb)
        if not e[3]<spec.mu_eV<e[4]:
            raise ValueError("chemical potential leaves the gap")
        Pexact=U[:,:4]@U[:,:4].conj().T
        P=np.zeros((8,8),complex)
        derivatives={x:np.zeros((8,8),complex) for x in
                     ('sx','velocity','sx_dressed','velocity_dressed')}
        maxres=0.
        for zs,ws in _occupied_doublet_contours(e,contour_n):
            for z,w in zip(zs,ws):
                source=feshbach(z,Hb,frame)
                Gdirect=np.linalg.inv(z*np.eye(8)-Hb)
                maxres=max(maxres,float(np.max(abs(Gdirect-source['G']))))
                P+=w*source['G']
                for O,key in ((sx,'sx'),(vel,'velocity')):
                    terms,eff=inserted_source(source,O)
                    total=sum(terms.values())
                    maxres=max(maxres,float(np.max(abs(total-Gdirect@O@Gdirect))))
                    derivatives[key]+=w*total
                    derivatives[key+'_dressed']+=w*terms['dressed']
        def chi(A,B):
            return float(np.real(1j*np.trace(Pexact@(A@B-B@A))))
        rows.append({'field_sign':sign,
            'chi_eigen':eigen_static_chi(Hb,sx,vel),
            'chi_exact_schur':chi(derivatives['sx'],derivatives['velocity']),
            'chi_dressed_only':chi(derivatives['sx_dressed'],
                                  derivatives['velocity_dressed']),
            'projector_error':float(np.max(abs(P-Pexact))),
            'resolvent_and_source_error':maxres,
            'remote_source_norm':float(np.linalg.norm(
                  derivatives['sx']-derivatives['sx_dressed'])),
            'min_mu_gap_eV':float(min(e[4]-spec.mu_eV,spec.mu_eV-e[3]))})
    numerical=lambda name:(rows[0][name]-rows[1][name])/(2*step_T)
    full=numerical('chi_exact_schur')
    dressed=numerical('chi_dressed_only')
    return {'k_inv_nm':k.tolist(),'u_fractional':list(spec.shift),
        't_perp_eV':spec.tunneling_eV,'Bstep_T':step_T,
        'contour_n_per_doublet':contour_n,
        'gamma_eigen':numerical('chi_eigen'),
        'gamma_feshbach_exact':full,
        'gamma_dressed_only_diagnostic':dressed,
        'gamma_missing_contacts':full-dressed,
        'plus_minus':rows,
        'protocol':'8band T=0 insulating equilibrium static, not driven DC'}

def feshbach_on_shell_vertices(spec:SlidingSpec,k):
    """Exact projected metric Z and energy-pair-dependent S/v/Zeeman vertices.

    For an exact eigenstate with P component p, q=r(E)B p and
    p†[I+B†r(E)^2B]p=1. For two energies Em,En the dressed vertex is
    O_PP+O_PQ*r(En)*B+B†*r(Em)*O_QP+B†*r(Em)*O_QQ*r(En)*B.
    The complete outer-band Kubo processes are NOT thereby eliminated.
    """
    obj=SlidingBilayer(spec)
    k=np.asarray(k,float)
    frame=reference_basis(obj,k)
    H=obj.matrices(k[None,:])[0]
    vx=obj.matrices(k[None,:],axis=0)[0]
    e,W=np.linalg.eigh(H)
    h=frame.conj().T@H@frame
    B=h[4:,:4];D=h[4:,4:]
    projected=frame.conj().T@W
    p=projected[:4,:]
    metrics=[];recon=[];resolvents=[]
    for n in range(2,6):
        r=np.linalg.inv(e[n]*np.eye(4)-D)
        Z=np.eye(4)+B.conj().T@r@r@B
        metrics.append(float(np.real(p[:,n].conj()@Z@p[:,n])))
        recon.append(float(np.linalg.norm(projected[4:,n]-r@B@p[:,n])))
        resolvents.append(r)
    result={}
    for name,O in [('spin_x',obj.spins['total_x']),
                   ('Zeeman_spin_z',obj.spins['total_z']),
                   ('velocity_x',vx)]:
        o=frame.conj().T@O@frame
        exact=W[:,2:6].conj().T@O@W[:,2:6]
        dressed=np.zeros((4,4),complex)
        bare=np.zeros((4,4),complex)
        for m in range(4):
            rm=resolvents[m]
            for n in range(4):
                rn=resolvents[n]
                oe=(o[:4,:4]+o[:4,4:]@rn@B+
                    B.conj().T@rm@o[4:,:4]+
                    B.conj().T@rm@o[4:,4:]@rn@B)
                dressed[m,n]=p[:,m+2].conj()@oe@p[:,n+2]
                bare[m,n]=p[:,m+2].conj()@o[:4,:4]@p[:,n+2]
        result[name]={'dressed_error':float(np.max(abs(dressed-exact))),
                      'bare_error':float(np.max(abs(bare-exact))),
                      'exact_norm':float(np.linalg.norm(exact))}
    return {'metric_normalization':metrics,
            'Q_amplitude_reconstruction_error':recon,'vertices':result,
            'scope':'four central exact eigenstates only; outer-band response remains'}
