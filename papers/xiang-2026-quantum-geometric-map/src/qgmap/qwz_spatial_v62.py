"""v62 source-independent spatial NEGF observable for the *same* QWZ lattice.

Finite 0-T linear-bias response of a coherent QWZ device. The four semi-infinite
QWZ strips have a controllable physical transverse width. The width changes
contact geometry: only full-edge leads produce a translational continuation
across that entire cut. No bulk dc transport limit or orbital EB is claimed.
"""
import numpy as np
from .qwz_finite_marker_v60 import qwz_open
from .qwz_device_v61 import TX,TY,strip_principal_layer
from .matrix_leads_v60 import MatrixLead,selfenergy_matrix
from .lead_transport_v57 import four_terminal


def qwz_width_leads(Lx,Ly,mass,width,eta=1e-4,offsets=(0,0,0,0)):
    """Lead order left, top, right, bottom; width is strip-site count.

    Physical QWZ hoppings remain unchanged; half-infinite strips couple only
    to selected *contiguous* sites on the device boundary, no wideband lead.
    offsets shifts contact within allowed boundary span (0 = centered).
    """
    if not isinstance(Lx,int) or not isinstance(Ly,int) or min(Lx,Ly)<3:
        raise ValueError('Lx,Ly must be integers >=3')
    if not isinstance(width,int) or not 2<=width<=min(Lx,Ly):
        raise ValueError('width must fit both sides, >=2')
    if len(offsets)!=4 or any(not isinstance(x,int) for x in offsets):
        raise ValueError('four integer side offsets')
    if not np.isfinite([mass,eta]).all() or eta<=0:
        raise ValueError('finite mass, eta>0')
    starts=[]
    for length,off in zip((Ly,Lx,Ly,Lx),offsets):
        start=(length-width)//2+off
        if start<0 or start+width>length: raise ValueError('contact outside boundary')
        starts.append(start)
    a,b,c,d=starts
    specs=((TY,TX.conj().T,[(0,y) for y in range(a,a+width)]),
           (TX,TY,[(x,Ly-1) for x in range(b,b+width)]),
           (TY,TX,[(Lx-1,y) for y in range(c,c+width)]),
           (TX,TY.conj().T,[(x,0) for x in range(d,d+width)]))
    n=2*Lx*Ly
    leads=[]
    for transverse,forward,sites in specs:
        h00,h01=strip_principal_layer(width,mass,transverse,forward)
        C=np.zeros((n,2*width),complex)
        for j,(x,y) in enumerate(sites):
            C[2*(y*Lx+x):2*(y*Lx+x)+2,2*j:2*j+2]=forward
        leads.append(MatrixLead(h00,h01,C,eta))
    return tuple(leads)


def spatial_charge_response(L=5,mass=-1.,width=3,mu=0.,eta=1e-4,offsets=(0,0,0,0),orbital_phases=None):
    """Independent local spectral and bond-current diagnostics, source at Fermi E.

    d_rho/dmu = sum_a V_a G Gamma_a G^dag (omits common 1/(2pi)).
    Internal bond current i->j / (q²/h per unit voltage) is
    -2 Im Tr[H_ij (d_rho)_{ji}]. A sign-fixed velocity operator follows
    directly from commutator [H, n_i], not a fitted edge-profile factor.
    """
    if not np.isfinite(mu):raise ValueError('finite mu')
    H,xy=qwz_open(L,L,mass)
    leads=qwz_width_leads(L,L,mass,width,eta,offsets)
    if orbital_phases is not None:
        phase=np.asarray(orbital_phases,float)
        if phase.shape!=(len(H),) or not np.isfinite(phase).all():
            raise ValueError('one finite phase per device orbital')
        U=np.diag(np.exp(1j*phase))
        H=U@H@U.conj().T
        leads=tuple(MatrixLead(l.h00,l.h01,U@l.coupling,l.eta) for l in leads)
    parts=[selfenergy_matrix(mu,l) for l in leads]
    Sigma=sum((p[0] for p in parts),np.zeros_like(H));Gam=[p[2] for p in parts]
    Gr=np.linalg.inv(mu*np.eye(len(H))-H-Sigma)
    Ga=Gr.conj().T
    inject=[Gr@G@Ga for G in Gam]
    transmissions=np.zeros((4,4))
    for a in range(4):
        for b in range(4):
            if a!=b:transmissions[a,b]=np.trace(Gam[a]@inject[b]).real
    o=four_terminal(transmissions,np.zeros_like(transmissions))
    bias_density=sum((v*a for v,a in zip(o['voltages'],inject)),np.zeros_like(H))
    # Identity: A = i(G^r-G^a) = sum_a G Gamma_a G^a
    A=1j*(Gr-Ga)
    spectral_error=float(np.linalg.norm(A-sum(inject)))
    LDOS=np.real(np.diag(A)).reshape(L,L,2).sum(axis=2)/(2*np.pi)
    injected=np.real(np.diag(bias_density)).reshape(L,L,2).sum(axis=2)
    Jx=np.zeros((L,L-1));Jy=np.zeros((L-1,L))
    def jbond(i,j):
        sl=slice(2*i,2*i+2);sr=slice(2*j,2*j+2)
        return float(-2*np.imag(np.trace(H[sl,sr]@bias_density[sr,sl])))
    for y in range(L):
        for x in range(L):
            i=y*L+x
            if x<L-1:Jx[y,x]=jbond(i,i+1)
            if y<L-1:Jy[y,x]=jbond(i,i+L)
    divergence=np.zeros((L,L))
    for y in range(L):
        for x in range(L):
            divergence[y,x]=(Jx[y,x] if x<L-1 else 0)-(Jx[y,x-1] if x else 0)+(Jy[y,x] if y<L-1 else 0)-(Jy[y-1,x] if y else 0)
    return {'L':L,'mass':mass,'width':width,'eta':eta,'mu':mu,
            'transmissions':transmissions,'voltages':o['voltages'],
            'drive_current_q2_over_h':o['longitudinal_conductance_q2_over_h'],
            'Rxy_h_over_q2':float(o['hall_ratio']/o['longitudinal_conductance_q2_over_h']) if abs(o['longitudinal_conductance_q2_over_h'])>1e-7 else None,
            'Jx':Jx,'Jy':Jy,'bond_divergence':divergence,'spectral_ldos':LDOS,
            'injection_site_weight':injected,'spectral_identity_error':spectral_error,
            'x_cut_currents':Jx.sum(axis=0),'y_cut_currents':Jy.sum(axis=1),
            'floating_current_error':o['max_floating_probe_current'],
            'total_lead_current_error':o['max_current_conservation_residual'],
            'status':'MODEL_NE_G_F_SPATIAL_CHARGE_ONLY_NO_BULK_ORBITAL_EB'}
