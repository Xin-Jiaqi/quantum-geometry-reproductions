"""v61: SAME Fourier QWZ tight binding in central cross and semi-infinite strip leads.

This is a coherent, *charge* four-terminal edge measurement, not bulk orbital EB.
Lead principal layers exactly inherit the same H00 and nearest layer H01 as the
open QWZ model. We compare the measured Hall resistance against the *inverse*
quantized bulk conductivity in the clean, insulating, single-chiral-channel limit.
Finite eta and the cross contact geometry are explicit conditional assumptions.
"""
import numpy as np
from .periodic_kubo_v59 import sx,sy,sz,periodic_chern
from .qwz_finite_marker_v60 import qwz_open
from .matrix_leads_v60 import MatrixLead,matrix_caroli
from .lead_transport_v57 import four_terminal

TX=(sz-1j*sx)/2
TY=(sz-1j*sy)/2


def strip_principal_layer(width,mass,transverse,forward):
    """Strip onsite block for width sites with exact QWZ transverse hoppings."""
    if not isinstance(width,int) or width<2 or not np.isfinite(mass):
        raise ValueError('finite mass and width>=2')
    h=np.kron(np.eye(width),mass*sz).astype(complex)
    for p in range(width-1):
        i=2*p
        h[i:i+2,i+2:i+4]=transverse
        h[i+2:i+4,i:i+2]=transverse.conj().T
    v=np.kron(np.eye(width),forward)
    return h,v


def qwz_cross_leads(Lx,Ly,mass,eta=1e-5,contact_scale=1.):
    """Lead order L,Top,R,Bottom; exact continuation at contact_scale=1.

    Our Hamiltonian convention has H_{r,r+ex}=TX, H_{r,r+ey}=TY.
    Left surface->outer lead block is TX^dag; right is TX.
    Bottom is TY^dag; top is TY. Device->surface uses same block.
    """
    if not isinstance(Lx,int) or not isinstance(Ly,int) or min(Lx,Ly)<3:
        raise ValueError('Lx,Ly >=3')
    if not np.isfinite([mass,eta,contact_scale]).all() or eta<=0 or contact_scale<=0:
        raise ValueError('finite mass, eta>0 and scale>0')
    n=2*Lx*Ly
    contacts=[]
    setups=((Ly,TY,TX.conj().T,[(0,y) for y in range(Ly)]),
            (Lx,TX,TY,[(x,Ly-1) for x in range(Lx)]),
            (Ly,TY,TX,[(Lx-1,y) for y in range(Ly)]),
            (Lx,TX,TY.conj().T,[(x,0) for x in range(Lx)]))
    for width,transverse,forward,sites in setups:
        h00,h01=strip_principal_layer(width,mass,transverse,forward)
        c=np.zeros((n,2*width),complex)
        for j,(x,y) in enumerate(sites):
            row=2*(y*Lx+x)
            c[row:row+2,2*j:2*j+2]=contact_scale*forward
        contacts.append(MatrixLead(h00,h01,c,eta))
    return tuple(contacts)


def qwz_four_terminal(L=5,mass=-1,mu=0.,eta=1e-5,contact_scale=1.,Lx=None,Ly=None,min_reportable_current=1e-7):
    """At T=0, energy mu, compute coherent four-terminal dimensionless Rxy.

    Rxy/(h/q^2)=(V_top-V_bottom)/I_left for drive VL-VR=1.
    The periodic Chern computes sigma_xy/(q^2/h), distinct from Rxy.
    """
    Lx=L if Lx is None else Lx
    Ly=L if Ly is None else Ly
    if not np.isfinite([mu,min_reportable_current]).all() or min_reportable_current<=0:
        raise ValueError('finite chemical potential, positive current threshold')
    h,_=qwz_open(Lx,Ly,mass)
    leads=qwz_cross_leads(Lx,Ly,mass,eta,contact_scale)
    T,_=matrix_caroli(mu,h,np.zeros_like(h),leads)
    obs=four_terminal(T,np.zeros_like(T))
    il=obs['longitudinal_conductance_q2_over_h']
    reportable=abs(il)>=min_reportable_current
    rxy=float(obs['hall_ratio']/il) if reportable else None
    chern=periodic_chern(mass,31)
    return {'Lx':Lx,'Ly':Ly,'mass':mass,'mu':mu,'eta':eta,
            'contact_scale':contact_scale,
            'transmissions':T,'hall_ratio':obs['hall_ratio'],
            'drive_current_q2_over_h':il,
            'hall_resistance_h_over_q2':rxy,
            'hall_resistance_reportable':bool(reportable),
            'min_reportable_current_q2_over_h':float(min_reportable_current),
            'bulk_sigma_xy_q2_over_h':chern['chern_wilson'],
            'bulk_gap_min':chern['min_direct_gap'],
            'current_error':obs['max_current_conservation_residual'],
            'floating_error':obs['max_floating_probe_current'],
            'status':'CONDITIONAL_ZERO_T_COHERENT_4_TERMINAL_QWZ_EDGE_NOT_BULK_ORBITAL_EB'}
