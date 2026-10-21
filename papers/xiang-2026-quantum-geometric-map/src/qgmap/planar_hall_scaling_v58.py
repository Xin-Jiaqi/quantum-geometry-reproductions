"""v58: Fixed-spacing nearest-neighbour planar Hall bars with explicit contact coverage.

Generalize v57's diagonal ideal wide-band lead self energy from one site per
reservoir to multiple independent local channels per reservoir, without
renormalizing the per-site Gamma. Analytic B derivative includes full Green
resolvent and self-consistent floating probe voltages.

A finite coherent device is NOT a bulk orbital Hall conductivity. Even a
fixed-density planar family may fail to converge in the coherent limit.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
from numpy.polynomial.legendre import leggauss
from .embedding_magnetic_v55 import peierls_field
from .lead_transport_v57 import four_terminal, laplacian, _weights, LeadSpec

@dataclass(frozen=True)
class PlanarSpec:
    nx:int=5
    ny:int=5
    lattice_spacing:float=1.
    hopping:float=1.
    onsite:float=0.
    coverage:float=1.
    gamma_per_site:float=.6
    gamma_mode:str='per_site'
    def __post_init__(self):
        if not isinstance(self.nx,int) or not isinstance(self.ny,int) or self.nx<3 or self.ny<3:
            raise ValueError('nx,ny integers >=3 for disjoint four-edge contacts')
        vals=(self.lattice_spacing,self.hopping,self.onsite,self.coverage,self.gamma_per_site)
        if not np.isfinite(vals).all() or self.lattice_spacing<=0 or self.hopping<=0 or not (0<self.coverage<=1) or self.gamma_per_site<=0 or self.gamma_mode not in ('per_site','fixed_reservoir_sum'):
            raise ValueError('finite lattice; spacing/hopping/gamma>0, 0<coverage<=1')

def planar_bar(spec:PlanarSpec=PlanarSpec()):
    """Same lattice sites per area a^-2, same NN hopping as nx,ny increase.

    Non-overlapping four reservoirs on edge interiors L,T,R,B. Each contacted
    edge site adds a physically separate local wide-band channel with gamma.
    """
    nx,ny=spec.nx,spec.ny
    pos=np.array([[x*spec.lattice_spacing,y*spec.lattice_spacing]
                  for y in range(ny) for x in range(nx)],float)
    H=np.eye(nx*ny,dtype=complex)*spec.onsite
    for y in range(ny):
        for x in range(nx):
            i=y*nx+x
            for xx,yy in ((x+1,y),(x,y+1)):
                if xx<nx and yy<ny:
                    j=yy*nx+xx
                    H[i,j]=H[j,i]=-spec.hopping
    def centered_coords(dim):
        available=np.arange(1,dim-1,dtype=int)
        # Select a centered mirror-symmetric contact set: same parity as edge.
        target=spec.coverage*len(available)
        counts=[c for c in range(1,len(available)+1) if (len(available)-c)%2==0]
        count=min(counts,key=lambda c:(abs(c-target),c))
        offset=(len(available)-count)//2
        return available[offset:offset+count]
    ys=centered_coords(ny);xs=centered_coords(nx)
    leads=(tuple(int(y*nx) for y in ys),
           tuple(int((ny-1)*nx+x) for x in xs),
           tuple(int(y*nx+nx-1) for y in ys),
           tuple(int(x) for x in xs))
    return pos,H,leads

def gamma_profiles(sites,number_sites,gamma_per_site,mode='per_site'):
    """Four separate positive local Γ_alpha diagonal profiles (4,N).

    Disallow overlapping leads to retain unambiguous reservoir groups.
    """
    if not isinstance(number_sites,int) or number_sites<4 or len(sites)!=4 or not np.isfinite(gamma_per_site) or gamma_per_site<=0 or mode not in ('per_site','fixed_reservoir_sum'):
        raise ValueError('4 terminal profiles, positive coupling and >=4 sites')
    result=np.zeros((4,number_sites),float)
    used=set()
    for alpha,lead in enumerate(sites):
        if len(lead)==0:raise ValueError('each lead requires >=1 distinct site')
        for s in lead:
            if not isinstance(s,(int,np.integer)) or not 0<=s<number_sites or s in used:
                raise ValueError('leads require disjoint valid sites')
            used.add(s)
            result[alpha,s]=gamma_per_site if mode=='per_site' else gamma_per_site/len(lead)
    return result

def caroli_multisite(E,H,M,gamma):
    """T_ab=sum_(i in a,j in b) Γ_ai Γ_bj |G_ij|², exact for diagonal Γ.

    The orientation follows Tr Γ_a G Γ_b G†. Magnetic derivative uses
    dG=G M G valid only if Γ is B-independent.
    """
    gamma=np.asarray(gamma,float)
    H=np.asarray(H,complex)
    n=H.shape[0]
    if H.shape!=(n,n) or gamma.shape!=(4,n) or M.shape!=H.shape:
        raise ValueError('four independent diagonal self energy profiles')
    if not np.isfinite(gamma).all() or np.any(gamma<0) or np.any(np.sum(gamma,axis=1)<=0):
        raise ValueError('positive lead Γ and valid entries')
    G=np.linalg.inv(E*np.eye(n,dtype=complex)-H+.5j*np.diag(gamma.sum(axis=0)))
    dG=G@M@G
    td=np.abs(G)**2
    dt=2*np.real(G.conj()*dG)
    T=gamma@td@gamma.T
    derivative=gamma@dt@gamma.T
    np.fill_diagonal(T,0)
    np.fill_diagonal(derivative,0)
    return T,derivative

def planar_response(planar=PlanarSpec(),lead=LeadSpec(mu=.65,temperature=.12,field=.12,
                                                    charge_over_hbar=-.87,energy_points=80),*,gauge='symmetric'):
    positions,H0,contact_sites=planar_bar(planar)
    coupling=gamma_profiles(contact_sites,len(H0),planar.gamma_per_site,planar.gamma_mode)
    p=peierls_field(H0,positions,lead.field,lead.charge_over_hbar,gauge=gauge)
    energies,weights=_weights(lead)
    T=np.zeros((4,4),float);dT=np.zeros((4,4),float)
    for E,w in zip(energies,weights):
        t,d=caroli_multisite(E,p.hamiltonian,p.magnetic_vertex,coupling)
        T+=w*t;dT+=w*d
    observed=four_terminal(T,dT)
    Ix=observed['longitudinal_conductance_q2_over_h']
    Vh=observed['hall_ratio']
    if abs(Ix)<1e-12:raise ValueError('device source current nearly zero; Hall resistance unstable')
    observed['hall_resistance_h_over_q2']=float(Vh/Ix)
    observed['dHall_resistance_dB']=float((observed['dHall_dB']*Ix-Vh*observed['dLongitudinal_dB'])/Ix**2)
    observed['longitudinal_resistance_h_over_q2']=float(1/Ix)
    observed['longitudinal_sheet_proxy_h_over_q2']=float((planar.ny-1)/(planar.nx-1)/Ix)
    aspect=(planar.nx-1)/(planar.ny-1)
    return {'T':T,'dT_dB':dT,'observables':observed,'positions':positions,
            'contact_sites':contact_sites,'contact_site_counts':tuple(map(len,contact_sites)),
            'gamma_total_by_reservoir':tuple(map(float,coupling.sum(axis=1))),
            'gamma_mode':planar.gamma_mode,
            'sample_sites':len(H0),'lattice_density':1/planar.lattice_spacing**2,
            'device_area':(planar.nx-1)*(planar.ny-1)*planar.lattice_spacing**2,
            'aspect_ratio':aspect,'plaquette_flux':lead.charge_over_hbar*lead.field*planar.lattice_spacing**2,
            'status':'VERIFIED_COHERENT_PLANAR_CONTACT_MODEL_NOT_BULK'}

def full_B_derivative(planar=PlanarSpec(),lead=LeadSpec(mu=.65,temperature=.12,field=.12,
                                                         charge_over_hbar=-.87,energy_points=80),*,gauge='symmetric',step=2e-4):
    if not np.isfinite(step) or step<=0:raise ValueError('positive finite field step')
    f=lambda B: planar_response(planar,replace(lead,field=B),gauge=gauge)
    p2,p1,n1,n2=[f(lead.field+d*step) for d in (2,1,-1,-2)]
    deriv=lambda fun:(-fun(p2)+8*fun(p1)-8*fun(n1)+fun(n2))/(12*step)
    return {'dT_dB':deriv(lambda d:d['T']),
            'dHall_dB':float(deriv(lambda d:d['observables']['hall_ratio'])),
            'dLongitudinal_dB':float(deriv(lambda d:d['observables']['longitudinal_conductance_q2_over_h'])),
            'dHall_resistance_dB':float(deriv(lambda d:d['observables']['hall_resistance_h_over_q2']))}

def compare_planar_family(configs, lead, *,gauge='symmetric'):
    """No automatic thermodynamic convergence claim; surface raw device data."""
    rows=[]
    for plan in configs:
        x=planar_response(plan,lead,gauge=gauge)
        rows.append({'nx':plan.nx,'ny':plan.ny,'sites':x['sample_sites'],
                     'coverage':plan.coverage,'contact_counts':x['contact_site_counts'],
                     'gamma_mode':x['gamma_mode'],
                     'gamma_total_by_reservoir':x['gamma_total_by_reservoir'],
                     'aspect_ratio':x['aspect_ratio'],'density':x['lattice_density'],
                     'plaquette_flux':x['plaquette_flux'],
                     'hall_ratio':x['observables']['hall_ratio'],
                     'dHall_dB':x['observables']['dHall_dB'],
                     'longitudinal':x['observables']['longitudinal_conductance_q2_over_h'],
                     'hall_resistance_h_over_q2':x['observables']['hall_resistance_h_over_q2'],
                     'dHall_resistance_dB':x['observables']['dHall_resistance_dB'],
                     'floating_current_error':x['observables']['max_floating_probe_current']})
    return {'status':'CONDITIONAL_FIXED_DENSITY_COHERENT_FINITE_SAMPLE_SCAN_NOT_BULK',
            'rows':rows,'thermodynamic_claim':False,
            'caveats':['Finite coherent sample with ideal local wideband electrodes',
                       'Fixed lattice density and fixed per-site gamma do not make contacts irrelevant',
                       'Hall voltage is a contact/device observable, not a bulk conductivity',
                       'No mean free path control or disorder/self-energy thermodynamic limit']}


def adaptive_planar_response(planar=PlanarSpec(), lead=LeadSpec(mu=.65,temperature=.12,field=.12,
                                                           charge_over_hbar=-.87,energy_points=80),
                             *,gauge='symmetric',epsrel=3e-8,epsabs=1e-10):
    """Independently adaptive Fermi-window integral for sharp mesoscopic resonances.

    Use eigenenergy splitting points and scipy quad_vec error estimate; never
    assume nE=80 Gaussian nodes are converged for a coherent finite device.
    """
    from scipy.integrate import quad_vec
    if not np.isfinite([epsrel,epsabs]).all() or epsrel<=0 or epsabs<=0:
        raise ValueError('positive adaptive tolerances')
    positions,H0,contacts=planar_bar(planar)
    gamma=gamma_profiles(contacts,len(H0),planar.gamma_per_site,planar.gamma_mode)
    p=peierls_field(H0,positions,lead.field,lead.charge_over_hbar,gauge=gauge)
    a=lead.mu-lead.window_in_kT*lead.temperature
    b=lead.mu+lead.window_in_kT*lead.temperature
    eig=np.linalg.eigvalsh(p.hamiltonian)
    points=sorted(set(float(v) for v in eig if a<v<b))
    evals=[0]
    def function(E):
        evals[0]+=1
        f=1/(4*lead.temperature*np.cosh((E-lead.mu)/(2*lead.temperature))**2)
        t,dt=caroli_multisite(E,p.hamiltonian,p.magnetic_vertex,gamma)
        return np.r_[t.ravel(),dt.ravel()]*f
    vec,error,info=quad_vec(function,a,b,epsabs=epsabs,epsrel=epsrel,points=points,limit=500,full_output=True)
    if not info.success:
        raise RuntimeError('adaptive Fermi integration did not converge; increase numerical limit or change physical parameters')
    T=vec[:16].reshape(4,4);dT=vec[16:].reshape(4,4)
    obs=four_terminal(T,dT)
    current=obs['longitudinal_conductance_q2_over_h']
    if abs(current)<1e-12:raise ValueError('nearly zero source current')
    obs['hall_resistance_h_over_q2']=float(obs['hall_ratio']/current)
    obs['dHall_resistance_dB']=float((obs['dHall_dB']*current-obs['hall_ratio']*obs['dLongitudinal_dB'])/current**2)
    return {'T':T,'dT_dB':dT,'observables':obs,
            'contact_sites':contacts,'contact_site_counts':tuple(map(len,contacts)),
            'gamma_total_by_reservoir':tuple(map(float,gamma.sum(axis=1))),
            'gamma_mode':planar.gamma_mode,
            'sample_sites':len(H0),'lattice_density':1/planar.lattice_spacing**2,
            'device_area':(planar.nx-1)*(planar.ny-1)*planar.lattice_spacing**2,
            'aspect_ratio':(planar.nx-1)/(planar.ny-1),
            'plaquette_flux':lead.charge_over_hbar*lead.field*planar.lattice_spacing**2,
            'quadrature_error_estimate':float(error),'integrand_evaluations':evals[0],
            'quadrature_success':bool(info.success),
            'source':'adaptive_fermi_window_split_at_Bloch_eigenenergies',
            'status':'VERIFIED_COHERENT_PLANAR_CONTACT_MODEL_NOT_BULK'}
