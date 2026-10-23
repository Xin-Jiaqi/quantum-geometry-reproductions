"""v76: source-independent two-band Dirac elastic collision and coherent EB.

STRICT DOMAIN: beta-tilted two-band *population* on true isoenergy shells;
UNTILTED beta=0 only for the same-k 2x2 interband GKSL coherence block.
This is NOT the original Xiang SM length-gauge two-color/selected tau pole:
its electric vertex below includes only the projected *interband* r_y,
not electric momentum drift. A collision superoperator can have the same
transport lifetime but different EB interband responses.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.linalg import expm

I = np.eye(2,dtype=complex)
SX = np.array([[0,1],[1,0]],complex)
SY = np.array([[0,-1j],[1j,0]],complex)
SZ = np.diag([1.,-1.]).astype(complex)

@dataclass(frozen=True)
class Spec:
    n:int=16
    energies:tuple=(.3,.6)
    beta:float=.125
    gamma:float=.07
    mu:float=.19
    temperature:float=.16
    omega_e:float=.13
    omega_b:float=.18
    eta_e:float=.04
    eta_b:float=.03
    pure_phi:float=0.
    def __post_init__(self):
        vals=(self.beta,self.gamma,self.mu,self.temperature,self.omega_e,
              self.omega_b,self.eta_e,self.eta_b,self.pure_phi)
        if not np.isfinite(vals).all() or abs(self.beta)>=1 or self.gamma<=0 or self.temperature<=0 or min(self.eta_e,self.eta_b)<=0 or self.pure_phi<0:
            raise ValueError('finite type-I tilt, rates and temperature; gamma, eta positive')
        if not isinstance(self.n,int) or self.n<12 or self.n%4:
            raise ValueError('n must be >=12 divisible by 4')
        if not isinstance(self.energies,tuple) or len(self.energies)<1 or any((not np.isfinite(e)) or e<=0 for e in self.energies):
            raise ValueError('strictly positive finite energy shells required')

def fermi(x,mu,T):return 1/(1+np.exp(np.clip((np.asarray(x)-mu)/T,-700,700)))

def shell(s:Spec,E,band,B=0.):
    """Exact tilted Dirac isoenergy contour E_band = band*E even under B sigma_x.

    H(k,B)=beta*kx I+kx sy-(ky+B) sx; lab ky=r sin(theta)-B.
    Spinor momentum q=(kx,ky+B) is a rigid translation. The impurity
    short-range scalar overlap depends only on q-hat separation.
    """
    if band not in (-1,1) or not np.isfinite(B):raise ValueError('band +/-1 and finite B required')
    theta=2*np.pi*np.arange(s.n)/s.n
    co,si=np.cos(theta),np.sin(theta)
    denom=1+band*s.beta*co
    r=E/denom
    kx=r*co;ky=r*si-B
    h=s.beta*kx[:,None,None]*I+kx[:,None,None]*SY-(ky+B)[:,None,None]*SX
    evals=np.linalg.eigvalsh(h)
    en=evals[:,int(band==1)]
    weight=denom**-2
    weight=weight/np.sum(weight)
    mat=s.gamma*weight[None,:]*(1+np.cos(theta[:,None]-theta[None,:]))
    np.fill_diagonal(mat,0.)
    L=np.diag(np.sum(mat,axis=1))-mat
    vx=s.beta+band*co;vy=band*si
    return dict(theta=theta,radius=r,kx=kx,ky=ky,weights=weight,transition=mat,
                L=L,band_energy=en,dos_angle=denom**-2,vx=vx,vy=vy,
                B=B,gamma=s.gamma,band=band)

def band_population(s:Spec,band,omega=0.,B=0.):
    if not np.isfinite(omega):raise ValueError('finite omega required')
    entries=[]
    total=0j
    for E in s.energies:
        d=shell(s,E,band,B)
        w=d['weights'];v=d['vy'];v0=v-np.dot(w,v)
        P=np.outer(np.ones(s.n),w)
        g=np.linalg.solve(d['L']-1j*omega*I_n(s.n)+P,v0)
        f=fermi(band*E,s.mu,s.temperature)
        fac=f*(1-f)/s.temperature*E*np.mean(d['dos_angle'])
        current=fac*np.dot(w,v*g)
        total+=current
        entries.append(dict(E=E,band=band,B=B,source_mean=np.dot(w,v0),
            carrier_weight=fac,current=current,distribution=g,
            number=np.dot(w,g),collision_charge=np.linalg.norm(w@d['L']),
            linear_residual=np.linalg.norm((d['L']-1j*omega*I_n(s.n))@g-v0)))
    return dict(total_current=total,entries=entries)

def I_n(n):return np.eye(n,dtype=complex)

def two_band_populations(s=Spec(),omega=0.,B=0.):
    ans={b:band_population(s,b,omega,B) for b in (-1,1)}
    return dict(bands=ans,current=sum(a['total_current'] for a in ans.values()),
        B=B,status='CONSERVING_BAND_AND_ENERGY_RESOLVED_BOLTZMANN_NOT_FULL_EB')

@dataclass(frozen=True)
class CoherentSpec:
    n:int=16
    energy:float=.4
    mass:float=.21
    mu:float=.08
    temperature:float=.14
    gamma:float=.07
    pure_phi:float=0.
    omega_e:float=.13
    omega_b:float=.18
    eta_e:float=.04
    eta_b:float=.03
    def __post_init__(self):
        vals=(self.energy,self.mass,self.mu,self.temperature,self.gamma,self.pure_phi,
            self.omega_e,self.omega_b,self.eta_e,self.eta_b)
        if not np.isfinite(vals).all() or self.energy<=0 or self.temperature<=0 or self.gamma<=0 or self.pure_phi<0 or min(self.eta_e,self.eta_b)<=0:
            raise ValueError('positive energy/rates and finite parameters required')
        if not isinstance(self.n,int) or self.n<12 or self.n%4:raise ValueError('n divisible by 4, >=12')

def coherent_data(s=CoherentSpec()):
    """Exact beta=0 equal-radius two-band spin texture and GKSL jump rates.

    For each k, 2x2 band density. Jump |theta_j,b><theta_i,b| with
    positive rates W_ij; momentum-offdiagonal density is not included.
    Optional positive local pure dephasing leaves population scattering L.
    """
    th=2*np.pi*np.arange(s.n)/s.n
    h=s.energy*(np.cos(th)[:,None,None]*SY-np.sin(th)[:,None,None]*SX)+s.mass*SZ
    e,u=np.linalg.eigh(h)
    ud=np.swapaxes(u.conj(),-1,-2)
    vx=ud@SY@u
    vy=ud@(-SX)@u
    spin=ud@SX@u
    ff=fermi(e,s.mu,s.temperature)
    rho0=np.zeros_like(h)
    rho0[:,0,0]=ff[:,0];rho0[:,1,1]=ff[:,1]
    weight=np.ones(s.n)/s.n
    # Scalar-impurity spinor overlap of the MASSIVE Dirac spin texture.
    # 2|<u_i,b|u_j,b>|^2 = 1+n_i.n_j, same for both bands.
    mm=s.mass*s.mass;rr=s.energy*s.energy
    overlap=1+(rr*np.cos(th[:,None]-th[None,:])+mm)/(rr+mm)
    W=s.gamma*weight[None,:]*overlap
    np.fill_diagonal(W,0)
    L=np.diag(W.sum(axis=1))-W
    rate=W.sum(axis=1)
    # Hermitian off-diagonal position matrix from v_y = i (eps_m-eps_n) r_mn.
    gap=e[..., :,None]-e[...,None,:]
    ry=np.zeros_like(vy)
    for m,n in ((0,1),(1,0)):
        ry[:,m,n]=vy[:,m,n]/(1j*gap[:,m,n])
    return dict(h=h,e=e,u=u,vx=vx,vy=vy,spin=spin,ry=ry,rho0=rho0,
                W=W,L=L,rate=rate,theta=th)

def apply_generator(s:CoherentSpec,rho,data=None):
    """-i[H,rho]+D_GKSL[rho], preserves total band-population trace."""
    d=coherent_data(s) if data is None else data
    r=np.asarray(rho,complex)
    if r.shape!=(s.n,2,2):raise ValueError('rho shape n,2,2')
    out=np.zeros_like(r)
    for b in (0,1):
        out[:,b,b]=-d['L']@r[:,b,b]
    rate=d['rate']+s.pure_phi
    gap=2*np.hypot(s.energy,s.mass)
    out[:,0,1]=(-rate+1j*gap)*r[:,0,1]
    out[:,1,0]=(-rate-1j*gap)*r[:,1,0]
    return out

def vertex(v,rho):return -1j*(v@rho-rho@v)

def _resolvent_spectral(s,source,z,data=None):
    d=coherent_data(s) if data is None else data
    if np.real(z)<=0 or not np.isfinite([z.real,z.imag]).all():raise ValueError('retarded z real >0')
    src=np.asarray(source,complex)
    sol=np.zeros_like(src)
    L=d['L'].astype(complex)
    # eigenbasis of circulant real positive population Markov Laplacian
    lam,u=np.linalg.eigh((L+L.T)/2)
    for b in (0,1):
        sol[:,b,b]=u@((u.T@src[:,b,b])/(z+lam))
    rate=d['rate']+s.pure_phi
    gap=2*np.hypot(s.energy,s.mass)
    sol[:,0,1]=src[:,0,1]/(z+rate-1j*gap)
    sol[:,1,0]=src[:,1,0]/(z+rate+1j*gap)
    return sol

def _resolvent_dense(s,source,z,data=None):
    d=coherent_data(s) if data is None else data
    n=s.n*4
    basis=np.eye(n,dtype=complex)
    A=np.column_stack([(z*basis[:,i]-apply_generator(s,basis[:,i].reshape(s.n,2,2),d).reshape(n)) for i in range(n)])
    return np.linalg.solve(A,source.reshape(n)).reshape(s.n,2,2)

def coherent_two_color(s=CoherentSpec(),method='spectral'):
    """Two-color COHERENCE-ONLY projected r_y electric and sigma_x spin sources.

    This is a fully specified frozen Lindblad protocol, but the electric
    vertex lacks the uniform length-gauge population *drift* D_y rho and
    its E interaction with rho_B. It CANNOT reproduce SM Eq20 tau residue.
    """
    if method not in ('spectral','dense'):raise ValueError('spectral or dense')
    d=coherent_data(s);res=_resolvent_spectral if method=='spectral' else _resolvent_dense
    ze=s.eta_e-1j*s.omega_e;zb=s.eta_b-1j*s.omega_b;zs=ze+zb
    rE=res(s,vertex(d['ry'],d['rho0']),ze,d)
    rB=res(s,vertex(d['spin'],d['rho0']),zb,d)
    BE=res(s,vertex(d['ry'],rB),zs,d)
    EB=res(s,vertex(d['spin'],rE),zs,d)
    val=lambda x:complex(np.mean(np.trace(d['vx']@x,axis1=-2,axis2=-1)))
    charge=lambda x:complex(np.sum(np.trace(x,axis1=-2,axis2=-1)))
    local=np.trace(d['vx']@(BE+EB),axis1=-2,axis2=-1)
    harmonic=complex(np.mean(np.exp(1j*d['theta'])*local))
    return dict(b_then_e=val(BE),e_then_b=val(EB),total=val(BE+EB),
        local_harmonic_1=harmonic,local_mixed_norm=float(np.linalg.norm(BE+EB)),
        first_e=val(rE),first_b=val(rB),trace_first_e=charge(rE),
        trace_first_b=charge(rB),trace_second=charge(BE+EB),
        two_color_order_residual=float(abs(val(BE+EB)-val(BE)-val(EB))),
        gamma=s.gamma,pure_phi=s.pure_phi,
        status='CONSERVING_TWO_BAND_COHERENT_FINITE_VERTEX_NOT_AUTHOR_LENGTH_GAUGE')

def coherent_static_spin(s=CoherentSpec(),method='spectral'):
    d=coherent_data(s);res=_resolvent_spectral if method=='spectral' else _resolvent_dense
    rB=res(s,vertex(d['spin'],d['rho0']),s.eta_b,d)
    zval=complex(np.mean(np.trace(d['spin']@rB,axis1=-2,axis2=-1)))
    return dict(spin_susceptibility=zval,pure_phi=s.pure_phi,
        diagonal_source_norm=float(np.linalg.norm(np.diagonal(vertex(d['spin'],d['rho0']),axis1=1,axis2=2))))

def decoherence_nonuniqueness(s=CoherentSpec()):
    a=coherent_two_color(s)
    b=coherent_two_color(replace(s,pure_phi=.14))
    x=coherent_static_spin(s)
    y=coherent_static_spin(replace(s,pure_phi=.14))
    da=coherent_data(s);db=coherent_data(replace(s,pure_phi=.14))
    return dict(population_L_error=float(np.linalg.norm(da['L']-db['L'])),
        eb_zero_phi=a['total'],eb_phi=b['total'],
        eb_difference=complex(b['total']-a['total']),
        local_harmonic_zero_phi=a['local_harmonic_1'],
        local_harmonic_phi=b['local_harmonic_1'],
        local_harmonic_difference=b['local_harmonic_1']-a['local_harmonic_1'],
        spin_zero_phi=x['spin_susceptibility'],spin_phi=y['spin_susceptibility'],
        spin_difference=complex(y['spin_susceptibility']-x['spin_susceptibility']),
        message='IDENTICAL_POPULATION_COLLISION_TAU_DISTINCT_INTERBAND_COHERENCE_AND_EB')

def time_evolution(s=CoherentSpec(),initial=None,t=20.):
    """Direct independent finite-time GKSL: population exp(-L t), coherence exp(-Gamma t)."""
    if not np.isfinite(t) or t<0:raise ValueError('nonnegative finite t')
    d=coherent_data(s)
    if initial is None:
        initial=np.zeros((s.n,2,2),complex)
        for k in range(s.n):
            initial[k]=np.array([[.7,.04j],[-.04j,.23]],complex)
    rho=np.array(initial,complex)
    if rho.shape!=(s.n,2,2):raise ValueError('bad rho')
    ans=np.zeros_like(rho)
    P=expm(-d['L']*t)
    for b in (0,1):ans[:,b,b]=P@rho[:,b,b]
    decay=np.exp(-(d['rate']+s.pure_phi)*t)
    ans[:,0,1]=rho[:,0,1]*decay*np.exp(1j*2*np.hypot(s.energy,s.mass)*t)
    ans[:,1,0]=rho[:,1,0]*decay*np.exp(-1j*2*np.hypot(s.energy,s.mass)*t)
    return ans


def equilibrium_B_collision_contact(s=CoherentSpec(),B=1e-3):
    """Static B Ward: frozen scattering FAILS to thermalize f(H0+B sigma_x).

    A B-dependent collision must obey D_B' F0 = -D0 F_B' to make
    f(H(B)) stationary. This is a NECESSARY equilibrium contact, not a
    constructed complete physically elastic dissipator at B != 0.
    """
    if not np.isfinite(B) or B<=0:raise ValueError('finite B positive')
    d=coherent_data(s)
    theta=d['theta'];U=d['u'];Ud=np.swapaxes(U.conj(),-1,-2)
    def f_B(x):
        # same fixed laboratory momentum, change only the external spin Hamiltonian
        H=s.energy*(np.cos(theta)[:,None,None]*SY-np.sin(theta)[:,None,None]*SX)+s.mass*SZ+x*SX
        ev,V=np.linalg.eigh(H)
        weights=fermi(ev,s.mu,s.temperature)
        FF=(V*weights[:,None,:])@np.swapaxes(V.conj(),-1,-2)
        return Ud@FF@U
    F0=d['rho0']
    FD_numeric=(f_B(B)-f_B(-B))/(2*B)
    spin=d['spin']
    f=fermi(d['e'],s.mu,s.temperature)
    de=d['e'][:,:,None]-d['e'][:,None,:]
    df=f[:,:,None]-f[:,None,:]
    ff=-f*(1-f)/s.temperature
    div=np.divide(df,de,out=np.zeros_like(de),where=np.abs(de)>1e-12)
    idx=np.arange(2)
    div[:,idx,idx]=ff
    FD=div*spin
    # Differentiate [H(B),f(H(B))]=0 in the unperturbed eigenbasis
    derivative_coherent=-1j*((d['e'][:,:,None]-d['e'][:,None,:])*FD +spin@F0-F0@spin)
    # Full generator on equilibrium derivative includes pure -i[H0,F_B']
    whole=apply_generator(s,FD,d)+vertex(spin,F0)
    scattered=whole-derivative_coherent
    # C0 f'_B is nonzero; B derivative of C_B must cancel it for
    # equilibrium preparation. Number trace must vanish.
    needed_contact=-scattered
    return dict(analytic_coherent_identity_residual=float(np.linalg.norm(derivative_coherent)),
        frozen_B_stationarity_defect=float(np.linalg.norm(whole)),
        required_collision_contact_norm=float(np.linalg.norm(needed_contact)),
        required_contact_total_trace=complex(np.sum(np.trace(needed_contact,axis1=-2,axis2=-1))),
        equilibrium_linear_B_derivative=FD,
        finite_B_derivative_residual=float(np.linalg.norm(FD_numeric-FD)),
        necessary_collision_contact=needed_contact,
        status='NECESSARY_EQ_B_COLLISION_VERTEX_NOT_CONSTRUCTED_FULL_B_DEPENDENT_BATH')


def _to_orb(d,r):return d['u']@r@np.swapaxes(d['u'].conj(),-1,-2)
def _from_orb(d,r):return np.swapaxes(d['u'].conj(),-1,-2)@r@d['u']

def orbital_dissipator_mass_field(s:CoherentSpec,B,rho_orb):
    """EXPLICIT GKSL scalar-impurity collision at spin Bz: H=Dirac+(m+B)sz.

    All momenta on the constant-radius shell remain exactly energy-degenerate
    at each B (unlike source spin Bx on a fixed-radius shell). Jumps are
    band-selective and use the TRUE field-dependent spinor overlap. This
    is a controlled Bz mass-field demonstration, not original Bx PHE.
    """
    if not np.isfinite(B):raise ValueError('finite mass-field B')
    dat=coherent_data(replace(s,mass=s.mass+B))
    r=_from_orb(dat,rho_orb)
    out=apply_generator(replace(s,mass=s.mass+B),r,dat)
    # remove coherent -i[H_B,rho], retain only collision D_B
    for i in range(s.n):
        eps=dat['e'][i]
        out[i]+=1j*(eps[:,None]-eps[None,:])*r[i]
    return _to_orb(dat,out)

def _linear_orb_B(s:CoherentSpec,B,vertex_e=SX,observable=SY):
    dat=coherent_data(replace(s,mass=s.mass+B))
    r0=dat['rho0']
    ev=_from_orb(dat,np.broadcast_to(vertex_e,(s.n,2,2)))
    current=_from_orb(dat,np.broadcast_to(observable,(s.n,2,2)))
    z=s.eta_e-1j*s.omega_e
    re=_resolvent_spectral(replace(s,mass=s.mass+B),vertex(ev,r0),z,dat)
    J=complex(np.mean(np.trace(current@re,axis1=-2,axis2=-1)))
    return J,_to_orb(dat,r0),_to_orb(dat,re)

def mass_zeeman_collision_vertex(s=CoherentSpec(),delta=5e-4,observable=SY):
    """Full Bz steady-state collision Ward and mixed E*B optical response.

    Analytic linearized Liouville in B with independently FIVE-POINT
    full B rediagonalization. Bz is an intentionally distinct controlled
    physical model, NOT the source Bx helicity, Dirac tau pole or PHE.
    """
    if not 0<delta<.02 or not np.isfinite(delta):raise ValueError('small positive delta')
    dat=coherent_data(s);z=s.eta_e-1j*s.omega_e
    J0,F0,RE0=_linear_orb_B(s,0,observable=observable)
    j2p,f2p,_=_linear_orb_B(s,2*delta,observable=observable)
    jp,fp,_=_linear_orb_B(s,delta,observable=observable)
    jm,fm,_=_linear_orb_B(s,-delta,observable=observable)
    j2m,f2m,_=_linear_orb_B(s,-2*delta,observable=observable)
    fd=(j2m-8*jm+8*jp-j2p)/(12*delta)
    Fprime=(f2m-8*fm+8*fp-f2p)/(12*delta)
    # Exact partial_B H is sigma_z in the orbital basis.
    def collision_diff(X):
        return (orbital_dissipator_mass_field(s,delta,X)-
                orbital_dissipator_mass_field(s,-delta,X))/(2*delta)
    dC_re=collision_diff(RE0)
    dC_eq=collision_diff(F0)
    C0_Fprime=orbital_dissipator_mass_field(s,0.,Fprime)
    eq_ward=dC_eq+C0_Fprime
    # Exact derivative of the Liouvillian L_B = -i[H_B,.]+C_B.
    zin=vertex(np.broadcast_to(SZ,(s.n,2,2)),RE0)+dC_re
    drive=vertex(np.broadcast_to(SX,(s.n,2,2)),Fprime)+zin
    re_deriv=_resolvent_spectral(s,_from_orb(dat,drive),z,dat)
    jband=_from_orb(dat,np.broadcast_to(observable,(s.n,2,2)))
    analytic=complex(np.mean(np.trace(jband@re_deriv,axis1=-2,axis2=-1)))
    no_vertex_drive=drive-dC_re
    re_wrong=_resolvent_spectral(s,_from_orb(dat,no_vertex_drive),z,dat)
    wrong=complex(np.mean(np.trace(jband@re_wrong,axis1=-2,axis2=-1)))
    return dict(linear_at_B0=J0,analytic_magnetic_derivative=analytic,
        independent_five_point=fd,full_difference_error=float(abs(analytic-fd)),
        frozen_collision_derivative=wrong,
        missing_collision_vertex_error=float(abs(wrong-fd)),
        collision_vertex_drive_norm=float(np.linalg.norm(dC_re)),
        equilibrium_collision_ward_residual=float(np.linalg.norm(eq_ward)),
        field_equilibrium_collision_defect_at_zero=float(np.linalg.norm(orbital_dissipator_mass_field(s,0,F0))),
        source='STATIC_SPIN_BZ_MASS_MODEL_NOT_ORIGINAL_BX_PHE')


def population_weighted_spectral(s:Spec,band,omega=0.,B=0.):
    """Independent detailed-balance eigenmode solver, not direct L inverse."""
    if band not in (-1,1) or not np.isfinite(omega):raise ValueError('band and omega')
    total=0j
    for E in s.energies:
        d=shell(s,E,band,B);w=d['weights'];r=np.sqrt(w)
        v=d['vy'];v0=v-np.dot(w,v)
        sym=r[:,None]*d['L']/r[None,:]
        ev,U=np.linalg.eigh((sym+sym.T)*.5)
        projected=U.T@(r*v0)
        value=np.divide(projected,ev-1j*omega,out=np.zeros_like(projected,dtype=complex),where=np.abs(ev)>1e-11)
        g=U@value/r
        f=fermi(band*E,s.mu,s.temperature)
        fac=f*(1-f)/s.temperature*E*np.mean(d['dos_angle'])
        total+=fac*np.dot(w,v*g)
    return complex(total)


def bx_fixed_radius_nonelastic(s:Spec,E=.3,B=.08):
    """Explicit rejection of fixed laboratory |k| for source spin Bx.

    True Bx isoenergy contours are translated in k_y. On a fixed lab ring,
    energies depend on angle at nonzero B and angular jumps are NOT elastic.
    """
    if E<=0 or not np.isfinite([E,B]).all():raise ValueError('positive E and finite B')
    th=2*np.pi*np.arange(s.n)/s.n
    kx=E*np.cos(th);ky=E*np.sin(th)
    h=s.beta*kx[:,None,None]*I+kx[:,None,None]*SY-(ky+B)[:,None,None]*SX
    e=np.linalg.eigvalsh(h)
    return dict(conduction_spread=float(np.ptp(e[:,1])),
                valence_spread=float(np.ptp(e[:,0])),
                fixed_radius_not_elastic=bool(np.ptp(e[:,1])>1e-5 and abs(B)>0),
                status='REJECT_FIXED_LAB_RING_FOR_SOURCE_BX_ELASTIC_SCATTERING')


def coherent_frame_covariance(s=CoherentSpec(),amplitude=.47):
    """Explicit nontrivial theta-dependent U(1)^2 band phase covariance."""
    if not np.isfinite(amplitude):raise ValueError('finite phase')
    d=coherent_data(s)
    theta=d['theta']
    phase=np.empty((s.n,2),complex)
    phase[:,0]=np.exp(1j*amplitude*np.sin(theta))
    phase[:,1]=np.exp(1j*amplitude*.7*np.cos(2*theta))
    V=np.zeros((s.n,2,2),complex)
    V[:,0,0]=phase[:,0];V[:,1,1]=phase[:,1]
    Ud=np.swapaxes(V.conj(),-1,-2)
    g={k:(Ud@d[k]@V) for k in ('ry','spin','vx','rho0')}
    ze=s.eta_e-1j*s.omega_e;zb=s.eta_b-1j*s.omega_b;zs=ze+zb
    e=_resolvent_spectral(s,vertex(g['ry'],g['rho0']),ze,d)
    b=_resolvent_spectral(s,vertex(g['spin'],g['rho0']),zb,d)
    be=_resolvent_spectral(s,vertex(g['ry'],b),zs,d)
    eb=_resolvent_spectral(s,vertex(g['spin'],e),zs,d)
    measure=lambda x:complex(np.mean(np.trace(g['vx']@x,axis1=-2,axis2=-1)))
    local=np.trace(g['vx']@(be+eb),axis1=-2,axis2=-1)
    return dict(b_then_e=measure(be),e_then_b=measure(eb),total=measure(be+eb),
       local_harmonic_1=complex(np.mean(np.exp(1j*theta)*local)))
