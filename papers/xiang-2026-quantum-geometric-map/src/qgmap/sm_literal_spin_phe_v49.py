"""v49 literal-source spin PHE audit: original main-text/16p-SM PDFs.

Two *distinct* published helicities and two *distinct* published derivations:
  SM Eq20 mixed-EB offdiagonal density times current operator,
  SM Eq31/32 symmetric Zeeman-curvature PHE.
Uses epsilon_mn = epsilon_m-epsilon_n and Bbar=-gmuB B; e=hbar=vF=tau=gmuB=1.
These literal source expressions have an unclosed SIGN difference. Do not silently
repair the paper, equate a finite-eta RTA to its adiabatic pole, or infer Fig2 units.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss

SX=np.array([[0.,1.],[1.,0.]],complex)
SY=np.array([[0.,-1j],[1j,0.]],complex)
ID=np.eye(2,dtype=complex)

@dataclass(frozen=True)
class SourceConfig:
    beta:float=.125
    mu:float=.3
    temperature:float=.03
    nr:int=100
    nphi:int=128
    cutoff:float=2.5
    def __post_init__(self):
        if not np.isfinite([self.beta,self.mu,self.temperature,self.cutoff]).all():
            raise ValueError('finite parameters required')
        if not 0<self.beta<1 or self.temperature<=0 or self.cutoff<=0:
            raise ValueError('0<beta<1, temperature>0, cutoff>0 required')
        if not isinstance(self.nr,int) or self.nr<20 or not isinstance(self.nphi,int) or self.nphi<32:
            raise ValueError('sufficient angular/radial quadrature required')


def geometry(beta):
    if not np.isfinite(beta) or not 0<beta<1:raise ValueError('0<beta<1')
    return -beta/(4*np.pi*(1+np.sqrt(1-beta*beta))**2)


def model(config=SourceConfig(),helicity='main',gauge_phases=(0.,0.)):
    """Evaluate literal Eq20 and Eq31/32 in one common, independently diagonalized basis.

    helicity=main   Eq10: +kx*sigma_y -ky*sigma_x
    helicity=sm     Eq62: -kx*sigma_y +ky*sigma_x
    Electric/current permutations xy and yx both computed. Unit output is
    conductivity/(tau*gmuB) with Bbar=-gmuB*B. No fitted signs.
    """
    if helicity not in ('main','sm'):raise ValueError('helicity=main or sm')
    if len(gauge_phases)!=2 or not np.isfinite(gauge_phases).all():raise ValueError('gauge phases')
    h=1 if helicity=='main' else -1
    q,w=leggauss(config.nr)
    r=config.cutoff*(1+q)*.5
    wr=w*config.cutoff*.5
    angles=2*np.pi*(np.arange(config.nphi)+.327)/config.nphi
    x=(r[:,None]*np.cos(angles)).ravel()
    y=(r[:,None]*np.sin(angles)).ravel()
    wt=np.broadcast_to((wr*r)[:,None]/(2*np.pi*config.nphi),(config.nr,config.nphi)).ravel()
    H=config.beta*x[:,None,None]*ID+h*(x[:,None,None]*SY-y[:,None,None]*SX)
    energies,U=np.linalg.eigh(H)
    U=U*np.exp(1j*np.asarray(gauge_phases))[None,None,:]
    Udag=U.conj().swapaxes(-1,-2)
    vx=config.beta*ID+h*SY
    vy=-h*SX
    Vx=Udag@vx@U
    Vy=Udag@vy@U
    Sx=Udag@SX@U
    f=1/(1+np.exp(np.clip((energies-config.mu)/config.temperature,-500,500)))
    fp=-f*(1-f)/config.temperature
    v={'x':Vx,'y':Vy}
    z={}
    for alpha,va in v.items():
        zz=np.zeros((len(x),2),float)
        for n,m in ((0,1),(1,0)):
            # Eq(31) uses v_nm=i(eps_n-eps_m)r_nm
            rnm=1j*va[:,n,m]/(energies[:,m]-energies[:,n])
            rmn=1j*va[:,m,n]/(energies[:,n]-energies[:,m])
            znm=1j*(rnm*Sx[:,m,n]-rmn*Sx[:,n,m])
            zz[:,n]=np.real(znm)
        z[alpha]=zz
    results={}
    for out,inp in (('x','y'),('y','x')):
        a=v[out]; b=v[inp]
        eq20_density_current=np.zeros(len(x),complex)
        eq31_curvature_current=np.zeros(len(x),float)
        for m,n in ((0,1),(1,0)):
            # SM Eq20 rho_mn/(tau E Bbar) = -sigma_mn * d_b(f_n-f_m)/eps_mn.
            grad_nm=fp[:,n]*b[:,n,n].real-fp[:,m]*b[:,m,m].real
            rho20=-Sx[:,m,n]*grad_nm/(energies[:,m]-energies[:,n])
            # Multiply rho_mn by v_nm and Bbar=-gmuB B (physical B).
            eq20_density_current+=-rho20*a[:,n,m]
            # Printed Eq31 after Bbar=-gmuB B gives +sum_n Z_nm d_b f_n.
            eq31_curvature_current+=z[out][:,n]*fp[:,n]*b[:,n,n].real
        literal20=float(np.dot(wt,np.real(eq20_density_current)))
        literal31=float(np.dot(wt,eq31_curvature_current))
        results[out+inp]={'SM_Eq20_current_physical_B':literal20,
                          'SM_printed_Eq31_current_physical_B':literal31,
                          'Eq20_plus_Eq31_residual':literal20+literal31}
    eq32_direct=.5*(results['xy']['SM_printed_Eq31_current_physical_B']+
                    results['yx']['SM_printed_Eq31_current_physical_B'])
    literal_formula=geometry(config.beta)*np.tanh(config.mu/(2*config.temperature))
    return {'protocol':'v49_source_literal_SM_Eq20_vs_Eq31_Eq32_main_Eq10_SM_Eq62',
            'scientific_status':'PARTIAL',
            'source_pdf_verified':True,'helicity':helicity,
            'beta':config.beta,'mu':config.mu,'temperature':config.temperature,
            'cutoff':config.cutoff,'nr':config.nr,'nphi':config.nphi,
            'xy':results['xy'],'yx':results['yx'],
            'SM_printed_Eq32_symmetric':eq32_direct,
            'main_Eq11_geometry_tanh_extension':literal_formula,
            'Eq32_to_Eq11_ratio':eq32_direct/literal_formula if abs(literal_formula)>1e-13 else None,
            'warning':['SM Eq20 and printed Eq31 give opposite physical B signs using the printed velocity identity',
                       'Main Eq10 and SM Eq62 have opposite spin helicity at fixed physical sigma_x',
                       'Tanh continuation is a model result; author Eq11 is at T=0',
                       'Author Eq15 regulates a selected electric-field adiabatic pole, not uniform eta on all matrix channels',
                       'Original full physics and Fig2 published SI ordinate remain PARTIAL/OPEN']}


def compare_source(config=SourceConfig()):
    main=model(config,'main')
    sm=model(config,'sm')
    return {'protocol':'v49_authentic_uploaded_main_SM_source_audit',
            'scientific_status':'PARTIAL',
            'main':main,'supplementary':sm,
            'main_to_SM_helicity_curvature_sum':main['SM_printed_Eq32_symmetric']+sm['SM_printed_Eq32_symmetric'],
            'source_conclusion':'Literal SM Eq20/31 and main/SM Hamiltonians do not have a unique consistent signed mapping without explicit convention changes.'}
