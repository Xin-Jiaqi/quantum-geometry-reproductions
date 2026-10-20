"""SM Eqs.(48),(49),(52) orbital intrinsic MNHE, genuine 3D two-band test.

Keep the *literal* Eq.52 factor of 1 beside the Eq.49-derived integrated
by parts result (factor 2). Numeric values are dimensionless (e=hbar=1),
BZ=[-pi,pi)^3. This evaluates the *two-band terms retained by the authors*
but not omitted three-band processes or a microscopic disorder/bath protocol.
"""
from __future__ import annotations
import numpy as np
from ._spin_grid import fermi

EPS=np.zeros((3,3,3),int)
for abc in ((0,1,2),(1,2,0),(2,0,1)):EPS[abc]=1
for abc in ((0,2,1),(2,1,0),(1,0,2)):EPS[abc]=-1


def gapped_3d(N=41, offset=(.19,.27,.31),tilt_scale=1.):
    if N<9:raise ValueError('N>=9 required for 3D differences')
    k=[(np.arange(N)+off)*2*np.pi/N-np.pi for off in offset]
    x,y,z=np.meshgrid(*k,indexing='ij')
    d=np.stack([np.sin(x)+.23*np.cos(y)+.11*np.sin(z),
                np.sin(y)+.15*np.cos(z)+.06*np.sin(x),
                2.7+.32*np.cos(x)+.28*np.cos(y)+.3*np.cos(z)
                  +.11*np.sin(x+y)+.07*np.sin(x-z)])
    dd=np.stack([np.stack([np.cos(x),.06*np.cos(x),-.32*np.sin(x)+.11*np.cos(x+y)+.07*np.cos(x-z)]),
                 np.stack([-.23*np.sin(y),np.cos(y),-.28*np.sin(y)+.11*np.cos(x+y)]),
                 np.stack([.11*np.cos(z),-.15*np.sin(z),-.3*np.sin(z)-.07*np.cos(x-z)])])
    d0=.33*np.sin(x)+.24*np.sin(y)+.18*np.sin(z)+.12*np.cos(x-y)
    d0der=np.stack([.33*np.cos(x)-.12*np.sin(x-y),
                     .24*np.cos(y)+.12*np.sin(x-y), .18*np.cos(z)])
    d0*=tilt_scale
    d0der*=tilt_scale
    D=np.linalg.norm(d,axis=0)
    dot=lambda a,b:np.sum(a*b,axis=0)
    dr=np.stack([dot(d,dd[i]) for i in range(3)])
    metric=np.empty((3,3,*D.shape))
    for i in range(3):
        for j in range(3):
            metric[i,j]=(dot(dd[i],dd[j])-dr[i]*dr[j]/D**2)/(4*D**2)
    E=np.stack([d0-D,d0+D])
    vel=np.stack([d0der-dr/D,d0der+dr/D])
    kxy=dot(d,np.cross(dd[0],dd[1],axis=0))
    Omega_xy=np.stack([kxy/(2*D**3),-kxy/(2*D**3)])
    moment_z=-kxy/(2*D**2)
    if D.min()<1e-4:raise ValueError('non-Abelian treatment needed')
    return {'E':E,'velocity':vel,'metric':metric,'Omega_xy':Omega_xy,
            'moment_z':moment_z,'min_direct_gap':float(2*D.min()),
            'd0der':d0der,'D':D, 'N':N}


def periodic_derivative(a,axis,dx):
    """Central periodic derivative; axis is one of the 3 momentum directions."""
    ax=a.ndim-3+axis
    return (np.roll(a,-1,axis=ax)-np.roll(a,1,axis=ax))/(2*dx)


def orbital_mnhe_3d(N=41,mu=1.6,T=.27,Bstep=0.0005,tilt_scale=1.,a=0,b=1):
    if T<=0 or Bstep<=0:raise ValueError('positive T and Bstep required')
    if (a,b) not in ((0,1),(1,0)):raise ValueError('only xy/yx with Bz currently implemented')
    p=gapped_3d(N,tilt_scale=tilt_scale);dk=2*np.pi/N
    metric=p['metric']; v=p['velocity'];E=p['E'];D=p['D'];d0der=p['d0der'];omega=p['Omega_xy'];m=p['moment_z']
    f,fp=fermi(E,mu,T)
    dg=np.stack([periodic_derivative(metric,axis=i,dx=dk) for i in range(3)])
    c=2
    F=np.zeros((2,3,N,N,N))
    geom_pref=np.zeros((2,3,N,N,N))
    for n,sign in enumerate((-1,1)):
        for direction in range(3):
            for i in range(3):
                for j in range(3):
                    t=EPS[direction,i,j]
                    if not t:continue
                    piece=sign*d0der[i]/D*metric[j,c]
                    geom_pref[n,direction] += t*piece
                    F[n,direction] -= t*(2*piece +dg[i,j,c])
    curl=periodic_derivative(F[:,b],a,dk)-periodic_derivative(F[:,a],b,dk)
    term_curl=float(np.mean(np.sum(curl*f,axis=0)))
    term_sbc=float(np.mean(np.sum(-EPS[a,b,c]*omega*m*fp,axis=0)))
    geometric_integrand=(v[:,a]*geom_pref[:,b]-v[:,b]*geom_pref[:,a])*fp
    geometric=float(np.mean(np.sum(geometric_integrand,axis=0)))
    metric_second=np.zeros((N,N,N))
    for i in range(3):
        for j in range(3):
            if EPS[b,i,j]:metric_second-=EPS[b,i,j]*periodic_derivative(dg[i,j,c],a,dk)
            if EPS[a,i,j]:metric_second+=EPS[a,i,j]*periodic_derivative(dg[i,j,c],b,dk)
    metric_term=float(np.mean(np.sum(f*metric_second,axis=0)))
    literal=geometric+metric_term+term_sbc
    ibp=2*geometric+metric_term+term_sbc
    direct=term_curl+term_sbc
    ibp_exact_integrand=(-F[:,b]*v[:,a]+F[:,a]*v[:,b])*fp
    ibp_direct=float(np.mean(np.sum(ibp_exact_integrand,axis=0)))+term_sbc
    def hall(B):
        occ,_=fermi(E-B*m,mu,T)
        return float(np.mean(np.sum((EPS[a,b,c]*omega+B*curl)*occ,axis=0)))
    finite_field=(hall(Bstep)-hall(-Bstep))/(2*Bstep)
    return {'protocol':'3D_2band_orbital_intrinsic_SM48_49_vs_literal_SM52',
        'N':N,'mu':mu,'T':T,'Bstep':Bstep,'tilt_scale':tilt_scale,'a':'xy'[a],'b':'xy'[b],'min_direct_gap':p['min_direct_gap'],
        'eq48_positional_curl':term_curl,'eq48_orbital_moment':term_sbc,
        'eq48_total':direct,'eq48_ibp_direct':ibp_direct,
        'eq52_printed_literal':literal,'eq52_eq49_derived_factor2':ibp,
        'factor_one_geometric_piece':geometric,'second_derivative_metric_piece':metric_term,
        'eq52_printed_minus_eq48':literal-direct,
        'eq49_derived_minus_eq48':ibp-direct,
        'finite_field_wavepacket_derivative':finite_field,
        'finite_field_minus_eq48':finite_field-direct,
        'units':'e=hbar=1; 3D BZ mean; B couples orbital moment in normalized energy units',
        'limits':'Eq49 and Eq52 omit three-band processes; finite-B check uses Eq49 positional shift'}


def orbital_mnhe_2d(N=141,mu=1.15,T=.10,Bstep=.0003):
    """Strictly 2D B_z: g^{iz}=0, so positional F^{bz} vanishes."""
    from .orbital_phe_v4 import periodic_data
    if Bstep<=0 or T<=0:raise ValueError('positive Bstep and T')
    q=periodic_data(N)
    E=q['E'];omega=q['Omega'];m=q['m'];f,fp=fermi(E,mu,T)
    sbc=float(np.mean(np.sum(-omega*m*fp,axis=-1)))
    eq52=float(np.mean(np.sum(-.5*(E-np.flip(E,axis=-1))*omega**2*fp,axis=-1)))
    def ah(B):
        ff,_=fermi(E-B*m,mu,T)
        return float(np.mean(np.sum(omega*ff,axis=-1)))
    fd=(ah(Bstep)-ah(-Bstep))/(2*Bstep)
    return dict(protocol='strict_2D_Bz_orbital_MNHE_SBC_only',N=N,mu=mu,T=T,
      Bstep=Bstep,min_direct_gap=q['minimum_direct_gap'],
      eq48_moment_Fermi_surface=sbc,eq52_squared_curvature=eq52,
      independent_B_derivative=fd,eq52_minus_eq48=eq52-sbc,finite_field_minus_eq48=fd-sbc,
      limit='No z momentum dependence, F_xz=F_yz=0. This does not test QMQ.')
