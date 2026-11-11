"""Stage 7j: gauge-covariant physical overlaps in mixed (k, spin-source) space.

A changing Wannier frame U(k,lambda) is a REPRESENTATION change only.
For occupied eigenvector coordinates y in that frame, the physical
Wilson link is y_a^dag U_a^dag U_b y_b, NOT y_a^dag y_b.

The source gauge used here is deliberately artificial and is a TEST
of operator completeness, not an ab-initio Wannierization or prediction.
Physical scope: gapped synthetic 8-band T=0 equilibrium/static NMEE only.
"""
import numpy as np
from .sliding_tight_binding import SlidingBilayer, SlidingSpec
from nmee_lab.clean_core.core import MU_B

# Artificial gauge phase lengths in nm, lambda phase slopes in eV^-1.
# Use spin-equal spatial lengths but unequal lambda slopes.
ORBITAL_X=np.array([.8,.8,-.5,-.5,.35,.35,-.3,-.3])
ORBITAL_Y=np.array([.25,.25,.7,.7,-.45,-.45,.5,.5])
SOURCE_SLOPE=np.array([20.,-15.,9.,-25.,12.,-7.,-18.,24.])


def embedding_unitary(k,lam=0.):
    k=np.asarray(k,float)
    if k.shape!=(2,):
        raise ValueError("Single Cartesian k=(kx,ky) required")
    return np.diag(np.exp(1j*(ORBITAL_X*k[0] +
                            ORBITAL_Y*k[1]+SOURCE_SLOPE*lam)))


def physical_link(ya,yb,Ua,Ub):
    """Link of occupied frames written in *different* orbital bases."""
    return ya.conj().T @ (Ua.conj().T@Ub) @ yb


def wilson_phase_at_field(model,k,field_T,dk,dl,*,gauge="original",
                          spin_axis="x",k_axis=0):
    if gauge not in ("original","naive","covariant"):
        raise ValueError("gauge must be original/naive/covariant")
    if spin_axis not in ("x","y") or k_axis not in (0,1):
        raise ValueError("Only x/y axes supported")
    if min(dk,dl)<=0:raise ValueError("positive plaquette source steps required")
    k=np.asarray(k,float)
    offset=np.eye(2)[k_axis]*dk/2
    corners=((k-offset,-dl/2),(k+offset,-dl/2),
             (k+offset,+dl/2),(k-offset,+dl/2))
    frames=[];Uframes=[];min_gap=float("inf")
    for kk,lam in corners:
        H=model.matrices(kk[None,:])[0]
        source=model.spins["total_"+spin_axis]
        h=H+lam*source-field_T*model.mono.g_factor*MU_B*model.spins["total_z"]
        eig,U=np.linalg.eigh(h)
        min_gap=min(min_gap,float(min(
            eig[4]-model.spec.mu_eV,model.spec.mu_eV-eig[3])))
        if min_gap<=0:raise ValueError("Occupied gap closes")
        unit=embedding_unitary(kk,lam) if gauge!="original" else np.eye(8)
        frames.append(unit.conj().T@U[:,:4])
        Uframes.append(unit)
    product=1.+0j
    for i in range(4):
        j=(i+1)%4
        if gauge=="covariant":
            link=physical_link(frames[i],frames[j],Uframes[i],Uframes[j])
        else:
            link=frames[i].conj().T@frames[j]
        det=np.linalg.det(link)
        if abs(det)<1e-10:raise ArithmeticError("Ill-conditioned occupied link")
        product*=det/abs(det)
    return np.angle(product),min_gap


def local_static_gamma(model,k,*,gauge="original",dk=.006,
                       dl=.00015,B_step_T=.25,spin_axis="x",k_axis=0):
    """Equilibrium STATIC dB mixed Berry curvature at a local k point."""
    if B_step_T<=0:raise ValueError("positive field step required")
    args=dict(gauge=gauge,dk=dk,dl=dl,spin_axis=spin_axis,k_axis=k_axis)
    plus,gp=wilson_phase_at_field(model,k,+B_step_T,**args)
    minus,gm=wilson_phase_at_field(model,k,-B_step_T,**args)
    result=np.angle(np.exp(1j*(plus-minus)))/(2*B_step_T*dk*dl)
    return result,min(gp,gm)


def covariant_source_derivative(model,k,axis="kx",step=1e-5):
    """Compare naive and connection-corrected source H derivatives.

    With Omega_a=U^dag d_a U, D_a H'=d_a H'+[Omega_a,H']
    equals U^dag (d_a H) U for a=kx,ky,lambda.
    """
    if axis not in ("kx","ky","lambda") or step<=0:
        raise ValueError("Unsupported source or step")
    k=np.asarray(k,float)
    H=model.matrices(k[None,:])[0]
    if axis=="lambda":
        slope=SOURCE_SLOPE
        def transform(s):
            U=embedding_unitary(k,s)
            return U.conj().T@(H+s*model.spins["total_x"])@U
        target=model.spins["total_x"]
    else:
        j=0 if axis=="kx" else 1
        slope=ORBITAL_X if j==0 else ORBITAL_Y
        def transform(s):
            kk=k+np.eye(2)[j]*s
            U=embedding_unitary(kk,0.)
            return U.conj().T@model.matrices(kk[None,:])[0]@U
        target=model.matrices(k[None,:],axis=j)[0]
    U0=embedding_unitary(k,0.)
    hp=U0.conj().T@H@U0
    derivative=(transform(step)-transform(-step))/(2*step)
    omega=1j*np.diag(slope)
    corrected=derivative+omega@hp-hp@omega
    physical=U0.conj().T@target@U0
    return dict(naive_error=float(np.max(abs(derivative-physical))),
                covariant_error=float(np.max(abs(corrected-physical))))


def gauge_stress_case(spec:SlidingSpec,k,*,dk=.006,dl=.00015):
    model=SlidingBilayer(spec)
    original,g0=local_static_gamma(model,k,dk=dk,dl=dl)
    corrected,g1=local_static_gamma(model,k,gauge="covariant",dk=dk,dl=dl)
    naive,g2=local_static_gamma(model,k,gauge="naive",dk=dk,dl=dl)
    deriv={axis:covariant_source_derivative(
        model,k,axis,step=1e-6 if axis=="lambda" else 1e-5)
           for axis in ("kx","ky","lambda")}
    return {"scope":"8band T=0 gapped STATIC mixed Berry curvature at one k",
            "u_fractional":list(spec.shift),"t_eV":spec.tunneling_eV,
            "k_inv_nm":list(np.asarray(k,float)),
            "original":float(original),"corrected":float(corrected),
            "naive":float(naive),
            "physical_gauge_error":float(abs(original-corrected)),
            "naive_artifact":float(naive-original),
            "minimum_mu_gap_eV":float(min(g0,g1,g2)),
            "derivative_errors":deriv,
            "warning":"u(k,lambda) is a pure synthetic gauge, not a new physical perturbation"}
