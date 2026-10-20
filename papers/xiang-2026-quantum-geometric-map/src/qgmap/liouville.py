"""Finite-dimensional two-color E/B Liouville resolvent.

ħ=1; E/B are independent explicit Hermitian operators, NOT an automatically
gauge-complete uniform electric field in a Bloch band. These are operator
response kernels, not universal physical DC conductivities.
"""
import numpy as np

def _validate(h,v_e,v_b,rho):
    n=h.shape[0]
    for name,x in [("H",h),("VE",v_e),("VB",v_b),("rho",rho)]:
        if x.shape!=(n,n) or not np.allclose(x,x.conj().T,atol=1e-10):
            raise ValueError(f"{name} must be Hermitian n×n")
    if not np.isfinite(h).all():raise ValueError("Nonfinite Hamiltonian")

def _comm(a,b):
    return a@b-b@a

def _resolvent(e,u,source,omega,eta):
    if eta<=0:raise ValueError("Positive eta required")
    s=u.conj().T@source@u
    d=e[:,None]-e[None,:]
    return u@(s/(omega+1j*eta-d))@u.conj().T

def two_color_response(H,V_E,V_B,rho0,omega_e=0.,omega_b=0.,
                       eta_e=.02,eta_b=.02,eta_total=None):
    """R_{ωE+ωB}([VE,R_B[VB,ρ0]]+[VB,R_E[VE,ρ0]]).

    Return both chronologies; they must not be silently mixed.
    Total broadening defaults to η_E+η_B (adiabatic switching convention).
    """
    h=np.asarray(H,complex);a=np.asarray(V_E,complex)
    b=np.asarray(V_B,complex);r=np.asarray(rho0,complex)
    _validate(h,a,b,r)
    e,u=np.linalg.eigh(h)
    rg=lambda source,omega,eta:_resolvent(e,u,source,omega,eta)
    r_e=rg(_comm(a,r),omega_e,eta_e)
    r_b=rg(_comm(b,r),omega_b,eta_b)
    total_eta=eta_e+eta_b if eta_total is None else eta_total
    be=rg(_comm(a,r_b),omega_e+omega_b,total_eta)
    eb=rg(_comm(b,r_e),omega_e+omega_b,total_eta)
    return {"B_then_E":be,"E_then_B":eb,"total":be+eb,
            "rho_E":r_e,"rho_B":r_b,"eta_total":float(total_eta)}

def equilibrium_population_derivative(energies,diagonal_perturbation,mu,kbt):
    """Separate thermodynamic f'(ε) V_nn population term; not in adiabatic kernel."""
    if kbt<=0:raise ValueError("Finite kBT required")
    e=np.asarray(energies,float);v=np.asarray(diagonal_perturbation,float)
    if e.shape!=v.shape:raise ValueError("Energy / V diagonal shape mismatch")
    z=np.clip((e-mu)/(2*kbt),-300,300)
    return -v/(4*kbt*np.cosh(z)**2)
