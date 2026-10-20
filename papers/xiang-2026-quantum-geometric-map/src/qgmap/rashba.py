"""T=0 Rashba Eq.(9) OHE, e=ℏ=1, constant τ; distinct Fermi-ring sectors."""
import math

def zero_temperature(m,lam,mu):
    if m<=0 or lam<=0:raise ValueError("m and lambda must be positive")
    floor=-m*lam*lam/2
    if abs(mu)<1e-13 or abs(mu-floor)<1e-13:
        raise ValueError("Crossing and band-edge order-of-limits not implemented")
    if mu<floor:return dict(regime="empty",bare=0.,metric=0.,total=0.,radii=[])
    v=math.sqrt(lam*lam+2*mu/m)
    if mu<0:
        bare=lam*v/(2*math.pi)
        return dict(regime="lower-band-two-rings",bare=bare,metric=-bare,total=0.,
                    radii=[m*(lam-v),m*(lam+v)])
    bare=v*v/(2*math.pi)
    return dict(regime="two-helicity-bands",bare=bare,metric=0.,total=bare,
                radii=[m*(v-lam),m*(v+lam)])

def direct_fermi_ring_sum(m,lam,mu):
    """Independent delta(ε-μ) evaluation; verifies analytic sector decomposition."""
    if m<=0 or lam<=0:raise ValueError("m and lambda must be positive")
    floor=-m*lam*lam/2
    if abs(mu)<1e-13 or abs(mu-floor)<1e-13:raise ValueError("Singular threshold")
    if mu<floor:return dict(bare=0.,metric=0.,total=0.,roots=[])
    v=math.sqrt(lam*lam+2*mu/m)
    bare=metric=0.;roots=[]
    for s in (-1,1):
        for k in (-m*s*lam-m*v,-m*s*lam+m*v):
            if k<=1e-12:continue
            velocity=k/m+s*lam
            bare+=k*velocity**2/(abs(velocity)*4*math.pi*m)
            metric+=s*lam*velocity**2/(abs(velocity)*4*math.pi)
            roots.append((s,k,velocity))
    return dict(bare=bare,metric=metric,total=bare+metric,roots=roots)
