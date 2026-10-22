"""B+D Gate 6: re-audit Gate 5 with EXCHANGE vs atomic-SOC torque.

Execute from scripts/ (depends ONLY on existing Gate4A, Gate4B and numpy):
  python bd_gate6_exchange_soc_torque.py
  python bd_gate6_exchange_soc_torque.py --dense

Independent SYNTHETIC orbital honeycomb bilayer, layer-A-type AFM;
NOT real MnPS3 honeycomb-site AFM, not DFT and not a new effect.
No jobs are submitted. e^{-i omega t}, eta=0, subgap.
"""
import json
import sys
import numpy as np
from bd_gate4_polar_honeycomb import B1, B2, fraction_delta
from bd_gate4b_atomic_spin import (
    model, K, I, I4, X, Y, Z, LX, LY, LZ, STAG_x, Tspin_y,
)

EXCHANGE_UNIT = K(Z, I, I4, Z)
SOC_UNIT = (K(I, I, LX, X) + K(I, I, LY, Y) + K(I, I, LZ, Z)) / 2

def susceptibility(frac=(1/3,1/3), N=11, omega=.06,
                   m=.34, soc=.30, out_scale=1.):
    shift = fraction_delta(frac)
    vals = np.zeros(7, dtype=complex)
    minimum_direct_gap = float("inf")
    max_operator_error = 0.
    for u in range(N):
        for v in range(N):
            k = (u*B1 + v*B2)/N
            H, velocity = model(k,shift,m=m,soc=soc,out_scale=out_scale)
            H_disconnected, _ = model(k,shift,m=m,soc=soc,out_scale=0.)
            H_inter = H - H_disconnected
            H_ex = m*EXCHANGE_UNIT
            H_soc = soc*SOC_UNIT
            H_kin_local = H - H_inter - H_ex - H_soc

            def torque(Hpart):
                return 1j*(Hpart@STAG_x-STAG_x@Hpart)

            Ointer = torque(H_inter)
            Ointra = torque(H_kin_local)
            Oex = torque(H_ex)
            Osoc = torque(H_soc)
            Ototal = torque(H)
            # For H_ex=m tau_z sigma_z and S_-^x=tau_z sigma_x/2:
            # i[H_ex,S_-^x]=-m sigma_y = -2m S_tot^y.
            max_operator_error=max(
                max_operator_error,
                np.linalg.norm(Ototal-Ointer-Ointra-Oex-Osoc),
                np.linalg.norm(Oex + 2*m*Tspin_y),
                np.linalg.norm(Ointra),
            )
            assert max_operator_error<1e-11

            # Generalized force for a *rigid, hypothetical A-type* Neel angle:
            # -dH_ex/dtheta=-m tau_z sigma_x=-2m S_-^x.
            Ftheta=-2*m*STAG_x
            eig,U=np.linalg.eigh(H)
            diff=eig[8:][None,:]-eig[:8,None]
            minimum_direct_gap=min(minimum_direct_gap,float(np.min(diff)))
            vel=U.conj().T@velocity@U
            for t,op in enumerate((STAG_x,Ointer,Ointra,Oex,Osoc,Ototal,Ftheta)):
                mat=U.conj().T@op@U
                r=-1j*mat[:8,8:]*vel[8:,:8].T/diff
                vals[t]+=np.sum(r/(omega-diff)-r.conj()/(omega+diff))/N**2
    return vals,minimum_direct_gap,max_operator_error

def test(N=11,omega=.06):
    a,gap,err=susceptibility(N=N,omega=omega)
    b,_,_=susceptibility(frac=(-1/3,-1/3),N=N,omega=omega)
    if not 0<=omega<gap: raise AssertionError("Not subgap")
    # Paper-derived MnPS3 L2 Eq13 is NOT used; these are only toy m_zT partners.
    assert np.max(np.abs(a[:6]+b[:6]))<1e-11
    assert abs(np.sum(a[1:5])-a[5])<1e-11
    assert abs(a[5]+1j*omega*a[0])<1e-11
    assert abs(a[6]+2*.34*a[0])<1e-11
    assert abs(a[2])<1e-12
    for kwargs in (dict(out_scale=0.), dict(m=0.), dict(soc=0.)):
        blank,_,_=susceptibility(N=7,omega=omega,**kwargs)
        assert np.max(np.abs(blank[:6]))<1e-10
    return {"mesh":N,"omega":omega,"min_direct_gap":gap,
            "chi_stag_spin":str(a[0]),"interlayer":str(a[1]),
            "exchange_electronic":str(a[3]),"SOC_electronic":str(a[4]),
            "net":str(a[5]),"Neel_generalized_force":str(a[6]),
            "exchange_inter_abs_ratio":float(abs(a[3]/a[1])),
            "force_inter_abs_ratio":float(abs(a[6]/a[1])),
            "ward_residual":float(abs(a[5]+1j*omega*a[0])),
            "operator_residual":float(err)}

def magnetic_texture_counterexample():
    """Four Mn sites illustrative A/B in top and bottom layers.
    Do NOT interpret a sample pattern as author verified L2 indexing.
    """
    layer_odd=np.array([1,1,-1,-1])
    for honeycomb_neel in ((1,-1,1,-1),(1,-1,-1,1)):
        eta=np.array(honeycomb_neel)
        assert eta.sum()==0
        assert np.dot(eta,layer_odd)==0
        assert not np.array_equal(eta,layer_odd)
    print("PASS: honeycomb Néel and layer-odd spin are distinct operators")

if __name__=="__main__":
    magnetic_texture_counterexample()
    dense="--dense" in sys.argv
    for N in ((11,19,23) if dense else (7,11)):
        print(json.dumps(test(N),indent=2))
    print("MODEL DECOMPOSITION PASS; MAGNETIC SWITCHING AND NOVELTY NOT ESTABLISHED")
