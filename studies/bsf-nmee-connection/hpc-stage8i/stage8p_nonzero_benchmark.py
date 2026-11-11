"""Stage8p: nonzero, exactly solvable 2-band mixed-geometry benchmark.

This is a numerical METHOD VALIDATION, not a WSe2 material model:
H(k,lambda,B) = sin(k)*sigma_y + (lambda/2)*sigma_x
                 + (mass - beta*B)*sigma_z
Source S_x = sigma_x/2 and magnetic derivative dH/dB=-beta sigma_z.
The independent nonzero analytic, resolvent/projector and determinant Wilson
answers must agree, including signs and independent B central differences.
No VASP, no external input, no GW, no new PBS jobs.
"""
import argparse
import json
import numpy as np

SIGX=np.array([[0,1],[1,0]],complex)
SIGY=np.array([[0,-1j],[1j,0]],complex)
SIGZ=np.diag([1,-1]).astype(complex)


def ham(k,lam,B,mass=1.2,beta=.4):
    return np.sin(k)*SIGY + .5*lam*SIGX + (mass-beta*B)*SIGZ


def occ_frame(k,lam,B):
    _,v=np.linalg.eigh(ham(k,lam,B))
    return v[:,[0]]


def occ_projector(k,lam,B):
    v=occ_frame(k,lam,B)
    return v@v.conj().T


def geom(P,dlP,dkP):
    return float((1j*np.trace(P@(dlP@dkP-dkP@dlP))).real)


def analytic_F(k,B,mass=1.2,beta=.4):
    z=mass-beta*B
    return float(.25*z*np.cos(k)/(np.sin(k)**2+z*z)**1.5)


def analytic_Gamma(k,mass=1.2,beta=.4):
    s=np.sin(k)
    return float(-.25*beta*np.cos(k)*(s*s-2*mass*mass)/(s*s+mass*mass)**2.5)


def source_derivative_resolvent(k,B):
    E,V=np.linalg.eigh(ham(k,0.,B))
    S=V.conj().T@(.5*SIGX)@V
    gap=E[0]-E[1]
    dP=np.zeros((2,2),complex)
    dP[1,0]=S[1,0]/gap
    dP[0,1]=S[0,1]/gap
    return V@dP@V.conj().T


def F_resolvent(k,B,dk=2e-4):
    P=occ_projector(k,0.,B)
    dPk=(occ_projector(k+dk,0.,B)-occ_projector(k-dk,0.,B))/(2*dk)
    return geom(P,source_derivative_resolvent(k,B),dPk)


def F_projector(k,B,dk=2e-4,dl=2e-4):
    P=occ_projector(k,0.,B)
    dPl=(occ_projector(k,dl,B)-occ_projector(k,-dl,B))/(2*dl)
    dPk=(occ_projector(k+dk,0.,B)-occ_projector(k-dk,0.,B))/(2*dk)
    return geom(P,dPl,dPk)


def wilson_phase(k,B,dk=2e-3,dl=2e-3,phase_shift=None):
    pts=[(k-dk/2,-dl/2),(k+dk/2,-dl/2),
         (k+dk/2,+dl/2),(k-dk/2,+dl/2)]
    fs=[]
    for idx,(kn,ln) in enumerate(pts):
        f=occ_frame(kn,ln,B)
        if phase_shift is not None:
            f=f*np.exp(1j*float(phase_shift[idx]))
        fs.append(f)
    prod=1.+0j
    for n in range(4):
        prod*=np.vdot(fs[n][:,0],fs[(n+1)%4][:,0])
    return float(np.angle(prod))


def F_wilson(k,B,dk=2e-3,dl=2e-3,phase_shift=None):
    return wilson_phase(k,B,dk,dl,phase_shift)/(dk*dl)


def gamma_centered(F,k,dB=.01):
    return (F(k,dB)-F(k,-dB))/(2*dB)


def validation(k=.4):
    exact=analytic_Gamma(k)
    values={
        "analytic":exact,
        "resolvent":gamma_centered(F_resolvent,k),
        "finite_projector":gamma_centered(F_projector,k),
        "wilson":gamma_centered(F_wilson,k),
    }
    err={key:abs(float(v-exact))/abs(exact) for key,v in values.items()}
    if abs(exact)<.005 or max(err.values())>5e-4:
        raise AssertionError(f"Nonzero static benchmark failed: {values}; rel={err}")
    # Gauge invariance: four independently rephased Bloch eigenvectors.
    phase=[.12,-.73,2.3,1.02]
    w0=F_wilson(k,0)
    wg=F_wilson(k,0,phase_shift=phase)
    if not np.isclose(w0,wg,atol=5e-12,rtol=1e-9):
        raise AssertionError("Wilson loop depends on arbitrary band phases")
    return {
        "stage":"8p","kind":"nonzero_synthetic_method_calibration_NOT_material",
        "physical_nmee":None,"new_PBS_jobs":0,
        "model":"sin(k)*sigma_y + lambda/2*sigma_x + (mass-beta*B)*sigma_z",
        "k":k,"mass_eV":1.2,"beta_eV_per_B_unit":.4,
        "gamma_dB_unit":values,
        "relative_errors":err,
        "all_three_nonzero_methods_agree":True,
        "warning":"This 2-band calibration model is not spinful WSe2, nor generally time-reversal invariant. It tests sign/factors only."
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",help="optional new JSON path; refuse overwrite")
    args=parser.parse_args()
    result=validation()
    if args.out:
        from pathlib import Path
        target=Path(args.out)
        if target.exists():raise FileExistsError("Refusing to overwrite")
        target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
