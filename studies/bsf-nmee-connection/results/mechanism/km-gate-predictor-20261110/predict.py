"""Analytic gate tangent of the existing eight-state insulating KM response.

No fitted bands or altered response definition. U is upper-minus-lower energy
(eV); observables and velocities are held in the original embedded basis.
The denominator/vertex pieces are derivative diagnostics, not new materials.
"""
from pathlib import Path
import sys, os, json, argparse, time
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np

ROOT=Path(__file__).resolve().parent
candidates=list(ROOT.parents)+[p/'quantum-geometry-reproductions' for p in ROOT.parents]
REPO=next(p for p in candidates if (p/'papers/jia-2026-nonlin-magnetoelectric-edelstein').is_dir())
STUDY=REPO/'studies/bsf-nmee-connection'
sys.path[:0]=[str(STUDY),str(REPO/'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from bsf_nmee.bilayer_audit import AuditBilayer,AuditSpec,integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift
from bsf_nmee.sliding_tight_binding import P_TOP,P_BOTTOM
from nmee_lab.clean_core.core import MU_B

adj=lambda a:a.conj().swapaxes(-1,-2)
W=(P_TOP-P_BOTTOM)/2

def kernel(h,v,spins,derivative=True,phase_seed=None):
    e,u=np.linalg.eigh(h)
    if phase_seed is not None:
        rng=np.random.default_rng(phase_seed);u=u*np.exp(1j*rng.uniform(-np.pi,np.pi,e.shape))[:,None,:]
    vb=adj(u)@v@u;wb=adj(u)@W@u
    gaps=e[:,:4,None]-e[:,None,4:]
    same_sep=np.concatenate([np.diff(e[:,:4],axis=1),np.diff(e[:,4:],axis=1)],axis=1)
    spacing=float(np.min(same_sep))
    # A finite B resolves Kramers partners. No derivative is assigned across
    # unresolved degeneracies; such points require a projector treatment.
    if derivative and spacing < 1e-10:
        raise ValueError(f'Unresolved same-sector degeneracy {spacing}')
    wbdiag=np.diagonal(wb,axis1=1,axis2=2).real
    gap_slope=wbdiag[:,:4,None]-wbdiag[:,None,4:]
    result=[]
    if derivative:
        de=e[:,None,:]-e[:,:,None] # e_n-e_m for K_mn
        inv=np.divide(1.,de,out=np.zeros_like(de),where=abs(de)>1e-10)
        K=wb*inv
        sectors={}
        for name,sl1,sl2 in [('occ_rotation',slice(0,4),slice(0,4)),('empty_rotation',slice(4,8),slice(4,8))]:
            q=np.zeros_like(K);q[:,sl1,sl2]=K[:,sl1,sl2];sectors[name]=q
        sectors['cross_gap_rotation']=K-sectors['occ_rotation']-sectors['empty_rotation']
    for spin in spins:
        sb=adj(u)@spin@u
        n=-2*np.imag(sb[:,:4,4:]*vb[:,4:,:4].swapaxes(1,2))
        vals={'chi':(n/gaps**2).sum(axis=(1,2))}
        if derivative:
            vals['energy_shift']=(-2*n*gap_slope/gaps**3).sum(axis=(1,2))
            for name,q in sectors.items():
                ds=sb@q-q@sb;dv=vb@q-q@vb
                dn=-2*np.imag(ds[:,:4,4:]*vb[:,4:,:4].swapaxes(1,2)+sb[:,:4,4:]*dv[:,4:,:4].swapaxes(1,2))
                vals[name]=(dn/gaps**2).sum(axis=(1,2))
            vals['matrix_change']=sum(vals[x] for x in sectors)
            vals['tangent']=vals['energy_shift']+vals['matrix_change']
        result.append(vals)
    return result, e[:,3].max(),e[:,4].min(),spacing

def run_case(lam,t,grid,b,u_mev=0.,derivative=True):
    model=AuditBilayer(AuditSpec(shift=ab_ba_shift(lam),staggered_eV=.045,tunneling_eV=t/1000))
    direction=np.array([.5,np.sqrt(3)/2])
    spins=[sum(direction[j]*model.spins[name+'_'+ax] for j,ax in enumerate(('x','y'))) for name in ('bottom','top','total')]
    k,w,_=integration_grid(model,*grid)
    totals={};vmax=-np.inf;cmin=np.inf;spacing=np.inf
    for start in range(0,len(k),512):
        ks=k[start:start+512];weights=w[start:start+512]
        h=model.matrices(ks)+(u_mev/1000)*W
        v=sum(direction[j]*model.matrices(ks,axis=j) for j in range(2))
        for sign in (1,-1):
            hs=h-sign*b*model.mono.g_factor*MU_B*model.spins['total_z']
            rows,vm,cm,sp=kernel(hs,v,spins,derivative)
            vmax=max(vmax,vm);cmin=min(cmin,cm);spacing=min(spacing,sp)
            for name in rows[0]:
                val=np.array([weights@r[name] for r in rows])*sign/(2*b)
                totals[name]=totals.get(name,np.zeros(3))+val
    assert cmin-vmax>0,'Insulating filling lost'
    return dict(fraction=lam,t_meV=t,U_meV=u_mev,grid=grid,nk=len(k),Bstep_T=b,
                order=['lower','upper','net'],C_units='(nm V T)^-1',G_units='(nm V T eV)^-1',
                common_gap_meV=(cmin-vmax)*1000,minimum_same_sector_spacing_eV=spacing,
                **{name:val.tolist() for name,val in totals.items()})

def audit():
    m=AuditBilayer(AuditSpec(shift=ab_ba_shift(.25),staggered_eV=.045,tunneling_eV=.01))
    k=np.array([[.13,.29],[8.1,-3.2],[-5.8,12.2]])
    h=m.matrices(k)-.0625*m.mono.g_factor*MU_B*m.spins['total_z']
    v=.5*m.matrices(k,axis=0)+np.sqrt(3)/2*m.matrices(k,axis=1)
    sp=[.5*m.spins['bottom_x']+np.sqrt(3)/2*m.spins['bottom_y']]
    exact=kernel(h,v,sp)[0][0]
    errors=[]
    for du in (1e-5,5e-6):
        plus=kernel(h+du*W,v,sp,False)[0][0]['chi'];minus=kernel(h-du*W,v,sp,False)[0][0]['chi']
        errors.append(float(np.max(abs((plus-minus)/(2*du)-exact['tangent']))))
    phased=kernel(h,v,sp,phase_seed=37)[0][0]
    phase_error=max(float(np.max(abs(phased[x]-exact[x]))) for x in exact)
    assert max(errors)<1e-7 and phase_error<1e-9,(errors,phase_error)
    return {'offgrid_tangent_errors':errors,'random_phase_error':phase_error}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--grid',nargs=3,type=int,default=[21,24,36]);p.add_argument('--b',type=float,default=.0625)
    p.add_argument('--t',nargs='+',type=float,default=[10,15,20]);p.add_argument('--fraction',type=float,default=.25)
    p.add_argument('--mode',choices=['predict','verify'],default='predict');p.add_argument('--output',required=True)
    args=p.parse_args();out=ROOT/args.output;start=time.monotonic();a=audit();rows=[]
    biases=[0.] if args.mode=='predict' else [-.5,-.25,.25,.5]
    for t in args.t:
        for U in biases:
            row=run_case(args.fraction,t,args.grid,args.b,U,args.mode=='predict');rows.append(row)
            print(json.dumps(row),flush=True)
            out.write_text(json.dumps(dict(mode=args.mode,audit=a,rows=rows,complete=False),indent=2)+'\n')
    out.write_text(json.dumps(dict(mode=args.mode,audit=a,rows=rows,complete=True,wall_seconds=time.monotonic()-start),indent=2)+'\n')
