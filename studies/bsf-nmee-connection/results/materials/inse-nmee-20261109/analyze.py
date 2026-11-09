"""Frozen-potential response: full KS space plus a selected d+p mixing intervention."""
from pathlib import Path
import json,os,socket
import numpy as np
from scipy.linalg import null_space
adj=lambda a:a.conj().swapaxes(-1,-2)
def polar(a):
    u,s,vh=np.linalg.svd(a,full_matrices=False)
    return u@vh,s
class Response:
    def __init__(self,d,nb=None):
        self.nk=len(d['k']);self.mesh=round(self.nk**.5);self.cell=d['cell'];self.nocc=int(d['nocc'])
        self.nb=nb or d['energies'].shape[1];self.links=d['links'][:,:,:self.nb,:self.nb];self.neighbors=d['neighbors']
        self.z=d['Z'][:,:self.nb,:self.nb];self.ops=np.array([d['spin_total'][:,:self.nb,:self.nb],d['spin'][:,:self.nb,:self.nb]])
        self.h=d['energies'][:,:self.nb,None]*np.eye(self.nb)
        self.factor=1e10*(2.00231930436256*5.7883818060e-5/2)/np.linalg.norm(np.cross(self.cell[0],self.cell[1]))
    def chi(self,h):
        en,u=np.linalg.eigh(h);occ=u[:,:,:self.nocc];emp=u[:,:,self.nocc:];deriv=np.zeros_like(occ,dtype=complex)
        smallest=1.
        for a in range(4):
            v=self.links[:,a]@occ[self.neighbors[:,a]];o=adj(occ)@v
            smallest=min(smallest,float(np.linalg.svd(o,compute_uv=False).min()))
            deriv+=(-1)**a*self.mesh*self.cell[a//2,0]/(4*np.pi)*v@np.linalg.inv(o)
        gap=float(en[:,self.nocc].min()-en[:,self.nocc-1].max())
        if gap<=0:raise ValueError('Gap closed; insulating response invalid')
        if smallest<1e-5:raise ValueError('Occupied neighboring spaces singular')
        dip=adj(emp)@(1j*deriv);delta=en[:,:self.nocc,None]-en[:,None,self.nocc:]
        cs=[float((2*((adj(occ)@s@emp)*dip.swapaxes(1,2)).real/delta).sum(axis=(1,2)).mean()*self.factor) for s in self.ops]
        return np.array([cs[0],cs[1],cs[0]-cs[1]]),gap,smallest
    def response(self,h=None,b=.002):
        h=self.h if h is None else h
        cp,gp,sp=self.chi(h+b*self.z);cm,gm,sm=self.chi(h-b*self.z)
        c=(cp-cm)/(2*b)
        return {'C_total_lower_upper':c.tolist(),'even_chi':((cp+cm)/2).tolist(),
                'gap_eV':min(gp,gm),'min_neighbor_singular':min(sp,sm),'b_eV':b}
def embedding(d,frozen_empty,nb):
    nc=int(d['ncore']);nv=int(d['nocc'])-nc;nf=nv+frozen_empty;ls=[];mins=[]
    for coeff,energy in zip(d['coeff'],d['energies'][:,:nb]):
        cf=coeff[:,:nb-nc].astype(complex);f,s=polar(cf[:,:nf]);rest=null_space(adj(f));emin=energy[nc+nf]
        tail,_=polar((adj(rest)@cf[:,nf:]/(1+((energy[nc+nf:]-emin)/4)**2)).conj().T)
        lift=np.vstack([np.zeros((nc,len(cf))),adj(f),tail@adj(rest)])
        if abs(adj(lift)@lift-np.eye(len(cf))).max()>1e-8:raise ValueError('Nonisometric lift')
        ls.append(lift);mins.append(float(s.min()))
    return np.array(ls),min(mins)
def run(root):
    if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute required')
    out={'units':'1e-9 (nm V T)^-1, Pauli spin density','scope':'Frozen-potential static spin-only EB derivative; no self-consistent screening',
         'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'materials':{}}
    for label in ('mono','A'):
        d=np.load(root/label/'cache.npz');r=Response(d);zero,gap,sv=r.chi(r.h)
        rec={'zero_chi':zero.tolist(),'gap_eV':gap,'baseline':[r.response(b=b) for b in (.001,.002,.008)],'cutoffs':{},'mixing':{}}
        out['materials'][label]=rec
        (root/'analysis.json').write_text(json.dumps(out,indent=2))
        for nb in ((40,48,56,64) if label=='mono' else (64,80,96,112)):
            rr=Response(d,nb);rec['cutoffs'][str(nb)]=rr.response()
        if label=='A':
            for nb in (96,112):
                rr=Response(d,nb)
                for frozen in (4,8,12):
                    lift,s=embedding(d,frozen,nb);hlocal=adj(lift)@rr.h@lift
                    bottom=d['orb_bottom'];off=bottom[:,None]!=bottom[None,:]
                    ht=lift@(hlocal*off)@adj(lift)
                    rows=[]
                    for alpha in (0,.5,1):
                        rows.append(dict(alpha=alpha,**rr.response(rr.h+(alpha-1)*ht)))
                    rec['mixing'][f'{nb}_{frozen}']={'min_projection_singular':s,'rows':rows,
                        'interpretation':'Scale selected d+p interlayer H only. Residual cross-layer vertices/states retained; alpha=0 is not pure electrostatics.'}
                    (root/'analysis.json').write_text(json.dumps(out,indent=2))
        (root/'analysis.json').write_text(json.dumps(out,indent=2));print(label,rec['baseline'],flush=True)
    print('ANALYSIS_DONE',flush=True)
if __name__=='__main__':run(Path.cwd())
