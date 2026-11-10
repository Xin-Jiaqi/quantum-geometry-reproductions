"""Fixed separated-layer VC interventions in unchanged eight-state KM."""
from pathlib import Path
import os,sys,json,argparse,time
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
ROOT=Path(__file__).resolve().parent
candidates=list(ROOT.parents)+[p/'quantum-geometry-reproductions' for p in ROOT.parents]
REPO=next(p for p in candidates if (p/'studies/bsf-nmee-connection/bsf_nmee').is_dir())
sys.path[:0]=[str(REPO/'studies/bsf-nmee-connection'),str(REPO/'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from bsf_nmee.bilayer_audit import AuditBilayer,AuditSpec,integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift
from bsf_nmee.sliding_tight_binding import P_TOP,P_BOTTOM
from nmee_lab.clean_core.core import MU_B
W=(P_TOP-P_BOTTOM)/2;adj=lambda a:a.conj().swapaxes(-1,-2)
DIRECTION=np.array([.5,np.sqrt(3)/2])
def channels(model,k):
 h=model.matrices(k);v=sum(DIRECTION[j]*model.matrices(k,axis=j) for j in range(2))
 p=np.zeros_like(h);dp=np.zeros_like(h);intra=np.zeros_like(h);vi=np.zeros_like(h)
 for sl in (slice(0,4),slice(4,8)):
  hl=h[:,sl,sl];vl=v[:,sl,sl];e,u=np.linalg.eigh(hl);occ=u[:,:,:2]
  assert np.min(e[:,2]-e[:,1])>1e-8
  p[:,sl,sl]=occ@adj(occ);intra[:,sl,sl]=hl;vi[:,sl,sl]=vl
  f=np.array([1,1,0,0]);df=f[:,None]-f[None,:];de=e[:,:,None]-e[:,None,:]
  rat=np.divide(df,de,out=np.zeros_like(de),where=df!=0)
  dp[:,sl,sl]=u@(rat*(adj(u)@vl@u))@adj(u)
 t=h-intra;dt=v-vi;q=np.eye(8)-p;pieces={};dpieces={}
 for name,lower,upper in [('Vb_Ct',P_BOTTOM,P_TOP),('Vt_Cb',P_TOP,P_BOTTOM)]:
  a=p@lower;b=q@upper;da=dp@lower;db=-dp@upper
  z=a@t@b;dz=da@t@b+a@dt@b+a@t@db
  pieces[name]=z+adj(z);dpieces[name]=dz+adj(dz)
 return h,v,pieces,dpieces

def chi(h,v,spins):
 e,u=np.linalg.eigh(h);occ=u[:,:,:4];emp=u[:,:,4:];ve=adj(emp)@v@occ
 de=e[:,:4,None]-e[:,None,4:]
 val=np.array([(-2*np.imag((adj(occ)@s@emp)*ve.swapaxes(1,2))/de**2).sum((1,2)) for s in spins])
 return val,float(e[:,3].max()),float(e[:,4].min())

def audit():
 m=AuditBilayer(AuditSpec(shift=ab_ba_shift(.25),staggered_eV=.045,tunneling_eV=.01))
 k=np.array([[.13,.29],[8.1,-3.2],[-5.8,12.2]])
 h,v,cs,ds=channels(m,k);err=[]
 for dk in (1e-5,5e-6):
  cp=channels(m,k+dk*DIRECTION)[2];cm=channels(m,k-dk*DIRECTION)[2]
  err.append(max(float(abs((cp[n]-cm[n])/(2*dk)-ds[n]).max()) for n in cs))
 assert max(err)<1e-7,err
 return {'channel_velocity_finite_difference_errors':err,'hermiticity_error':max(float(abs(c-adj(c)).max()) for c in cs.values())}

def run(args):
 out={'complete':False,'scope':'Unchanged 8-state KM, fixed separated-layer directional VC projector intervention. U is internal offset, not physical gate field.',
      'order':['lower','upper','net'],'C_units':'(nm V T)^-1','G_units':'(nm V T eV)^-1','audit':audit(),'rows':[],'parameters':vars(args).copy()}
 out['parameters']['output']=str(args.output)
 starttime=time.monotonic()
 for t in args.t:
  m=AuditBilayer(AuditSpec(shift=ab_ba_shift(args.fraction),staggered_eV=args.mass/1000,tunneling_eV=t/1000))
  spins=[sum(DIRECTION[j]*m.spins[name+'_'+ax] for j,ax in enumerate(('x','y'))) for name in ('bottom','top','total')]
  k,w,_=integration_grid(m,*args.grid)
  ctrls=[('full',1.)]+[(c,a) for c in args.channels for a in args.alpha]
  accum={(c,a,U,d):{'C':np.zeros(3),'vm':-np.inf,'cm':np.inf} for c,a in ctrls for U in args.U for d in (-1,0,1)}
  for start in range(0,len(k),512):
   ks=k[start:start+512];ww=w[start:start+512];h,v,cs,ds=channels(m,ks)
   for c,a in ctrls:
    hc=h if c=='full' else h-(1-a)*cs[c];vc=v if c=='full' else v-(1-a)*ds[c]
    for U in args.U:
     for d in (-1,0,1):
      hh=hc+(U+d*args.du)/1000*W;r=accum[c,a,U,d]
      for sign in (1,-1):
       z,vm,cm=chi(hh-sign*args.b*m.mono.g_factor*MU_B*m.spins['total_z'],vc,spins)
       r['C']+=sign*(z@ww)/(2*args.b);r['vm']=max(r['vm'],vm);r['cm']=min(r['cm'],cm)
  for c,a in ctrls:
   for U in args.U:
    rs=[accum[c,a,U,d] for d in (-1,0,1)]
    cm,c0,cp=[r['C'] for r in rs];g=(cp-cm)/(2*args.du/1000)
    gap=min(r['cm'] for r in rs)-max(r['vm'] for r in rs)
    out['rows'].append({'t_meV':t,'U_meV':U,'channel':c,'retained_alpha':a,'C0':c0.tolist(),'Cminus':cm.tolist(),'Cplus':cp.tolist(),'G':g.tolist(),'common_gap_eV':gap,'insulating':bool(gap>0),'nk':len(k)})
  args.output.write_text(json.dumps(out,indent=2)+'\n');print('DONE t',t,'nk',len(k),'seconds',time.monotonic()-starttime,flush=True)
 out['complete']=True;out['seconds']=time.monotonic()-starttime;args.output.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--grid',nargs=3,type=int,default=[21,24,36]);p.add_argument('--t',nargs='+',type=float,default=[10,20]);p.add_argument('--U',nargs='+',type=float,default=[-4,0,4]);p.add_argument('--mass',type=float,default=45);p.add_argument('--fraction',type=float,default=.25);p.add_argument('--du',type=float,default=.1);p.add_argument('--b',type=float,default=.0625);p.add_argument('--channels',nargs='+',default=['Vb_Ct','Vt_Cb']);p.add_argument('--alpha',nargs='+',type=float,default=[0]);p.add_argument('--output',type=Path,required=True);run(p.parse_args())
