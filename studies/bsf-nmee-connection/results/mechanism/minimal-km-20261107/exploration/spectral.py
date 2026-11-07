exec(open(__file__.replace('spectral.py','explore.py')).read().split('rows=[]')[0])
sys.path.insert(0,str(repo/'studies/bsf-nmee-connection/scripts'))
from compare_reversal_mechanism import transplant
from bsf_nmee.remote_empty_km import selective_extra
angle=np.array([.5,np.sqrt(3)/2]); rows={};errors=[]
def spectra(points,g):
 data=ingredients(points)
 h,v=selective_extra(data,.12,.01,g,.1)
 e,u=np.linalg.eigh(h)
 return e,u,h,v,data
for a in [.04,.12]:
 for b in [.04,.12]:
  val=np.zeros(3);vmax=-np.inf;cmin=np.inf
  for start in range(0,len(k),512):
   points=k[start:start+512];wt=w[start:start+512]
   ea,ua,ha,va,data=spectra(points,a);eb,ub,hb,vb,_=spectra(points,b)
   hc=transplant(ea,ua,eb)
   if a==b:vc=va
   else:
    derivatives=[]
    for dk in [1e-5,5e-6]:
     plus=spectra(points+dk*angle,a);ep=spectra(points+dk*angle,b)[0]
     minus=spectra(points-dk*angle,a);em=spectra(points-dk*angle,b)[0]
     derivatives.append((transplant(plus[0],plus[1],ep)-transplant(minus[0],minus[1],em))/(2*dk))
    errors.append(float(np.max(abs(derivatives[0]-derivatives[1]))));vc=derivatives[-1]
   c,ev,ec=response(hc,vc,data[-3],data[-2],wt,g=data[-1]);val+=c
   vmax=max(vmax,ev);cmin=min(cmin,ec)
  rows[f'states{a}_energies{b}']={'C':val.tolist(),'gap':cmin-vmax}
  print(a,b,val*1e5,cmin-vmax,flush=True)
assert max(errors)<1e-5
Path(__file__).with_name('spectral.json').write_text(json.dumps({'rows':rows,'derivative_error':max(errors),'grid':[21,20,32]},indent=2)+'\n')
