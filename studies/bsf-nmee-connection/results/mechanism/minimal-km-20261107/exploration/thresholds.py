exec(open(__file__.replace('thresholds.py','explore.py')).read().split('rows=[]')[0])
from scipy.optimize import brentq
data=ingredients(k,site=0)
rows=[];roots=[]
def calc(energy,t,rho=0.,site=0,tlow=0.):
 dd=data if site==0 else ingredients(k,site=1)
 channel={'VC':1.,'CC':rho,'VV':1.}
 h,v=assemble(dd,energy,tlow,t,channel)
 c,ev,ec=response(h,v,dd[-3],dd[-2],w,g=dd[-1])
 r=dict(energy=energy,t=t,rho=rho,site=site,tlow=tlow,C=c.tolist(),gap=ec-ev);rows.append(r)
 return c[0]
for energy in [.12,.25,.5]:
 root=brentq(lambda t:calc(energy,t),.06,.25,xtol=1e-5)
 roots.append(dict(energy=energy,t=root,scale=root**2/energy));print('ROOT',roots[-1],flush=True)
# Prediction made only from these training energies; retain before held-out checks.
scale=float(np.median([r['scale'] for r in roots]))
predictions=[dict(energy=e,predicted_t=np.sqrt(scale*e),criterion_g2_over_Delta=scale) for e in [.18,.35]]
Path(__file__).with_name('predictions-before-check.json').write_text(json.dumps({'training':roots,'predictions':predictions},indent=2)+'\n')
for pred in predictions:
 e=pred['energy'];t=pred['predicted_t']
 pred['below']=calc(e,.9*t);pred['above']=calc(e,1.1*t)
 pred['actual_root']=brentq(lambda x:calc(e,x),.06,.25,xtol=1e-5)
 print('PREDICTION',pred,flush=True)
for rho in [0.,.1,.25,.5,1.]:
 for t in [.1,.12,.15]:print('ROBUST',rho,t,calc(.12,t,rho),flush=True)
for tlow in [.005,.01]:
 print('LOW',tlow,calc(.12,.12,.1,tlow=tlow),flush=True)
Path(__file__).with_name('thresholds.json').write_text(json.dumps({'rows':rows,'training':roots,'prediction_checks':predictions},indent=2)+'\n')
