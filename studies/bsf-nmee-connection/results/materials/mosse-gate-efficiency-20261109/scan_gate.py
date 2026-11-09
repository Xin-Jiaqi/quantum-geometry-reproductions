"""Gate sensitivity of net NMEE, with matched controls and fixed zero-gate channels."""
from pathlib import Path
import os,socket,json
import numpy as np
from analyze import Response,embedding,adj
from mixing_channels import split_channels
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute node required')
root=Path.cwd();source=Path('/work/sunq/xjq/mosse-strain-20261108/e000')
d=np.load(source/'A/cache.npz');dm=np.load(source/'mono/cache.npz');ops=np.load(root/'operators.npz')
fields=np.array([-.5,-.25,-.1,-.05,-.025,0,.025,.05,.1,.25,.5])
result={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'field_grid':fields.tolist(),
        'units':'C: 1e-9 (nm V T)^-1, Pauli spin density; F_model: imposed potential slope V/nm',
        'scope':__doc__,'field_is_screened_external':False,'references':[],'bilayer':[],'diagnostics':[]}
def save(): (root/'scan.json').write_text(json.dumps(result,indent=2))
def response(r,h):
    try:
        ans=r.response(h)
        if ans['gap_eV']<.05:return {'invalid':'gap below 0.05 eV','gap_eV':ans['gap_eV']}
        return ans
    except ValueError as e:return {'invalid':str(e)}

# Real-space field acts within isolated layers too; upper is the mirrored lower.
for nb in (48,56,64):
    r=Response(dm,nb);v=ops['mono_gate'][:,:nb,:nb]
    one={float(f):response(r,r.h+f*v) for f in fields}
    gauge0=r.response();gauge1=r.response(r.h+.2*np.eye(nb))
    result['diagnostics'].append({'mono_nb':nb,'rigid_shift_invariance':float(np.max(abs(np.array(gauge0['C_total_lower_upper'])-gauge1['C_total_lower_upper'])))})
    for f in fields:
        lo=one[float(f)];up=one[float(-f)]
        if 'invalid' not in lo and 'invalid' not in up:
            cl=lo['C_total_lower_upper'][0];cu=-up['C_total_lower_upper'][0]
            result['references'].append({'kind':'isolated_pair','nb':nb,'F_model':float(f),
                'C':[cl+cu,cl,cu],'gap_eV':min(lo['gap_eV'],up['gap_eV'])})
        else:result['references'].append({'kind':'isolated_pair','nb':nb,'F_model':float(f),'invalid':[lo,up]})
        for i,end in enumerate((3.2,3.8)):
            lo=response(r,r.h+ops['lower_reference'][i,:,:nb,:nb]+f*v)
            up=response(r,r.h+ops['upper_reference'][i,:,:nb,:nb]-f*v)
            row={'kind':'electrostatic_pair','nb':nb,'taper_end_A':end,'F_model':float(f)}
            if 'invalid' not in lo and 'invalid' not in up:
                cl=lo['C_total_lower_upper'][0];cu=-up['C_total_lower_upper'][0]
                row.update(C=[cl+cu,cl,cu],gap_eV=min(lo['gap_eV'],up['gap_eV']))
            else:row['invalid']=[lo,up]
            result['references'].append(row)
    save();print('REFERENCES',nb,flush=True)

bottom=d['orb_bottom'];layers=[np.flatnonzero(bottom),np.flatnonzero(~bottom)]
mask_b=np.diag(bottom.astype(float));mask_t=np.eye(len(bottom))-mask_b
off=bottom[:,None]!=bottom[None,:]
for nb in (96,112):
    r=Response(d,nb);v=ops['A_gate'][:,:nb,:nb]
    for frozen in (4,8,12):
        lift,s=embedding(d,frozen,nb)
        hloc=adj(lift)@r.h@lift;vloc=adj(lift)@v@lift
        _,_,p,checks=split_channels(hloc,layers,14);q=np.eye(len(bottom))-p
        pb=p@mask_b;pt=p@mask_t;qb=q@mask_b;qt=q@mask_t
        def parts(x):
            t=x*off;vv=p@t@p;cc=q@t@q
            bc=pb@t@qt;bc=bc+adj(bc)
            tc=pt@t@qb;tc=tc+adj(tc)
            assert abs(vv+cc+bc+tc-t).max()<1e-8
            return {'VV_CC':vv+cc,'Vb_Ct':bc,'Vt_Cb':tc,'VC':bc+tc,'all':t}
        hparts=parts(hloc);vparts=parts(vloc)
        controls={'full':(), 'no_VC':('VC',), 'VC_only':('VV_CC',), 'no_selected_hopping':('all',)}
        if frozen==8:controls.update(no_Vb_Ct=('Vb_Ct',),no_Vt_Cb=('Vt_Cb',))
        for f in fields:
            h=r.h+f*v
            for name,removed in controls.items():
                delta=sum((hparts[k]+f*vparts[k] for k in removed),np.zeros_like(hloc))
                hc=h-lift@delta@adj(lift);ans=response(r,hc)
                row={'nb':nb,'frozen_empty':frozen,'F_model':float(f),'control':name,**ans}
                if name=='full' and frozen==8 and 'invalid' not in ans:
                    en,u=np.linalg.eigh(h);occ=u[:,:,:r.nocc]
                    rho=adj(lift)@occ@adj(occ)@lift
                    row['empty_layer_weights']=[float(np.trace(rho@qq,axis1=-2,axis2=-1).real.mean()) for qq in (qb,qt)]
                    row['VC_coherence_energies_eV']={k:float(np.trace(rho@(hparts[k]+f*vparts[k]),axis1=-2,axis2=-1).real.mean()) for k in ('Vb_Ct','Vt_Cb')}
                    el=[np.linalg.eigvalsh((hloc+f*vloc)[:,ids[:,None],ids]) for ids in layers]
                    row['decoupled_cross_gaps_eV']=[float((el[1][:,14]-el[0][:,13]).min()),float((el[0][:,14]-el[1][:,13]).min())]
                    half=r.response(h,b=.001)
                    row['half_probe_step_max_change']=float(abs(np.array(ans['C_total_lower_upper'])-half['C_total_lower_upper']).max())
                result['bilayer'].append(row)
            save()
        result['diagnostics'].append({'nb':nb,'frozen_empty':frozen,'projection_min_singular':s,**checks})
        print('CHANNELS',nb,frozen,flush=True)
result['complete']=True;save();print('GATE_SCAN_COMPLETE',flush=True)
