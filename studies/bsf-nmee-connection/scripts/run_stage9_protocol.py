"""B-prepared E response, independent mixed-Liouville limit and layer response."""
import json
from dataclasses import asdict,replace
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from bsf_nmee.bilayer_audit import AuditSpec,AuditBilayer,integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift,half_filled_spec,prepared_B_response
from bsf_nmee.static_observables import observables
from nmee_lab.ptbackend.liouville_general import eb_single


def main():
    k,w,meta=integration_grid(AuditBilayer(),41,45,72)
    out={'grid':meta,'scope':'gapped fixed-orbital-gauge; B-prepared low-frequency E and comparison ramp protocols'}
    def save():
        (ROOT/'results/bilayer_stage9_protocol.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    out['prepared']={}
    for name,shift in [('AA',(0.,0.)),('AB',ab_ba_shift(0)),('BA',ab_ba_shift(1)),
                       ('generic_small',(.023,.041)),('generic_large',(.092,.164))]:
        spec=half_filled_spec(AuditSpec(shift=shift,tunneling_eV=.005),k);model=AuditBilayer(spec)
        row={'parameters':asdict(spec),'static':observables(model,k,w,step_T=.0625)}
        row['frequency']=prepared_B_response(model,k,w,[0.,.0001,.001,.003],eta_eV=1e-6)
        row['top']=prepared_B_response(model,k,w,[0.],eta_eV=1e-6,observation='top')
        row['bottom']=prepared_B_response(model,k,w,[0.],eta_eV=1e-6,observation='bottom')
        out['prepared'][name]=row
        print('prepared',name,flush=True);save()
    out['mixed_limits']=[]
    for name,shift in [('AB',ab_ba_shift(0)),('generic_small',(.023,.041))]:
        spec=half_filled_spec(AuditSpec(shift=shift,tunneling_eV=.005),k);model=AuditBilayer(spec)
        for eta_e,eta_b in [(1e-4,1e-6),(1e-5,1e-7),(1e-5,1e-5),(1e-6,1e-6),(1e-7,1e-5)]:
            value=0j
            for start in range(0,len(k),512):
                ks=k[start:start+512]
                v=eb_single(model.matrices(ks),model.matrices(ks,axis=0),
                   model.spins['total_z'],model.spins['total_x'],spec.mu_eV,0.,
                   model.mono.g_factor,1j*eta_e,1j*eta_b)
                value+=w[start:start+512]@v
            out['mixed_limits'].append({'model':name,'eta_E_eV':eta_e,'eta_B_eV':eta_b,
                                         'Cxx_real':float(value.real),'Cxx_imag':float(value.imag)})
            print('mixed',name,eta_e,eta_b,value,flush=True)
        save()


if __name__=='__main__':main()
