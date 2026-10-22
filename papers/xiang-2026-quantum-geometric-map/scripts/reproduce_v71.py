"""Reproduce v71 source-observable separation in three frames, independent of frozen JSON."""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.source_observable_bridge_v71 import (
    Spec,finite_contact_identity,source_eq20_eq31_sign,pauli_helicity_gate)

KEYS=('dynamic_device','dynamic_contact','dynamic_full','static_isothermal_device',
      'static_isothermal_contact','static_isothermal_full','dynamic_ward_residual')

def value(x):
    if isinstance(x,complex):return {'real':float(x.real),'imag':float(x.imag)}
    return float(x)

def generate():
    spec=Spec(length=22,eta=.025,omega=.29)
    frames={k:{key:value(v[key]) for key in KEYS}
            for k in ('symmetric','landau_y','orbital')
            for v in [finite_contact_identity(spec,k)]}
    rng=np.random.default_rng(734)
    r=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));r=(r+r.conj().T)/2
    s=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));s=(s+s.conj().T)/2
    z=source_eq20_eq31_sign(np.array([-.8,-.21,.51,1.2]),r,s,
                            np.array([.3,-.11,.19,-.25]))
    return {'version':'v71','frame_current_responses':frames,
        'literal_sm_sign':{key:value(z[key]) for key in ('direct_eq20','inferred_Z','printed_eq31','algebra_residual','sign_discrepancy')},
        'source_helicity':pauli_helicity_gate(),
        'model_status':'VERIFIED finite full/device/contact and literal source-index algebra; Xiang bulk article PARTIAL',
        'source_sm_sha256':'9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4'}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='results/frozen/v71_source_observable_bridge.json')
    args=parser.parse_args();p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(generate(),indent=2,sort_keys=True)+'\n')
    print('v71 source-observable evidence',p)
