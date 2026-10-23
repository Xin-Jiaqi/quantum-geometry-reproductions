"""Reproduce v74 source-specific selected electric pole channel / path audit."""
from pathlib import Path
import argparse,json
from dataclasses import replace
from qgmap.selected_electric_pole_v74 import (
    Spec,channels,paths,finite_frequency_negative_control,literal_eq20,printed_eq31)

def as_pair(z):return {'real':float(z.real),'imag':float(z.imag)}

def convert(obj):
    if isinstance(obj,complex):return as_pair(obj)
    if isinstance(obj,dict):return {k:convert(v) for k,v in obj.items()}
    if isinstance(obj,list):return [convert(v) for v in obj]
    if isinstance(obj,tuple):return [convert(v) for v in obj]
    return obj

def generate():
    s=Spec(nr=64,nphi=88,eta_e=.003,eta_b=.001)
    o=channels(s)
    direct=literal_eq20(s)
    print31=printed_eq31(s)
    return convert({
       'version':'v74',
       'scientific_status':'PARTIAL',
       'model':'Homogeneous spin-Zeeman two-colour Dirac Liouville, independent Fermi Frechet derivative, source Eq20 literal contraction. No relaxation bath or transport dc claim.',
       'spec':{k:getattr(s,k) for k in s.__dataclass_fields__},
       'channels':o,
       'source_literal_eq20_physical_B':direct,
       'source_printed_eq31_physical_B':print31,
       'source_sign_discrepancy':direct-print31,
       'ratio_scaled_full_to_eq20':o['eta_e_scaled_full']/direct,
       'regulator_paths':paths(s),
       'finite_omega_control':finite_frequency_negative_control(s),
       'limits':{
          'electric_selected_pole_only':'E diagonal population contribution isolated algebraically; E-offdiagonal + B-first also have output-pole effects and must not be silently omitted',
          'joint_regulators':'Fixed finite eta_B path is not the joint eta_E/eta_B -> 0 source limit',
          'no_physical_bath':'adiabatic infinitesimals are not a uniform scattering rate; Eq15 phenomenological tau substitution is not derived here',
          'source_conflict':'SM Eq20 and printed Eq31 remain opposite with indexed v=i(e_m-e_n) r; source author intent unknown'
       }})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='results/frozen/v74_selected_electric_pole.json')
    args=ap.parse_args();dst=Path(args.output);dst.parent.mkdir(exist_ok=True,parents=True)
    dst.write_text(json.dumps(generate(),sort_keys=True,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('v74 frozen science data:',dst)
