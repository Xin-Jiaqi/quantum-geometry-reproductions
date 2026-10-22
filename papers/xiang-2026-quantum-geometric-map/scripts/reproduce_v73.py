"""Independent numerical frozen reproduction v73, two-color source + static magnetic prep."""
from pathlib import Path
import argparse,json
from dataclasses import replace
from qgmap.two_color_finite_q_spin_v73 import (Spec,two_color,two_color_dense_independent,
    static_magnetic_derivative,regulator_order_scan)

def pair(x):return {'real':float(x.real),'imag':float(x.imag)}
def generate():
    s=Spec(n=24)
    main=two_color(s)
    dense=two_color_dense_independent(s)
    gauge=two_color(replace(s,orbital_gauge=.71))
    static={}
    for eta in (.001,.0001,.00001):
        a=static_magnetic_derivative(replace(s,omega_b=0,eta_b=eta))
        static[str(eta)]={k:(pair(v) if isinstance(v,complex) else v)
                           for k,v in a.items()}
    scan=regulator_order_scan(s)
    scanconv={k:{kk:pair(vv) for kk,vv in item.items()}
              for k,item in scan['values'].items()}
    return dict(version='v73',source_status='PARTIAL',
        physical_scope='periodic finite-q scalar-potential/spin-Zeeman two-color; NOT original adiabatic selected tau-pole',
        parameters={key:getattr(s,key) for key in s.__dataclass_fields__},
        two_color={key:{k:pair(v) for k,v in item.items()}
                   for key,item in [(name,main[name]) for name in ('b_then_e','e_then_b','total')]},
        frequency_sum=pair(main['frequency_sum']),ward_residual=pair(main['ward_residual']),
        dense_total_E_transverse=pair(dense['total']['normalized_E_transverse']),
        gauged_total_E_transverse=pair(gauge['total']['normalized_E_transverse']),
        static_B_preparation=static,
        finite_q_zero_frequency_regulator_scan=scanconv)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='results/frozen/v73_two_color_finite_q_spin.json')
    args=parser.parse_args();p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(generate(),indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('v73 frozen evidence written:',p)
