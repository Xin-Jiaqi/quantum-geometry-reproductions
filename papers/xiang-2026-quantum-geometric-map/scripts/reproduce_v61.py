"""Reproduce v61 same-QWZ four-terminal charge edge versus periodic bulk Chern."""
import argparse,json
from pathlib import Path
from qgmap.qwz_device_v61 import qwz_four_terminal

def create():
    cases={}
    for m in (-1.,1.,3.):
        for L in (3,5,7):
            d=qwz_four_terminal(L,m,eta=1e-4)
            d.pop('transmissions')
            cases[f'm{m:+g}_L{L}']=d
    eta_scan={}
    for eta in (1e-3,1e-4,1e-5):
        d=qwz_four_terminal(5,-1,eta=eta)
        eta_scan[str(eta)]={'hall_resistance_h_over_q2':d['hall_resistance_h_over_q2'],
                            'drive_current_q2_over_h':d['drive_current_q2_over_h']}
    return {'version':'v61','model':'same-QWZ-Fourier-Hamiltonian-in-bulk-device-and-four-semi-infinite-leads',
            'units':'q2_over_h and h_over_q2',
            'status':'CONDITIONAL coherent zero-T clean charge edge; NOT intrinsic orbital mixed EB',
            'zero_current_guard':'Rxy=null when |Ileft|<1e-7 for normalized unit voltage bias',
            'cases':cases,'eta_scan_mminus1_L5':eta_scan}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='results/frozen/v61_same_qwz_four_terminal.json')
    args=parser.parse_args()
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(create(),sort_keys=True,indent=2)+'\n')
    print('v61 results saved:',out)
