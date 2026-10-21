"""Deterministic finite-size / probe / eta ordered-limit evidence for v64."""
import argparse,json
from pathlib import Path
from qgmap.qwz_virtual_probe_v64 import virtual_probe_limit
from qgmap.qwz_probe_limit_v63 import probe_strength_response


def compact(v):
    return {'Rxy_lambda_to_0_h_over_q2':v['Rxy_lambda_to_0_h_over_q2'],
            'two_terminal_T_lr':v['two_terminal_T_lr'],
            'two_terminal_T_rl':v['two_terminal_T_rl'],
            'probe_spectral_weights':v['probe_spectral_weights'].tolist(),
            'spectral_identity_error':v['spectral_identity_error'],
            'limit_defined':v['limit_defined']}


def create():
    result={'version':'v64','status':'VERIFIED algebraic virtual-probe limit; CONDITIONAL finite-device; original orbital E/B UNRESOLVED',
            'scope':'finite L, fixed positive eta and mu; lambda->0+ before L/eta variation',
            'size':{},'position':{},'eta':{},'lambda':{},'metal':{}}
    for L in (5,7,9,11):
        result['size'][str(L)] = compact(virtual_probe_limit(L))
    for offsets in ((0,0,0,0),(0,1,0,-1),(0,-1,0,1)):
        result['position'][str(offsets)] = compact(virtual_probe_limit(7,offsets=offsets))
    for eta in (1e-3,1e-4,1e-5):
        result['eta'][str(eta)] = compact(virtual_probe_limit(7,eta=eta))
    ref=virtual_probe_limit(7)
    for strength in (.1,.01,.001,.0001):
        d=probe_strength_response(7,probe_scale=strength)
        result['lambda'][str(strength)]={
            'Rxy_finite_h_over_q2':d['Rxy_h_over_q2'],
            'Rxy_minus_virtual':d['Rxy_h_over_q2']-ref['Rxy_lambda_to_0_h_over_q2'],
            'current_minus_virtual':d['drive_current_q2_over_h']-ref['two_terminal_T_lr'],
            'max_xcut_minus_current':float(max(abs(d['x_cut_currents']-d['drive_current_q2_over_h'])))}
    result['metal']=compact(virtual_probe_limit(5,mu=1.5))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='results/frozen/v64_virtual_probe_ordered_limits.json')
    a=parser.parse_args();f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(create(),sort_keys=True,separators=(',',':'))+'\n')
    print('Wrote',f)
