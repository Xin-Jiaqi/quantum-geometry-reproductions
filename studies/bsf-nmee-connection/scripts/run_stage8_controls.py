"""Targeted asymptotic and convergence controls accompanying the main scan."""
import json
from dataclasses import asdict, replace
from run_stage8_audit import ROOT, at
from bsf_nmee.bilayer_audit import AuditBilayer, integration_grid, symmetry_residuals
from bsf_nmee.static_observables import observables


def run():
    grids = {name: integration_grid(AuditBilayer(), *shape)[:2]
             for name, shape in [('base',(61,60,96)),('fine',(101,100,144))]}
    def evaluate(spec, grid='base', flat=None, step=.125):
        k,w = grids[grid]
        return {'parameters':asdict(spec),'grid':grid,'step_T':step,
                'flatten_gap_eV':flat,
                **observables(AuditBilayer(spec),k,w,step_T=step,flatten_gap_eV=flat)}
    result = {'scope':'targeted controls; same conventions as bilayer_stage8_audit.json'}
    result['symmetry'] = symmetry_residuals(at(.1))
    result['broken_C3_asymptotic'] = []
    for u in (.003125,.00625,.0125):
        result['broken_C3_asymptotic'].append(
            {'u':u,**evaluate(at(u,bond_anisotropy_eV=.008))})
        print('broken C3',u,flush=True)
    result['broken_C3_fine'] = evaluate(at(.003125,bond_anisotropy_eV=.008),'fine')
    result['same_Rashba_origin'] = evaluate(at(0.,stacking='direct'))
    result['flat_convergence'] = []
    for t in (.005,.02):
        spec=replace(at(.1),tunneling_eV=t)
        result['flat_convergence'].append({'t_eV':t,'u':.1,
            'fine':evaluate(spec,'fine',flat=.03),
            'half_B':evaluate(spec,flat=.03,step=.0625)})
        print('flat convergence',t,flush=True)
    output=ROOT/'results/bilayer_stage8_controls.json'
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    run()
