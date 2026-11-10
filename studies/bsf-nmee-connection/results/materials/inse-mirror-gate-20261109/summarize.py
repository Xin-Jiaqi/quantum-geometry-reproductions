"""Audit independent mirror-domain gate data; no fitting or synthetic points."""
from pathlib import Path
import csv, json
import numpy as np

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'inse-nmee-20261109'
POINTS = {
    'A06': OLD/'A0/result.json', 'Am6': OLD/'Am/result.json',
    'Ap6': OLD/'Ap/result.json', 'B06': OLD/'B0/result.json',
    'A08': OLD/'A0_mesh8/result.json',
    **{s: ROOT/s/'result.json' for s in ('Bm6','Bp6','B08','Bp8','Ap8')},
}

def build():
    missing = [str(p) for p in POINTS.values() if not p.exists()]
    if missing:
        raise SystemExit('Missing independent results: ' + ', '.join(missing))
    data = {s: json.loads(p.read_text()) for s,p in POINTS.items()}
    for s,d in data.items():
        assert d['complete'], s
        if s=='A08':
            # Earlier mesh_check.py reused A0 INCAR (NSW=0); JSON predates fixed_ions.
            assert d['input_sha256']['POSCAR']==data['A06']['input_sha256']['POSCAR']
        else:
            assert d['fixed_ions'], s
        assert d['label'] == s[0], s
        assert d['export']['nk'] == int(s[-1])**2, s
        assert abs(d['F_external_V_nm'] - {'0':0.,'m':-.1,'p':.1}[s[1]]) < 1e-12, s
    for domain in ('A','B'):
        hashes={d['input_sha256']['POSCAR'] for d in data.values() if d['label']==domain}
        assert len(hashes)==1, ('Geometry changed',domain,hashes)
    rows, deltas, mirrors, steps = [], [], [], []
    for s,d in data.items():
        for cutoff in ('128','144','160'):
            for step in ('standard','half_b'):
                r = d['response'][cutoff][step]
                c = np.array(r['C_total_lower_upper'])
                assert np.isfinite(c).all() and r['gap_eV'] > .05, (s,cutoff,step)
                assert abs(c[0]-c[1]-c[2]) < 1e-7, (s,cutoff,step)
                rows.append(dict(point=s,domain=s[0],mesh=int(s[-1]),F_V_nm=d['F_external_V_nm'],
                                 cutoff=int(cutoff),step=step,C_net=c[0],C_lower=c[1],C_upper=c[2],
                                 C_staggered=(c[2]-c[1])/2,gap_eV=r['gap_eV']))
            cs=np.array(d['response'][cutoff]['standard']['C_total_lower_upper'])
            ch=np.array(d['response'][cutoff]['half_b']['C_total_lower_upper'])
            steps.append({'point':s,'cutoff':int(cutoff),'C_standard_minus_half_b':(cs-ch).tolist()})
    def C(s,cutoff='160',step='standard'):
        return np.array(data[s]['response'][cutoff][step]['C_total_lower_upper'])
    for mesh in (6,8):
        for domain in ('A','B'):
            for cutoff in ('128','144','160'):
                for step in ('standard','half_b'):
                    c0=C(f'{domain}0{mesh}',cutoff,step); cp=C(f'{domain}p{mesh}',cutoff,step)
                    dc=cp-c0
                    deltas.append(dict(domain=domain,mesh=mesh,cutoff=int(cutoff),step=step,
                        delta_C_net=dc[0],delta_C_lower=dc[1],delta_C_upper=dc[2],
                        delta_abs_C_lower=abs(cp[1])-abs(c0[1]),delta_abs_C_upper=abs(cp[2])-abs(c0[2]),
                        delta_C_staggered=(dc[2]-dc[1])/2,
                        fractional_abs_change_lower=(abs(cp[1])-abs(c0[1]))/abs(c0[1]),
                        fractional_abs_change_upper=(abs(cp[2])-abs(c0[2]))/abs(c0[2])))
    for b,a in (('B06','A06'),('Bp6','Am6'),('Bm6','Ap6'),('B08','A08')):
        for cutoff in ('128','144','160'):
            cb=C(b,cutoff); ca=C(a,cutoff); err=cb+ca[[0,2,1]]
            mirrors.append(dict(B=b,A=a,cutoff=int(cutoff),absolute_error=err.tolist(),
                layer_relative_max=float(np.max(abs(err[1:])/abs(ca[[2,1]])))))
    # Difference of changes, not the large raw local C, determines confidence in tuning.
    derivative_step_errors=[]
    for mesh in (6,8):
        for domain in ('A','B'):
            for cutoff in ('128','144','160'):
                dd=(C(f'{domain}p{mesh}',cutoff)-C(f'{domain}0{mesh}',cutoff))-(C(f'{domain}p{mesh}',cutoff,'half_b')-C(f'{domain}0{mesh}',cutoff,'half_b'))
                derivative_step_errors.append(dict(domain=domain,mesh=mesh,cutoff=int(cutoff),delta_step_error=dd.tolist()))
    slopes={domain:{cutoff:((C(f'{domain}p6',cutoff)-C(f'{domain}m6',cutoff))/.2).tolist()
                    for cutoff in ('128','144','160')} for domain in ('A','B')}
    forces={s:{'layers':d.get('layer_forces_eV_A'), 'scf_seconds':d.get('scf_seconds',d.get('vasp_seconds')),
               'iterations':d.get('electronic_iterations'), 'job':d['job'],'node':d['node']} for s,d in data.items()}
    summary={'units_C':'1e-9 (nm V T)^-1, Pauli-spin areal density coefficient',
      'units_G':'1e-9 (V^2 T)^-1','positive_gate_V_nm':.1,
      'results':rows,'positive_gate_changes':deltas,'mirror_errors':mirrors,
      'central_slopes_6':slopes,'step_checks_C':steps,'step_checks_delta_C':derivative_step_errors,
      'metadata':forces,'all_layers_lower_negative_upper_positive':all(r['C_lower']<0<r['C_upper'] for r in rows),
      'A_both_weaken_all_settings':all(d['delta_abs_C_lower']<0 and d['delta_abs_C_upper']<0 for d in deltas if d['domain']=='A'),
      'B_both_strengthen_all_settings':all(d['delta_abs_C_lower']>0 and d['delta_abs_C_upper']>0 for d in deltas if d['domain']=='B'),
      'minimum_gap_eV':min(r['gap_eV'] for r in rows),
      'max_mirror_local_absolute_error':max(max(abs(x) for x in d['absolute_error'][1:]) for d in mirrors),
      'max_delta_C_probe_step_error':max(max(abs(x) for x in d['delta_step_error']) for d in derivative_step_errors),
      'limits':['Fixed ionic coordinates, no switching threshold or relaxation calculated.',
                'Two meshes test trend; no fully converged enhancement-factor claim.',
                'Different cutoffs/steps are deterministic checks, not statistical replicates.',
                'No 8x8 negative-field points; central G only available at 6x6.',
                'Layer sum, mirror and step checks do not by themselves establish completeness of the physical response.']}
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    for name,records in (('source_data',rows),('positive_gate_changes',deltas)):
        with (ROOT/f'{name}.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator="\n");writer.writeheader();writer.writerows(records)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['results','positive_gate_changes','mirror_errors','step_checks_C','step_checks_delta_C','metadata']},indent=2))
    print('Representative positive-gate changes:')
    for d in deltas:
        if d['cutoff']==160 and d['step']=='standard':print(d)
    print('Maximum mirror layer error:',max(d['layer_relative_max'] for d in mirrors))
    return summary

if __name__=='__main__':build()
