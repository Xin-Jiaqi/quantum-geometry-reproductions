"""Compare independently gate-self-consistent, fixed-geometry NMEE responses.
Runs on local, collected small JSON only; not on the Japan login frontend.
"""
from pathlib import Path
import csv,json
import numpy as np
root=Path(__file__).resolve().parent
data={case:json.loads((root/case/'result.json').read_text()) for case in ('fm010','fp010')}
zero=json.loads((root/'zero_reference.json').read_text())
assert all(x['complete'] and x['stage']=='complete' for x in data.values())
assert data['fm010']['F_external_V_nm']==-.1 and data['fp010']['F_external_V_nm']==.1
points=[];slopes={};audit={'half_probe_errors':[],'geometry_fixed_max_error':{},'checks':[]}
for material in ('A','mono','upper'):
    lo=data['fm010']['materials'][material];hi=data['fp010']['materials'][material]
    assert lo['scf_converged'] and hi['scf_converged']
    error=float(np.max(abs(np.array(lo['positions_A'])-hi['positions_A'])))
    audit['geometry_fixed_max_error'][material]=error;assert error<1e-10
    assert lo['input_sha256']['POSCAR']==hi['input_sha256']['POSCAR']
    electrons=48 if material=='A' else 24
    charges=[float(np.mean(x['planar_charge_raw'])) for x in (lo,hi)]
    assert max(abs(q-electrons) for q in charges)<1e-4
    audit.setdefault('total_electrons',{})[material]=charges
    slopes[material]={}
    for nb in lo['response']:
        cm=np.array(lo['response'][nb]['standard']['C_total_lower_upper'])
        cp=np.array(hi['response'][nb]['standard']['C_total_lower_upper'])
        for case,row,c in [('fm010',lo,cm),('fp010',hi,cp)]:
            ans=row['response'][nb]['standard'];assert ans['gap_eV']>.05
            assert abs(c[0]-c[1]-c[2])<1e-10
            half=np.array(row['response'][nb]['half_b']['C_total_lower_upper'])
            physical_error=abs(c-half) if material=='A' else abs(c[:1]-half[:1])
            audit['half_probe_errors'].append(float(physical_error.max()))
            # In an isolated layer the z=cell/2 partition cuts THROUGH the layer.
            # Only its total response is a physical isolated-layer reference.
            points.append({'material':material,'nb':int(nb),'F_V_nm':data[case]['F_external_V_nm'],
                           'C_net':c[0],'C_lower_layer':c[1] if material=='A' else None,
                           'C_upper_layer':c[2] if material=='A' else None,'gap_eV':ans['gap_eV']})
        slope=(cp-cm)/.2
        slopes[material][nb]=slope.tolist()
        if material in zero:
            c0=np.array(zero[material]['results'][nb]['standard']['C_total_lower_upper'])
            audit.setdefault('central_even_shift',{}).setdefault(material,{})[nb]=((cp+cm)/2-c0).tolist()
    # Scalar potential's induced slope is measured, not assumed from VASP input.
    dv=(np.array(hi['planar_potential_eV'])-lo['planar_potential_eV'])/.2
    L=lo['cell_A'][2][2];z=np.arange(len(dv))*L/len(dv)
    # Exclude the periodic sawtooth boundary and any atomic region.
    vac=(z>2)&(z<4)
    fit=np.polyfit(z[vac],dv[vac],1)
    audit.setdefault('potential',{})[material]={'vacuum_dV_dz_per_external_V_nm':float(fit[0]),
        'vacuum_fit_residual_eV_per_V_nm':float(np.std(dv[vac]-np.polyval(fit,z[vac])))}
    for atom,pos in zip(lo['symbols'],lo['positions_A']):
        if atom=='Mo':
            mask=abs(z-pos[2])<.25
            audit['potential'][material].setdefault('Mo_plane_dV_per_F',[]).append(float(dv[mask].mean()))
    if material=='A':
        zz=[p[2] for atom,p in zip(lo['symbols'],lo['positions_A']) if atom=='Mo']
        mv=audit['potential'][material]['Mo_plane_dV_per_F']
        ordered=sorted(zip(zz,mv));lever=(ordered[1][1]-ordered[0][1])
        audit['potential'][material]['interlayer_potential_lever_eV_per_V_nm']=lever
        audit['potential'][material]['bare_linear_lever_eV_per_V_nm']=(ordered[1][0]-ordered[0][0])/10

comparisons=[]
for nb in ('96','112'):
    g=np.array(slopes['A'][nb])
    for mono_nb in ('48','56','64'):
        gl=slopes['mono'][mono_nb][0];gu=slopes['upper'][mono_nb][0];ref=gl+gu
        assert abs(ref)>1e-5,'Reference slope too small for relative gain'
        comparisons.append({'bilayer_nb':int(nb),'mono_nb':int(mono_nb),'G_net':g[0],
            'G_lower':g[1],'G_upper':g[2],'G_isolated_lower':gl,'G_isolated_upper':gu,
            'G_isolated_pair':ref,'absolute_gain_fraction':float(abs(g[0])/abs(ref)-1),
            'upper_absolute_share':float(abs(g[2])/(abs(g[1])+abs(g[2]))),
            'layer_slopes_same_sign':bool(g[1]*g[2]>0)})
for nb in ('48','56','64'):
    gl=slopes['mono'][nb][0];gu=slopes['upper'][nb][0]
    audit.setdefault('exact_upper_vs_mirrored_lower_slope',{})[nb]={
        'difference':gu-gl,'fraction_of_pair':float(abs(gu-gl)/abs(gu+gl))}
best=[r for r in comparisons if r['mono_nb']==64]
result={'complete':True,'scope':'perpendicular-gate SCF; fixed-ion, spin-only frozen probe response',
        'units_C':'1e-9 (nm V T)^-1','units_G':'1e-9 (V^2 T)^-1',
        'fields_external_V_nm':[-.1,.1],'slopes':slopes,'comparisons':comparisons,
        'gain_high_cutoff_range':[min(r['absolute_gain_fraction'] for r in best),max(r['absolute_gain_fraction'] for r in best)],
        'upper_share_range':[min(r['upper_absolute_share'] for r in best),max(r['upper_absolute_share'] for r in best)],
        'min_gap_eV':min(r['gap_eV'] for r in points),'audit':audit,
        'same_layer_reversal_target':False,'k_mesh':'6x6 inherited quick mechanism test, not final mesh convergence'}
with (root/'points.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(points[0]),lineterminator='\n');w.writeheader();w.writerows(points)
with (root/'slopes.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(comparisons[0]),lineterminator='\n');w.writeheader();w.writerows(comparisons)
(root/'summary.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k] for k in ['gain_high_cutoff_range','upper_share_range','min_gap_eV','slopes']},indent=2))
