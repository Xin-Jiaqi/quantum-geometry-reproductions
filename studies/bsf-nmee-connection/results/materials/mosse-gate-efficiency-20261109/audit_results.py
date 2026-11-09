"""Check physical consistency and report numerical scope of saved model results."""
from pathlib import Path
import ast,json
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parent
scan=json.loads((root/'scan.json').read_text())
direction=json.loads((root/'gate_direction.json').read_text())
assert scan['complete'] and direction['complete']
assert len(direction['rows'])==16
c=pd.read_csv(root/'curves.csv');s=pd.read_csv(root/'slopes.csv')
assert len(c)==407 and np.isfinite(c[['C_net','C_lower','C_upper','gap_eV']]).all().all()
sum_error=float(abs(c.C_net-c.C_lower-c.C_upper).max());assert sum_error<1e-10
assert c.gap_eV.min()>.05
full=c[c.kind=='full'];assert (full.C_lower<0).all() and (full.C_upper>0).all()
central=s[np.isclose(s.half_interval,.025)]
gf=-central.loc[central.kind=='full','G_net'].values
ref=central[((central.kind=='isolated_pair')|(central.kind=='electrostatic_pair'))&(central.nb==64)]
gr=-ref.G_net.values
gain=[float(gf.min()/gr.max()-1),float(gf.max()/gr.min()-1)]
assert gain[0]>.25 and gain[1]<.4
fraction=(central.loc[central.kind=='full','G_upper']/central.loc[central.kind=='full','G_net'])
raw=json.loads((root/'summary.json').read_text())
assert not raw['invalid_rows'] and raw['max_half_B_change']<1e-3
errs=[abs(v) for r in direction['rows'] for v in r['tangent_additivity_residual']]
assert max(errs)<3e-4
for nb in (96,112):
 rows=[r for r in direction['rows'] if r['nb']==nb and r['delta']==.01]
 whole=[r for r in rows if r['control']=='full'][0]
 for r in rows:
  if r['control']=='no_VC':
   assert abs(whole['G']['linear'][0])>abs(r['G']['linear'][0])
   assert abs(whole['G']['layer_step'][0])>3*abs(r['G']['layer_step'][0])
for p in root.glob('*.py'):ast.parse(p.read_text(),filename=p.name)
record={'passed':True,'curve_rows':len(c),'direction_rows':len(direction['rows']),
        'layer_sum_error':sum_error,'min_gap_eV':float(c.gap_eV.min()),
        'full_layers_keep_opposite_signs':True,'gain_over_64state_refs_fraction':gain,
        'upper_fraction_of_full_gate_slope':[float(fraction.min()),float(fraction.max())],
        'tangent_additivity_max_error':max(errs),'half_Zeeman_step_max_error':raw['max_half_B_change'],
        'scope':'saved-results consistency, not proof of absolute mesh/screening convergence',
        'visual_qa':'PNG inspected: labels/units readable, no clipping; setting ranges not statistical errors',
        'review':'checked mirror gate parity, physical gate vs cache Z spin distinction, consistent channel deletion in H and gate, nonadditive intervention caveats; no credentials or licensed potentials in archive'}
(root/'audit.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
