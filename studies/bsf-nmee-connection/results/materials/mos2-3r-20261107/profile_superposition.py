"""Compare the interface with a rigid sum of isolated-layer response profiles."""
import json
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parent
a=json.loads((root/'profiles/548484.json').read_text())
centers=np.array(a['ion_z_A'][:2])
actual=np.array(a['slice_response'])
rows=[]
for jid in ['548482','548488','548489']:
    m=json.loads((root/f'profiles/{jid}.json').read_text())
    z=np.array(m['z_A']);v=np.array(m['slice_response'])
    dz=z[1]-z[0];length=dz*len(z)
    w=(z<0).astype(float);w[[0,len(z)//2]]=.5
    linear=sum(np.interp(z-c,z,v,period=length) for c in centers)
    frequency=np.fft.fftfreq(len(z),d=dz)
    fourier=sum(np.fft.ifft(np.fft.fft(v)*np.exp(-2j*np.pi*frequency*c)).real for c in centers)
    def integrate(x):return [float(x.sum()),float(w@x),float((1-w)@x)]
    row={'gate_physical_Vnm':m['gate_physical_Vnm'],
         'rigid_sum_total_bottom_top':integrate(fourier),
         'linear_translation_total_bottom_top':integrate(linear)}
    if m['gate_physical_Vnm']==0:
        row.update(actual_total_bottom_top=integrate(actual),
                   interface_change_total_bottom_top=integrate(actual-fourier))
    rows.append(row)
result={'Mo_centers_A':centers.tolist(),'rows':rows,
        'units':'1e-9 (nm V T)^-1',
        'note':'Rigidly translated isolated-layer response densities, not a self-consistent decoupled-bilayer calculation. Half-space values include density tails from both layers; the total is partition independent. Deviations include all interface effects and do not uniquely isolate hybridization.'}
(root/'profile_superposition.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
