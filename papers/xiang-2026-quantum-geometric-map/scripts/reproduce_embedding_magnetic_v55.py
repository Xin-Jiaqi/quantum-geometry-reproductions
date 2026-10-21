#!/usr/bin/env python3
"""v55 orbital connection/H(k) or finite open graph explicit Peierls B coupling."""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.embedding_magnetic_v55 import (
 orbital_center_connection,connected_band_geometry,peierls_field,
 field_gauge_unitary,magnetic_vertex_frame_law)
from qgmap.orbital_geometry_v54 import periodic_three_band_model

def evaluate(mode='reference',bloch_npz=None,graph_npz=None,field=.32,charge_over_hbar=-.87):
    if mode not in ('reference','bloch','peierls'):raise ValueError('mode')
    if mode in ('reference','bloch'):
        if mode=='reference':
            k,H,V=periodic_three_band_model(12)
            H=H[4,7];dH=V[:,4,7]
            A=orbital_center_connection([[0.,0.],[.19,-.25],[-.32,.47]])
            source='analytic_three_band_with_point_orbital_centers'
        else:
            if bloch_npz is None:raise ValueError('--bloch-npz required')
            with np.load(bloch_npz,allow_pickle=False) as f:
                if not {'hamiltonian','derivatives','orbital_connection'}<=set(f.files):
                    raise ValueError('NPZ requires hamiltonian derivatives orbital_connection')
                H=f['hamiltonian'];dH=f['derivatives'];A=f['orbital_connection']
            source='external_complete_H_dH_A'
        g=connected_band_geometry(H,dH,A)
        return {'status':g.acceptance,'source':source,'min_direct_gap':g.min_direct_gap,
                'energies':np.asarray(g.energies).tolist(),
                'berry_curvature':np.asarray(g.curvature).tolist(),
                'orbital_moment_energy':np.asarray(g.moment_energy).tolist(),
                'metric_xy':np.asarray(g.metric)[...,0,1].tolist(),
                'scope':'Nondegenerate complete finite-orbital connection required; no total magnetic transport'}
    if graph_npz is None:raise ValueError('--graph-npz required')
    with np.load(graph_npz,allow_pickle=False) as f:
        if not {'hoppings','positions'}<=set(f.files):
            raise ValueError('NPZ requires hoppings and positions')
        H0=f['hoppings'];positions=f['positions']
    s=peierls_field(H0,positions,field,charge_over_hbar)
    l=peierls_field(H0,positions,field,charge_over_hbar,gauge='landau_y')
    G,dG=field_gauge_unitary(positions,field,charge_over_hbar)
    compared=magnetic_vertex_frame_law(s.hamiltonian,s.magnetic_vertex,G,dG)
    return {'status':s.acceptance,'nodes':len(positions),'field':field,'q_over_hbar':charge_over_hbar,
            'max_gauge_spectrum_error':float(np.max(abs(np.linalg.eigvalsh(s.hamiltonian)-np.linalg.eigvalsh(l.hamiltonian)))),
            'max_magnetic_vertex_gauge_residual':float(np.max(abs(compared-l.magnetic_vertex))),
            'max_magnetic_vertex_element':float(np.max(abs(s.magnetic_vertex))),
            'scope':'Finite open graph and straight link path; NOT periodic Bloch B or EB conductivity'}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['reference','bloch','peierls'],default='reference')
    p.add_argument('--bloch-npz',type=Path)
    p.add_argument('--graph-npz',type=Path)
    p.add_argument('--B',type=float,default=.32)
    p.add_argument('--q-over-hbar',type=float,default=-.87)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    result=evaluate(a.mode,a.bloch_npz,a.graph_npz,a.B,a.q_over_hbar)
    text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text,encoding='utf8')
    else:print(text,end='')
if __name__=='__main__':main()
