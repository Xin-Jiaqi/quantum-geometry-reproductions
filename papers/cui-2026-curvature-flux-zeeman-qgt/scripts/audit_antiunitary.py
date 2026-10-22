"""Standalone synthetic antiunitary microscopic Kubo/Zeeman audit (no DFT)."""
import json
import numpy as np
from flux_lab.antiunitary_response import (examples,toy_bloch,operator_residual,
 qgt_symmetry_residual,kubo_symmetry_residual,wrong_retarded_branch_error,
 fourband_antiunitary_residual)
from flux_lab.causal_response import interband_kubo

def bz(kind,n=16,w=.035):
    ks=-np.pi+(np.arange(n)+.5)*2*np.pi/n
    positive=np.zeros((2,3),complex);negative=positive.copy()
    for x in ks:
        for y in ks:
            H,V,S=toy_bloch([x,y],kind)
            positive+=interband_kubo(H,V,S,z=w,mu=.1,temperature=.16)[0]
            negative+=interband_kubo(H,V,S,z=-w,mu=.1,temperature=.16)[0]
    odd=(positive-negative)/(2*n*n);even=(positive+negative)/(2*n*n)
    return dict(odd_real=odd.real.tolist(),odd_imag=odd.imag.tolist(),
                even_real=even.real.tolist(),even_imag=even.imag.tolist(),
                max_odd=float(np.max(abs(odd))),max_even=float(np.max(abs(even))))

def run(n=16):
    k=np.array([.43,-.38]);output={}
    for name,spec in examples().items():
        q=kubo_symmetry_residual(name,k,z=.07+.014j)
        l=kubo_symmetry_residual(name,k,z=.07+.014j,method='liouville')
        output[name]=dict(group=spec.name,spin_parity=spec.spin_parity.tolist(),
            momentum_map=spec.R.tolist(),operators=operator_residual(name,k),
            qgt_pair_max_residual=qgt_symmetry_residual(name,k),
            spectral_residual={a:q[a] for a in ('chi','alpha')},
            independent_liouville_residual={a:l[a] for a in ('chi','alpha')},
            fourband_spin_orbital_residual=fourband_antiunitary_residual(name,k),
            wrong_fixed_retarded_branch=wrong_retarded_branch_error(name,k),
            BZ=bz(name,n))
    return dict(source='Cui et al. 2026 v3 Eq.(16)-(19), model-only',
                momentum=k.tolist(),complex_drive_frequency=[.07,.014],
                BZ_kmesh=n,models=output,full_paper_status='PARTIAL')

if __name__=='__main__':print(json.dumps(run(),indent=2))
