"""Gapped T=0 bilayer charge dipole and full static EB spin response.

Convention: chi_ij=-2 Im sum_occ,empty S_i(vc) V_j(cv)/(ev-ec)^2;
C_ij=d_Bz chi_ij, H(B)=H(0)-g mu_B B sigma_z (legacy effective-g convention).
No Fermi-surface transport, scattering model, or finite-T DC claim.
"""
import numpy as np
from nmee_lab.clean_core.core import MU_B
from .sliding_tight_binding import P_TOP, P_BOTTOM

LAYER = P_TOP-P_BOTTOM


def flatten_hamiltonian(h, velocities, gap):
    """Counterfactual with identical +/-gap/2 spectrum at EVERY k and stacking.

Keep the occupied projector P and physical spin/layer operators. Set
Hflat=gap*(I/2-P), Vflat=-gap*d_k P. Add Zeeman coupling only AFTER flattening.
This is a different, generally long-range Hamiltonian, not a decomposition of
the original response and not an experimentally available operation.
"""
    if not np.isfinite(gap) or gap <= 0:
        raise ValueError('positive flattening gap required')
    e,u = np.linalg.eigh(h)
    if np.min(e[:,4]-e[:,3]) <= 1e-12:
        raise ValueError('cannot flatten a closed occupied-empty gap')
    ud = u.conj().swapaxes(1,2)
    p = u[:,:,:4]@u[:,:,:4].conj().swapaxes(1,2)
    vflat = []
    for v in velocities:
        ve = ud@v@u
        dp = np.zeros_like(ve)
        dp[:,:4,4:] = ve[:,:4,4:]/(e[:,:4,None]-e[:,None,4:])
        dp[:,4:,:4] = dp[:,:4,4:].conj().swapaxes(1,2)
        vflat.append(-gap*(u@dp@ud))
    return gap*(np.eye(8)/2-p), vflat


def observables(model, k, weights, *, step_T=.125, flatten_gap_eV=None,
                batch=1024, compute_response=True):
    k, w = np.asarray(k,float), np.asarray(weights,float)
    if k.ndim != 2 or k.shape[1] != 2 or w.shape != (len(k),) or not len(k):
        raise ValueError('k and weights must have shapes (N,2) and (N,)')
    if not np.isfinite(k).all() or not np.isfinite(w).all() or np.any(w < 0):
        raise ValueError('finite k and nonnegative finite weights required')
    if not np.isfinite(step_T) or step_T <= 0 or batch < 1:
        raise ValueError('positive field step and batch required')
    spin = [model.spins['total_x'],model.spins['total_y']]
    sz = model.spins['total_z']
    c = np.zeros((2,2)); charge = 0.; energy = 0.
    gap = np.inf; mu_distance = np.inf; valence_max = -np.inf; conduction_min = np.inf
    for start in range(0,len(k),batch):
        ks, ws = k[start:start+batch], w[start:start+batch]
        h = model.matrices(ks)
        e,u = np.linalg.eigh(h)
        gap = min(gap,float(np.min(e[:,4]-e[:,3])))
        valence_max = max(valence_max,float(e[:,3].max()))
        conduction_min = min(conduction_min,float(e[:,4].min()))
        mu_distance = min(mu_distance,float((model.spec.mu_eV-e[:,3]).min()),
                          float((e[:,4]-model.spec.mu_eV).min()))
        if mu_distance <= 0:
            raise ValueError('fixed chemical potential does not lie in a gap')
        occupied = u[:,:,:4]
        imbalance = np.einsum('nai,ab,nbi->n',occupied.conj(),LAYER,occupied).real
        charge += float(ws@imbalance)
        energy += float(ws@e[:,:4].sum(axis=1))
        if not compute_response:
            continue
        velocities = [model.matrices(ks,axis=j) for j in range(2)]
        mu = model.spec.mu_eV
        if flatten_gap_eV is not None:
            h, velocities = flatten_hamiltonian(h,velocities,flatten_gap_eV)
            mu = 0.
        chis = []
        for sign in (1,-1):
            eb,ub = np.linalg.eigh(h-sign*model.mono.g_factor*MU_B*step_T*sz)
            if not np.all((eb[:,3]<mu)&(mu<eb[:,4])):
                raise ValueError('Zeeman field changes insulating occupations')
            ud = ub.conj().swapaxes(1,2)
            ss = [(ud@s@ub)[:,:4,4:] for s in spin]
            vv = [(ud@v@ub)[:,4:,:4].swapaxes(1,2) for v in velocities]
            den = (eb[:,:4,None]-eb[:,None,4:])**2
            chis.append(np.array([[-2*np.imag(s*v/den).sum(axis=(1,2))
                                   for v in vv] for s in ss]))
        c += np.einsum('ijn,n->ij',(chis[0]-chis[1])/(2*step_T),ws)
    area = abs(float(np.linalg.det(model.a)))
    # Neutral rigid layers: equal ionic + background charges at +/-d/2 cancel.
    # D_z/(e*d) per unit cell = -(N_top-N_bottom)/2; no arbitrary thickness.
    return {
        'C': c.tolist() if compute_response else None,
        'trace': float(np.trace(c)) if compute_response else None,
        'difference': float(c[0,0]-c[1,1]) if compute_response else None,
        'offdiagonal_sum': float(c[0,1]+c[1,0]) if compute_response else None,
        'delta_N_per_cell': charge*area,
        'dipole_over_e_d_per_cell': -.5*charge*area,
        'dipole_over_e_d_per_nm2': -.5*charge,
        'occupied_energy_eV_per_cell': energy*area,
        'sampled_direct_gap_eV': gap,
        'sampled_indirect_gap_eV': conduction_min-valence_max,
        'sampled_mu_distance_eV': mu_distance,
    }
