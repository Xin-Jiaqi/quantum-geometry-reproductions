"""Validate same stacking in physical Cartesian units, including periodic wrapping."""
import numpy as np

def validate_fixed_registry(initial,relaxed,tolerance_A=1e-5):
    if initial.get_chemical_symbols()!=relaxed.get_chemical_symbols():
        raise ValueError('Atomic identities/order changed')
    if np.max(abs(initial.cell.array-relaxed.cell.array))>1e-7:
        raise ValueError('Prescribed strained cell changed')
    df=relaxed.get_scaled_positions()[:,:2]-initial.get_scaled_positions()[:,:2]
    df-=np.rint(df)
    dr=df@initial.cell.array[:2]
    maximum=float(np.linalg.norm(dr,axis=1).max())
    if maximum>tolerance_A:
        raise ValueError(f'Registry changed by {maximum:g} A (> {tolerance_A:g} A)')
    return {'max_lateral_displacement_A':maximum,'max_lateral_fractional_change':float(abs(df).max()),'tolerance_A':tolerance_A}
