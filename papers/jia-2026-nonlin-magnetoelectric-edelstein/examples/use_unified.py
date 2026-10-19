"""Zero-external-data demonstration; reproduces the validated calculation paths."""
from nmee_lab import CleanSpec, CuMnAsBGKSpec, evaluate_clean, evaluate_cumnas_bgk, mesh
from nmee_lab.clean_core.models import kane_mele

if __name__=='__main__':
    km=kane_mele()
    p=mesh(km.reciprocal,451)
    result=evaluate_clean(km,p,CleanSpec(alpha='x',gamma='z',mu_eV=.025,temperature_K=10))
    print('Kane-Mele clean full_static:',result['result'])
    print('Kane-Mele printed branch:   ',result['parts']['printed'])
    print('Diagonal term:              ',result['parts']['diagonal_disputed'])
    result=evaluate_cumnas_bgk(CuMnAsBGKSpec(
              observable='A_z',magnetic_component='x',electric_beta=0),n=51)
    print('CuMnAs A sublattice, BGK:  ',result['result'])
    print('CuMnAs units:',result['normalization'])
