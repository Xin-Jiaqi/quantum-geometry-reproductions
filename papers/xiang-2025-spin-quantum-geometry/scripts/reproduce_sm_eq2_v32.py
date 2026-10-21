"""Reproduce literal SM Eq.(2) source groups and Eq.(6) difference."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from methods.sm_eq2_v32 import probe
if __name__=='__main__':
    print(json.dumps({'paper':'PRL 135 146303 (2025)',
      'scientific_status':'PARTIAL',
      'protocol':'adiabatic zero-frequency SM Eq2, 1st eta, 2nd 2eta, nondegenerate',
      'cases':[probe(n) for n in (2,3,4,5)],
      'Eq2_minus_printed_Eq6':'-i*r_mn*(df_n-df_m)/(epsilon_m-epsilon_n)**2',
      'intrinsic_metallic_DC':'NOT_VALIDATED'},sort_keys=True,indent=2))
