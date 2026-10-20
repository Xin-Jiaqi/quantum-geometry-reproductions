#!/usr/bin/env python3
"""Compare independent twisted QWZ-SCBA finite differences with RR current ladder."""
import argparse,json
from qgmap.rr_spatial_ward_v40 import audit
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,default=21)
    p.add_argument('--energy',type=float,default=1.1)
    p.add_argument('--W',type=float,default=.12)
    p.add_argument('--eta',type=float,default=.12)
    p.add_argument('--step',type=float,default=2e-4)
    a=p.parse_args()
    print(json.dumps(audit(mesh=a.mesh,energy=a.energy,W=a.W,eta=a.eta,step=a.step),indent=2,allow_nan=False))
if __name__=='__main__':main()
