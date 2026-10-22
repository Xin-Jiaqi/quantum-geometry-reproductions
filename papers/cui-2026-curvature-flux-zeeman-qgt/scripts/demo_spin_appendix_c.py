"""Only a synthetic four-band conventional spin-current Appendix C test."""
import json, argparse
from flux_lab.spin_current import spin_orbital_benchmark,intrinsic_spin_geometry,integrate_spin_rta_bz

def run(mesh=40):
    H,v,K,s=spin_orbital_benchmark(.31,-.22,periodic=True)
    d=intrinsic_spin_geometry(H,v,s)
    bz=integrate_spin_rta_bz(mesh=mesh)
    return dict(k=[.31,-.22],spin='I_tau tensor sigma_z',
        intrinsic_yz_E_x_band0=float(d['intrinsic'][0,1,0]),
        kmesh=mesh,
        direct_fermi_surface=float(bz['direct_fermi_surface'][0,0]),
        full_fermi_sea=float(bz['full_fermi_sea'][0,0]),
        reduced_geometry_only=float(bz['geometry_only'][0,0]),
        missing_hessian=float(bz['missing_hessian'][0,0]),
        units='natural, tau=1, synthetic conventional spin current')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--kmesh',type=int,default=40)
    print(json.dumps(run(p.parse_args().kmesh),indent=2))
