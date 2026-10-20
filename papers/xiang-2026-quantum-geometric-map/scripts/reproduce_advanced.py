#!/usr/bin/env python3
"""Run named, protocol-separated model reproductions; not a real-material solver."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from qgmap.spin_multiband_v5 import four_band_audit
from qgmap.orbital_phe_v4 import orbital_phe_eq56
from qgmap.orbital_mnhe_v5 import orbital_mnhe_2d,orbital_mnhe_3d
from qgmap.warped_disorder_v17 import Lattice,ra_observable,density_ward,slope as ra_slope
from qgmap.spectral_kubo_v18 import spectral_bubble,slope as spectral_slope
from qgmap.finite_qomega_ward_v19 import finite_ward
from qgmap.real_frequency_continuity_v20 import ExactSystem,kinematic_ward,uniform_density_vs_compressibility
from qgmap.limit_order_v20 import scan_limit
from qgmap.electromagnetic_eb_v21 import EBSpec, liouville_slope, finite_zeeman_slope, mixed_ward, gauge_mixed
from qgmap.two_color_exact_v22 import TwoColorSpec, frequency_channels, independent_sylvester_channels, static_preparation
from qgmap.homogeneous_length_gauge_v23 import UniformSpec, source_and_channels, algebraic_channels, leading_fermi_pole, eq20_direct_band_current
from qgmap.continuum_ti_v24 import DiracSpec, projector_channels, eigen_channels, eq20_direct, angular_exact
from qgmap.si_fig2_v25 import FigureInputs, evaluate as fig2_SI_evaluate
from qgmap.physical_si_liouville_v26 import SIConfig, eq20_physical, eq11_SI_magnitude, frequency_hypotheses
from qgmap.fig2_source_audit_v27 import load_points, fit_branch, ratio_to_literal, calibration_sensitivity

from dataclasses import replace



def _spin(n): return four_band_audit(N=n)
def _phe(n): return orbital_phe_eq56(N=n)
def _mnhe2d(n): return orbital_mnhe_2d(N=n)
def _mnhe3d(n): return orbital_mnhe_3d(N=n)

def _scba_ra(n):
    p=Lattice(n=n)
    r=ra_observable(p,b=0.,protocol="onsite",vertex=True)
    onsite=ra_slope(p,step=.02,protocol="onsite",vertex=True)
    shift=ra_slope(p,step=.02,protocol="shift",vertex=True)
    ward=density_ward(p,b=0.,protocol="onsite")
    return {"protocol":"SCBA_RA_finite_eta_periodic_lattice",
            "model":"independent Wilson/warped periodic two-band regulator",
            "units":"dimensionless BZ mean; b is model energy or shift, not Tesla",
            "at_zero":{k:r[k] for k in ("symmetric","xy","yx","vertex_residual","scba_residual","scba_iterations")},
            "onsite_zeeman_slope":onsite,
            "shift_control_slope":shift,
            "q0_density_ward":{k:ward[k] for k in ("abs_error","residual","iterations","causality_max_eigenvalue")},
            "limit":"RA Fermi-surface bubble only; not full Kubo-Streda or published PHE"}

def _spectral(n):
    p=Lattice(n=n)
    r=spectral_bubble(p,b=0.,protocol="onsite",vertex=True)
    s=spectral_slope(p,step=.02,protocol="onsite",vertex=True)
    return {"protocol":"SCBA_symmetric_spectral_Kubo_four_channels",
            "model":"independent periodic two-band scalar-disorder regulator",
            "units":"dimensionless BZ mean; common e^2*hbar/pi and k measure omitted",
            "at_zero":r,"onsite_zeeman_slope":s,
            "limit":"RA+AR-RR-AA with matching noncrossing ladders; not full Kubo-Streda"}

def _ward(n):
    p=Lattice(n=n)
    r=finite_ward(p,q=.08,omega=.04,direction=1,b=0.,protocol="onsite")
    return {"protocol":r["protocol"],"model":"independent periodic two-band scalar-disorder regulator",
            "parameters":{k:r[k] for k in ("q","omega","eta","W","direction","field_protocol","b")},
            "residuals":r["residuals"],
            "vertex_norms":{k:r[k] for k in ("density_ladder_correction_norm","longitudinal_ladder_correction_norm","density_vertex_norm")},
            "iterations":r["iterations"],"warnings":r["warnings"],
            "limit":"finite-eta resolvent Ward identity; not physical continuity at eta>0"}


def _complex_dict(z):
    return {"real":float(z.real),"imag":float(z.imag)}

def _exact_kubo_v20(L):
    system=ExactSystem(L=L)
    axes={}
    for direction,name in ((0,"x"),(1,"y")):
        out=kinematic_ward(system,direction=direction)
        axes[name]={
            "ward_density_abs":out["ward_density_abs"],
            "ward_current_abs":out["ward_current_abs"],
            "commutator_abs":out["commutator_abs"],
            "contact_rho_j":_complex_dict(out["contact_rho_j"]),
            "missing_contact_abs":out["bare_missing_current_contact_abs"],
            "q":out["q"],"omega":out["omega"],"delta":out["delta"],
        }
    thermo=uniform_density_vs_compressibility(system)
    return {"protocol":"exact_finite_system_real_frequency_longitudinal_Kubo_continuity",
            "L":L,"hilbert_space_dim":2*L*L,
            "directions":axes,
            "thermodynamic":{
                "uniform_dynamic_density":_complex_dict(thermo["dynamic_uniform_density"]),
                "grand_canonical_compressibility":thermo["equilibrium_compressibility"],
                "compressibility_finite_difference":thermo["compressibility_finite_difference"],
            },
            "limits":"finite size, exact static disorder, longitudinal charge continuity; not PHE or generic open-system DC"}

def _limit_order_v20(n):
    result=scan_limit(n=n,eta_values=(.02,.005,.001),q=.06,omega=.04)
    return {"protocol":result["protocol"],"n":n,"q":result["q"],"omega":result["omega"],
            "values":result["values"],"warnings":result["warnings"],
            "limits":"noncrossing SCBA regulator order; external eta not a physical collision rate"}


def _complex_eb(z):
    return {"real":float(complex(z).real),"imag":float(complex(z).imag)}

def _eb_v21(L):
    spec=EBSpec(system=ExactSystem(L=L),magnetic_step=2e-4)
    analytic=liouville_slope(spec)
    fd=finite_zeeman_slope(spec)
    ward=mixed_ward(spec)
    gauge=gauge_mixed(spec)
    return {
      "protocol":"v21_finite_q_scalar_potential_exact_disorder_static_B_derivative",
      "L":L,"q":float(2*3.141592653589793/L),
      "liouville":_complex_eb(analytic),"finite_B_difference":_complex_eb(fd),
      "finite_difference_error":float(abs(analytic-fd)),
      "mixed_ward_density_abs":float(ward["density_source_ward_abs"]),
      "mixed_ward_current_abs":float(ward["current_source_ward_abs"]),
      "pure_gauge_magnetic_total_abs":float(abs(gauge["total_mixed"])),
      "limits":"finite-q finite-omega longitudinal scalar potential, not homogeneous DC PHE"}

def _eb_v22(L):
    spec=TwoColorSpec(system=ExactSystem(L=L))
    freq=frequency_channels(spec)
    independent=independent_sylvester_channels(spec)
    prepared=static_preparation(replace(spec,omega_b=0.,omega_e=.31,eta_e=.045))
    return {
      "protocol":"v22_two_color_finite_q_scalar_E_and_dynamic_onsite_Zeeman",
      "L":L,
      "B_then_E":_complex_eb(freq["b_then_e"]),
      "E_then_B":_complex_eb(freq["e_then_b"]),
      "total":_complex_eb(freq["total"]),
      "independent_sylvester_error":float(abs(freq["total"]-independent["total"])),
      "static_prepared_response":_complex_eb(prepared["total"]),
      "static_prepared_population":_complex_eb(prepared["population"]),
      "limits":"adiabatic eta E/B are regulators, not a microscopic scattering time; not homogeneous DC PHE"}

def _eb_v23(n):
    spec=UniformSpec(n=n,omega_e=.023,omega_b=.015,eta_e=.07,eta_b=.04)
    freq=source_and_channels(spec)
    direct=algebraic_channels(spec)
    poleSpec=UniformSpec(n=n)
    residue=leading_fermi_pole(poleSpec)
    eq20=eq20_direct_band_current(poleSpec)
    return {
      "protocol":"v23_periodic_homogeneous_length_gauge_two_color_EB",
      "n":n,
      "B_then_E":_complex_eb(freq["b_then_e"]),
      "E_then_B":_complex_eb(freq["e_then_b"]),
      "total":_complex_eb(freq["total"]),
      "independent_algebraic_error":float(abs(freq["total"]-direct["total"])),
      "sm20_pole":_complex_eb(residue),
      "sm20_direct":_complex_eb(eq20),
      "pole_difference":float(abs(residue-eq20)),
      "minimum_direct_gap":float(freq["gap_min"]),
      "limits":"independent periodic two-band model, dimensionless fields, no conserving collision or physical DC PHE"}


def _continuum_dirac_v24(n):
    spec=DiracSpec(nr=n,nphi=160,omega_e=.016,omega_b=.007,eta_e=.023,eta_b=.017)
    proj=projector_channels(spec)
    diag=eigen_channels(spec)
    zero=DiracSpec(nr=n,nphi=160)
    direct=eq20_direct(zero)
    analytic=angular_exact(zero)
    return {"protocol":"source_main_eq10_continuum_two_color_length_gauge_per_zeeman_energy",
            "model":"continuum tilted Dirac with UV cutoff, e=hbar=vF=1",
            "omega_e":spec.omega_e,"omega_b":spec.omega_b,
            "B_then_E":_complex_eb(proj["b_then_e"]),
            "E_then_B":_complex_eb(proj["e_then_b"]),
            "total":_complex_eb(proj["total"]),
            "projector_eigensolver_abs_diff":float(abs(proj["total"]-diag["total"])),
            "sm20_direct_per_bbar_tau":_complex_eb(direct),
            "analytic_angular_limit":float(analytic),
            "sm20_angular_abs_diff":float(abs(direct-analytic)),
            "limits":"Not a physical dc conductivity; bbar=-g muB B, distinct from Tesla and Fig2 SI."}

def _fig2_SI_v25(n):
    inp=FigureInputs()
    literal=fig2_SI_evaluate(inp,mu_mev=20.,temp_K=0.,hypothesis="literal_hbar")
    diagnostic=fig2_SI_evaluate(inp,mu_mev=20.,temp_K=0.,hypothesis="single_2pi_diagnostic")
    return {"protocol":"published_eq11_literal_SI_vs_single_2pi_diagnostic",
            "model":"original tilted Dirac closed form and conditional Hall-bar estimate",
            "literal_hbar":literal,"single_2pi_diagnostic":diagnostic,
            "ratio":literal["minus_1e4_sigma_over_e2_hbar"]/diagnostic["minus_1e4_sigma_over_e2_hbar"],
            "limits":"The 1/(2pi) branch is diagnostic, not an accepted correction; 1000ohm assumed sheet R."}

def _si_liouville_v26(n):
    cfg=SIConfig(nr=n,nphi=max(96,2*n))
    current=eq20_physical(cfg)
    sigma=current/cfg.E_Vm
    expected=eq11_SI_magnitude(cfg)
    return {"protocol":"direct_2band_SI_Liouville_vs_main_eq11",
            "model":"tilted Dirac, Hamiltonian in joules, momentum in inverse meters",
            "current_A_per_m":_complex_eb(current),
            "conductivity_S":_complex_eb(sigma),
            "literal_eq11_magnitude_S":float(expected),
            "relative_magnitude_error":float(abs(abs(sigma)/expected-1)),
            "field_T":cfg.B_T,"electric_field_V_per_m":cfg.E_Vm,
            "limits":"SI magnitude audit of adiabatic SM20; not microscopic steady-state DC PHE."}

def _source_fig2_v27(n):
    csv_path=Path(__file__).resolve().parents[1]/"results/raw/v27_fig2b_digitized.csv"
    points=load_points(csv_path)
    fits={key:fit_branch(value) for key,value in points.items()}
    plateau=float(sum(x["amplitude"] for x in fits.values())/len(fits))
    ratio,literal=ratio_to_literal(plateau)
    return {"protocol":"published_Fig2b_approximate_raster_digitization_temperature_fit",
            "model":"published figure approximate pixel read-off, NOT authors' source data",
            "fits":fits,"platform_mean":plateau,
            "literal_SI_paper_ordinate":float(literal),"ratio_to_literal":float(ratio),
            "axis_calibration_sensitivity":calibration_sensitivity(csv_path),
            "limits":"Pixel-origin ambiguity; inferred normalization remains OPEN. No publisher image redistributed."}

CASES={
    "spin-four-band":(_spin,23,"Four-band spin-MNHE density matrix versus SM Eq.29","Synthetic nondegenerate four-band model; no impurity or material Wannier operator."),
    "orbital-phe":(_phe,61,"2D Bz orbital PHE; literal SM Eq.56 versus finite-field derivative","Two-band periodic 2D model; printed Eq.56 discrepancy remains open."),
    "orbital-mnhe-2d":(_mnhe2d,61,"Strict 2D Bz orbital MNHE squared-curvature sector","Strict 2D SBC sector only; this does not test general QMQ."),
    "orbital-mnhe-3d":(_mnhe3d,17,"3D two-band positional shift versus literal SM Eq.52","Two-band 3D model; Eq.52 factor-two discrepancy remains open."),
    "scba-ra":(_scba_ra,41,"SCBA self-energy, RA bubble and noncrossing current ladder","Independent periodic regulator; RA Fermi-surface kernel only."),
    "kubo-spectral":(_spectral,31,"Symmetric four-channel spectral Kubo with RA, AR, RR, AA ladders","Zero-temperature finite-eta spectral sector; not complete Kubo-Streda."),
    "finite-qomega-ward":(_ward,31,"Finite-(q,omega) SCBA resolvent Ward identity","External eta is a regulator; real-frequency physical continuity and dc limits remain open."),
    "exact-kubo-v20":(_exact_kubo_v20,5,"Exact finite-system real-frequency continuity and contact terms","L is a finite periodic-system linear size, not a k-mesh. Exact longitudinal benchmark; not PHE."),
    "limit-order-v20":(_limit_order_v20,41,"SCBA regulator limit at fixed finite probe versus static origin","Numerical regulator-order evidence, not exact dc electromagnetic transport."),
    "eb-finite-q-v21":(_eb_v21,4,"Finite-q equilibrium-B mixed Kubo/Liouville derivative and gauge contact","Finite-q scalar-electric probe, not published homogeneous DC PHE."),
    "eb-two-color-v22":(_eb_v22,3,"Independently switched E then B and B then E density-matrix channels","Finite-q scalar-potential and onsite Zeeman, adiabatic regulators."),
    "eb-uniform-v23":(_eb_v23,64,"Uniform length-gauge EB and SM Eq.20 pole model benchmark","Independent periodic model, not original continuum surface or steady-state transport."),
    "continuum-dirac-v24":(_continuum_dirac_v24,96,"Original tilted Dirac two-color kernel, Eq20 current and angular limit","Model-level E field and Zeeman-energy protocol, no steady-state DC."),
    "fig2-si-v25":(_fig2_SI_v25,0,"Published Eq11 literal SI vs 1/(2pi) diagnostic and voltage estimate","Independent units check; sheet-resistance assumption and figure normalization OPEN."),
    "si-liouville-v26":(_si_liouville_v26,92,"Direct physical-SI band current vs Eq11 restored magnitude","Adiabatic Liouville coefficient, not measurable steady-state DC."),
    "fig2-curve-v27":(_source_fig2_v27,0,"Digitized approximate Fig2b curve, fit temperatures, amplitudes and calibration","Approximate PDF pixel data, not published raw dataset or accepted 2pi correction."),
}

def evaluate(case,mesh=None):
    if case not in CASES: raise ValueError("Unknown case; use one of "+",".join(CASES))
    fn,default_mesh,description,limits=CASES[case]
    n=default_mesh if mesh is None else mesh
    if case in ("fig2-si-v25","fig2-curve-v27"):
        if mesh is not None:raise ValueError("No mesh override for closed-form or frozen digitization")
    elif case in ("exact-kubo-v20","eb-finite-q-v21","eb-two-color-v22"):
        if not 3<=n<=12:raise ValueError("Finite exact system requires 3<=L<=12")
    elif case=="eb-uniform-v23":
        if n<12 or n%2:raise ValueError("Homogeneous FFT mesh must be even n>=12")
    elif case=="continuum-dirac-v24":
        if n<24:raise ValueError("Continuum radial quadrature requires >=24")
    elif case=="si-liouville-v26":
        if n<32:raise ValueError("Physical SI radial quadrature requires >=32")
    elif n<12: raise ValueError("Mesh must be >=12 for momentum-grid cases")
    result=fn(n)
    return {"schema":"qgmap-advanced-2","case":case,
            "science_status":"PARTIAL","validation_scope":"named synthetic model / method subset",
            "description":description,"limitations":limits,
            "source":"v30 offline independently audited implementation",
            "parameters":{"N":n},"result":result}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case",required=True,choices=sorted(CASES))
    parser.add_argument("--mesh",type=int,help="override the case's recommended diagnostic mesh")
    parser.add_argument("--output",type=Path,help="write JSON instead of stdout")
    opt=parser.parse_args()
    data=json.dumps(evaluate(opt.case,opt.mesh),indent=2,allow_nan=False)
    if opt.output is None:print(data)
    else:
        opt.output.parent.mkdir(parents=True,exist_ok=True)
        opt.output.write_text(data+"\n",encoding="utf-8")
if __name__=="__main__":main()
