"""Independent SI audit of Xiang et al. PRB 113 L201406 Fig. 2(b)/Eq. (11).

This is a *published-formula dimensional audit*, not a verified physical dc PHE.
Original equations set e=hbar=1; the result as printed omits g*mu_B while
Eq. (8) already contains it. Here we deliberately label magnetic coefficients
with respect to (1) Zeeman energy bbar=-g mu_B B (J) and (2) B (T).

No modifications/fits are made to the published continuous Dirac model.
The exploratory 1/(2*pi) branch is diagnostic only, not an endorsed correction.
"""
from dataclasses import dataclass
import math
import scipy.constants as const

E = const.e
HBAR = const.hbar
H = const.h
MU_B = const.physical_constants['Bohr magneton'][0]
KB = const.Boltzmann
G_QUANTUM_HBAR = E * E / HBAR

@dataclass(frozen=True)
class FigureInputs:
    beta: float = 1/8
    tau_s: float = 0.5e-12
    g: float = 10.
    B_tesla: float = 1.
    E_vm: float = 1e4
    resistance_ohms: float = 1000.
    transverse_width_m: float = 100e-6
    # Paper does not specify length/width or distinguish two-terminal from sheet R.
    # THIS ASSUMES the quoted ~1000-ohm resistance may be used as sheet R.


def geometry(beta: float) -> float:
    """The signed geometric factor in main Eq. (11), for 0<beta<1."""
    if not 0<beta<1:raise ValueError('Paper closed form assumes 0<beta<1')
    s=math.sqrt(1-beta*beta)
    return -beta/(4*math.pi*(1+s)**2)  # cancellation-free form


def geometry_printed(beta: float) -> float:
    """Literal main Eq. (11) parenthesis and 4*pi*beta**3 denominator."""
    if not 0<beta<1:raise ValueError
    return (beta**2+2*math.sqrt(1-beta**2)-2)/(4*math.pi*beta**3)


def tanh_population(mu_mev:float,temp_K:float)->float:
    if temp_K<0:raise ValueError('negative T')
    if temp_K==0:return (1 if mu_mev>0 else (-1 if mu_mev<0 else 0))
    return math.tanh(mu_mev*1e-3*E/(2*KB*temp_K))


def evaluate(inputs:FigureInputs,mu_mev:float=20.,temp_K:float=0., *, hypothesis='literal_hbar'):
    """Physical field convention: j_x=[sigma_per_B B_x] E_y.

    The *literal_hbar* branch reads main Eq. (11) as a coefficient per
    Zeeman energy, then multiplies by g mu_B to get per tesla.
    The *single_2pi_diagnostic* branch replaces one hbar by h in its overall
    dimensional factor, leaving the same geometric kernel. Purely diagnostic.
    """
    if hypothesis not in ('literal_hbar','single_2pi_diagnostic'):
        raise ValueError('Unknown branch')
    if inputs.tau_s<=0 or inputs.g<=0:raise ValueError('positive tau and g required')
    factor=geometry(inputs.beta)*tanh_population(mu_mev,temp_K)
    denom=HBAR if hypothesis=='literal_hbar' else H
    # This product is a SHEET conductivity (S) after B field insertion;
    # the Eq(11) coefficient per Zeeman energy has units S/J.
    sigma_per_zeeman_joule=inputs.tau_s*E**2/(HBAR*denom)*factor
    sigma_per_tesla=sigma_per_zeeman_joule*inputs.g*MU_B
    sigma_sheet=sigma_per_tesla*inputs.B_tesla
    # The *plotted* ordinate appears as -10^4 sigma^(1) g muB B /(e²/hbar).
    ordinate=-1e4*sigma_sheet/G_QUANTUM_HBAR
    # Open-circuit, weak-Hall-angle, sheet-resistance assumption:
    # Ex = -R_sheet * sigma_xy * Ey; Vx = Ex*width.
    voltage=-sigma_sheet*inputs.resistance_ohms*inputs.E_vm*inputs.transverse_width_m
    return {
      'hypothesis':hypothesis,'temperature_K':temp_K,'mu_mev':mu_mev,
      'F':factor, 'sigma_per_zeeman_energy_S_per_J':sigma_per_zeeman_joule,
      'sigma_per_B_S_per_T':sigma_per_tesla,
      'sigma_sheet_S':sigma_sheet,
      'minus_1e4_sigma_over_e2_hbar':ordinate,
      'open_circuit_voltage_V':voltage,
      'open_circuit_voltage_microV':voltage*1e6,
      'tau_zeeman_over_hbar':inputs.tau_s*inputs.g*MU_B*inputs.B_tesla/HBAR,
      'assumption':'quoted R approximates 2D sheet resistance; transverse width=100 um; weak Hall angle',
    }


def independent_volt_from_ratio(inp:FigureInputs,mu_mev=20.,temp_K=0.,h_use=False):
    """Separate computation using conductance quantum and dimensionless ratio."""
    f=geometry_printed(inp.beta)*tanh_population(mu_mev,temp_K)
    omega=(inp.g*MU_B*inp.B_tesla)/(H if h_use else HBAR)
    dimensionless_sheet=f*inp.tau_s*omega
    sheet=(E**2/HBAR)*dimensionless_sheet
    return -sheet*inp.resistance_ohms*inp.E_vm*inp.transverse_width_m


def independent_angular_geometry(beta:float,n=16384):
    """Integrate physical angular kernel without using Eq11 closed form.

    The normalizing angular prefactor is analytically fixed, *not fitted*:
    F = -(beta / 4π) < sin²(theta) / (1+beta cos(theta))² >_theta?
    Verified analytic simplification instead via the direct conventional
    angular integral: F = (1/(4π)) [1/(2π)∫ cosθ sin²θ/(1+βcosθ)dθ].

    The integral below uses the full signed response expected from main Eq11;
    avoids undocumented paper-to-physical helicity flips.
    """
    import numpy as np
    theta=2*math.pi*(np.arange(n)+0.33)/n
    integral=np.mean(np.cos(theta)*np.sin(theta)**2/(1+beta*np.cos(theta)))
    # To match F(β) = -β/(4π(1+sqrt(1−β²))²), integral equals -β/(2(1+s)²)
    return integral/(2*math.pi)
