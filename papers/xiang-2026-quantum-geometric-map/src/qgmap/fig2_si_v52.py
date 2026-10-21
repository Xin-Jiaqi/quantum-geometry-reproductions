"""v52 independent absolute SI audit of published PRB Eq11, Fig2b and 43 microV.
Literal Eq11: Delta_sigma/(e^2/hbar)=(tau*g*mu_B*Bx/hbar)*F(beta).
The alternative using e^2/h instead of e^2/hbar is separately reported,
NOT silently inserted as a manuscript correction or proof of author intent.
The original paper's main Eq20/31 sign and orbital terms remain OPEN.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
import math
from scipy import constants

SOURCE_MAIN_SHA256='b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730'
SOURCE_SM_SHA256='9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4'

@dataclass(frozen=True)
class Fig2Config:
    beta:float=1/8
    tau_ps:float=.5
    g_factor:float=10.
    bx_tesla:float=1.
    ey_v_per_m:float=1e4
    lateral_um:float=100.
    resistance_ohm:float=1e3
    def __post_init__(self):
        if not all(math.isfinite(x) for x in asdict(self).values()):
            raise ValueError('all SI parameters finite')
        if not 0<self.beta<1 or self.tau_ps<=0 or self.g_factor<=0 or self.bx_tesla<=0 or self.ey_v_per_m<=0 or self.lateral_um<=0 or self.resistance_ohm<=0:
            raise ValueError('0<beta<1 and positive tau,g,B,E,width,R')

def eq11_f(beta):
    if not math.isfinite(beta) or not 0<beta<1:raise ValueError('Eq11 beta')
    return -beta/(4*math.pi*(1+math.sqrt(1-beta*beta))**2)

def independent_eq11_raw(beta):
    if not math.isfinite(beta) or not 0<beta<1:raise ValueError('Eq11 beta')
    return (beta*beta+2*math.sqrt(1-beta*beta)-2)/(4*math.pi*beta**3)

def absolute_fig2_audit(cfg:Fig2Config=Fig2Config()):
    tau=cfg.tau_ps*1e-12
    width=cfg.lateral_um*1e-6
    magnetic_energy=cfg.g_factor*constants.physical_constants['Bohr magneton'][0]*cfg.bx_tesla
    magnetic_energy_ev=magnetic_energy/constants.e
    pref=tau*magnetic_energy/constants.hbar
    geometric=eq11_f(cfg.beta)
    e2hbar=constants.e**2/constants.hbar
    e2h=constants.e**2/constants.h
    normalized=pref*geometric
    literal=normalized*e2hbar
    alternative=normalized*e2h
    fig_height=1.75  # hand-read original published figure within +/- .03
    figure_tolerance=.03
    voltage_claim=43. # authors write approximately 43 microV
    vgeom=cfg.ey_v_per_m*width*cfg.resistance_ohm*1e6
    literal_height=-1e4*normalized
    alternative_height=-1e4*alternative/e2hbar
    return {
      'protocol':'v52_PRB_printed_Eq11_vs_Fig2b_and_43microV_absolute_SI_audit',
      'scientific_status':'PARTIAL',
      'source_main_SHA256':SOURCE_MAIN_SHA256,'source_sm_SHA256':SOURCE_SM_SHA256,
      'source_statement':'Fig2(b) prints -10^4 sigma^(1)_xy,x g mu_B B_x in units e^2/hbar; text claims ~43 microV',
      'params':asdict(cfg),
      'physics_constants':{
        'bohr_magneton_j_per_t':constants.physical_constants['Bohr magneton'][0],
        'bohr_magneton_ev_per_t':constants.physical_constants['Bohr magneton'][0]/constants.e,
        'planck_h_j_s':constants.h,'hbar_j_s':constants.hbar,
        'e_squared_over_h_siemens':e2h,'e_squared_over_hbar_siemens':e2hbar,
        'ratio_e2hbar_over_e2h':e2hbar/e2h},
      'geometry_F_beta':geometric,
      'zeeman_g_muB_B_ev':magnetic_energy_ev,
      'tau_zeeman_over_hbar':pref,
      'literal_Eq11':{
         'sigma_S':literal,'sigma_over_e2hbar':normalized,
         'figure_axis_minus_1e4':literal_height,
         'predicted_voltage_microV':abs(literal)*vgeom},
      'independent_e2h_variant':{
         'sigma_S':alternative,'sigma_over_e2hbar':alternative/e2hbar,
         'figure_axis_minus_1e4':alternative_height,
         'predicted_voltage_microV':abs(alternative)*vgeom},
      'published':{
         'figure_plateau_read_approx':fig_height,
         'figure_digitization_tolerance':figure_tolerance,
         'hall_voltage_claim_microV_approx':voltage_claim,
         'figure_axis_unit_printed':'e^2 / hbar',
         'fig2b_figure_temps_K':[0.,20.,40.]},
      'quantitative_comparisons':{
         'literal_plateau_divided_by_graph_readout':literal_height/fig_height,
         'literal_voltage_divided_by_quoted_43microV':abs(literal)*vgeom/voltage_claim,
         'literal_divided_by_e2h_variant':literal_height/alternative_height,
         'h_variant_plateau_minus_observed':alternative_height-fig_height,
         'h_variant_voltage_minus_quoted_microV':abs(alternative)*vgeom-voltage_claim},
      'not_proven':[
         'Whether the authors used e^2/h as their intended plotting convention',
         'How the original SM Eq20/31 signed PHE branch entered the published figure',
         'Original plotting code has not been independently inspected',
         'No arbitrary 2pi correction is accepted as an author erratum',
         'Full original physical steady-state PHE and orbital multiband remain OPEN']}

# Source data from original printed PDF page 4, 8x direct pixel crop:
# fitz clip (140,55,285,215), x(0 meV)=718.5, y(zero)=418.5, 170.5 pixels per axis unit.
# Green is T20K, blue is T40K. No OCR or fitted amplitude.
PUBLISHED_THERMAL_PIXEL_SAMPLES=(
  (20,-10,715.),(20,-5,683.),(20,5,151.),(20,10,122.),
  (40,-10,685.),(40,-5,601.),(40,5,232.),(40,10,150.))

def finite_temperature_shape_audit(cfg:Fig2Config=Fig2Config(),samples=PUBLISHED_THERMAL_PIXEL_SAMPLES):
    """Finite-T tanh is an independent *model continuation*, source Eq11 is T=0."""
    amplitude=absolute_fig2_audit(cfg)['independent_e2h_variant']['figure_axis_minus_1e4']
    out=[]
    for T,mu,pix in samples:
        if not all(math.isfinite(x) for x in (T,mu,pix)) or T<=0:
            raise ValueError('finite T>0 and source pixel required')
        observed=(418.5-pix)/170.5
        expected=amplitude*math.tanh(mu*.001*constants.e/(2*constants.Boltzmann*T))
        out.append({'T_K':T,'mu_meV':mu,'source_pixel_y':pix,
                    'figure_axis_digitized':observed,
                    'no_amplitude_fit_prediction':expected,
                    'digitized_minus_prediction':observed-expected})
    return {'protocol':'v52_published_Fig2b_temperature_curve_axis_digitization',
       'scientific_status':'PARTIAL','pixels_from_original_published_page4':True,
       'pixel_calibration':{'zero_x':718.5,'zero_y':418.5,'pixels_per_meV':8.5,'pixels_per_vertical_unit':170.5},
       'input_samples':out,'max_abs_residual':max(abs(x['digitized_minus_prediction']) for x in out),
       'limitation':'Eq11 is T=0; tanh continuation is independent model, not published formula. Pixel uncertainty approximately 0.03.'}

def full_report(cfg:Fig2Config=Fig2Config()):
    return {'absolute_SI':absolute_fig2_audit(cfg),
        'finite_T_source_figure_digitization':finite_temperature_shape_audit(cfg),
        'scientific_status':'PARTIAL'}
