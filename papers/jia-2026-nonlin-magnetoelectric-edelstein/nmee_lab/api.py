"""Explicit-protocol research NMEE API; not a single universal response theory."""
from __future__ import annotations
from dataclasses import dataclass,asdict
from typing import Literal
import numpy as np
from .clean_core.core import BlochModel,integrate_chunked
from .ptbackend.cumnas import Model as CuMnAsParams, points as cumnas_mesh,PA,PB,IX,IY,IZ
from .ptbackend.conserving_bgk import response as bgk_response


@dataclass(frozen=True)
class CleanSpec:
    """Finite-dimensional, nondegenerate, k-independent spin/observable operators."""
    alpha:str='x'
    gamma:str='z'
    beta:int=0
    mu_eV:float=.025
    temperature_K:float=10.
    branch:Literal['full_static','printed']='full_static'
    chunk_size:int=1024


@dataclass(frozen=True)
class CuMnAsBGKSpec:
    """CuMnAs only: B equilibrates first, then uniform dc E drives steady state.

    Results in the model's historical *program normalization*; no conversion
    to experimental SI response. Not the PRL's unique intrinsic prescription.
    """
    theta_rad:float=.15
    mu_eV:float=.03
    eta_pop_eV:float=.013164239138
    eta_coh_eV:float=.013164239138
    tau_anisotropy:float=0.
    observable:Literal['A_z','B_z','total_x']='A_z'
    magnetic_component:Literal['x','y']='x'
    electric_beta:int=0
    B_step_bar_eV:float=1.e-6


def evaluate_clean(model:BlochModel,points:np.ndarray, spec:CleanSpec):
    """The printed and full-static branches are both returned for provenance."""
    if spec.branch not in ('full_static','printed'):
        raise ValueError('invalid response branch')
    result=integrate_chunked(model,points,alpha=spec.alpha,gamma=spec.gamma,
                  beta=spec.beta,mu=spec.mu_eV,temperature_K=spec.temperature_K,
                  chunk_size=spec.chunk_size)
    return dict(protocol='nondegenerate_clean_band_kernel',response_branch=spec.branch,
                result=float(result[spec.branch]),parts=result,settings=asdict(spec),
                normalization='(nm T V)^-1 under legacy 2D Pauli/Bloch convention',
                warnings=['full_static includes disputed Eq(9) diagonal term; printed omits it',
                          'finite-temperature driven dc/metallic response not established by this clean band-sum',
                          'Wannier zero-connection mapping is exploratory unless covariant operators verified'])


def evaluate_cumnas_bgk(spec:CuMnAsBGKSpec, *,n=101,params=None):
    """PT-adapted conserving BGK, restricted to the validated CuMnAs model."""
    p=params if params is not None else CuMnAsParams()
    if not isinstance(p,CuMnAsParams):raise TypeError('params must be CuMnAsParams')
    if spec.electric_beta not in (0,1):raise ValueError('electric_beta must be 0 or 1')
    if n<2:raise ValueError('n>=2 required')
    spin={'x':IX,'y':IY}[spec.magnetic_component]
    obs={'A_z':PA@IZ,'B_z':PB@IZ,'total_x':IX}[spec.observable]
    data=bgk_response(cumnas_mesh(n),spec.theta_rad,spec.mu_eV,p,
                      eta_p0=spec.eta_pop_eV,eta_phi=spec.eta_coh_eV,
                      alpha=spec.tau_anisotropy,step_bar=spec.B_step_bar_eV,
                      spin=spin,beta=spec.electric_beta,obs=obs)
    return dict(protocol='cumnas_PT_BGK_static_B_equilibrated_then_DC_E',
                result=float(data['gamma_bgk']), rta_result=float(data['gamma_rta']),
                correction=float(data['gamma_correction']), settings=asdict(spec),
                mesh_n=n,model=asdict(p),normalization='CuMnAs original program normalization; not calibrated SI',
                warnings=['phenomenological globally charge-conserving BGK model, not microscopic impurity theory',
                          'Dirac closings require convergence and eta-dependent analysis',
                          'do not compare this number directly to the clean or printed branch'])
