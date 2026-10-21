"""Stage8j: non-destructive independently re-evaluate DeepSeek Stage8i JSON interpretations.

Only Python standard library. Never opens remote VASP output or submits jobs.
Emits correction metadata; it is not a replacement for raw WAVECAR validation.
"""
import json, math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
ORIGINAL=ROOT
A_A=3.319006624468

def reduced(k):
    return tuple(round((float(x)+0.5)%1-0.5,10) for x in k[:2])

def c3(k):
    u,v=k[:2]
    return reduced((-(u+v),u))

def k_length(k,a=A_A):
    # Real basis: a1=(a,0), a2=(-a/2,sqrt(3)*a/2).
    # Dual basis: b1=(2pi/a, 2pi/(sqrt(3)*a));
    #             b2=(0,4pi/(sqrt(3)*a)).
    u,v=k[:2]
    return 2*math.pi/a*math.hypot(u,(u+2*v)/math.sqrt(3))

def spin_conservation_counterexample():
    # H=epsilon*I conserves Sz. The Kramers-degenerate +x state
    # (up+down)/sqrt2 nonetheless has <Sx>=+1/2.
    return {"Sz_commutator_exactly_zero":True,
            "degenerate_state_Sx":0.5,
            "conclusion":"Gamma Kramers expectation does NOT prove Sz nonconservation"}

def evaluate():
    valley=json.loads((ORIGINAL/"STAGE8I_ISING_SPLIT_20261022.json").read_text())
    spin=json.loads((ORIGINAL/"STAGE8I_SPIN_TEXTURE_20261022.json").read_text())
    k=tuple(valley["jobs"]["soc_AB_s0.00"]["per_k"]["K"]["k"][:2])
    kp_wrong=tuple(valley["jobs"]["soc_AB_s0.00"]["per_k"]["Kp"]["k"][:2])
    true_kp=tuple(-x for x in k)
    assert reduced(c3(k))==reduced(k)
    assert reduced(kp_wrong)!=reduced(true_kp)
    # It is an interior Gamma--M point with M=(-1/2,+1/2).
    assert all(abs(x - 2*m/3)<1e-6 for x,m in zip(kp_wrong,(-.5,.5)))
    results={}
    for job,raw in valley["jobs"].items():
        states={d["band"]:d for d in spin["jobs"][job]["per_k"]["K"]["bands_near_gap"]}
        a,b=states[35],states[36]
        assert abs((b["energy_eV"]-a["energy_eV"])*1000-
                   (raw["per_k"]["K"]["bands_31_40_eV"][5]-
                    raw["per_k"]["K"]["bands_31_40_eV"][4])*1000)<1e-8
        if a["Sz"]*b["Sz"]<=0:
            raise RuntimeError("Source K band 35/36 pair no longer same Sz")
        results[job]={
            "K_true_fractional":list(k),
            "Kprime_true_fractional":list(true_kp),
            "misnamed_Kp_actual_fractional":list(kp_wrong),
            "misnamed_Kp_actual_location":"Gamma-to-M, two-thirds toward M",
            "reported_K_to_wrongKp_band36_delta_eV":
                 raw["per_k"]["K"]["VBM_band36_eV"]-raw["per_k"]["Kp"]["VBM_band36_eV"],
            "this_is_true_valley_split":False,
            "band35_band36_energy_separation_meV":
                 round((b["energy_eV"]-a["energy_eV"])*1000,5),
            "band35_Sz":a["Sz"],
            "band36_Sz":b["Sz"],
            "these_bands_form_verified_opposite_spin_pair":False,
            "original_report_spin_split_label":"NOT_SUPPORTED",
            "sampled_indirect_gap_eV":raw["global"]["sampled_indirect_gap_eV"]
        }
    return {
        "stage":"8j","source_commit":"6ffef840ca5ca5045e9073dbf781c20164b1a071",
        "scope":"read-only re-analysis of existing sanitized DeepSeek JSON, no external VASP data read",
        "reciprocal_C3_map":"(u,v)->(-u-v,u) modulo integers",
        "genuine_K":list(k),"genuine_time_reversed_Kprime":list(true_kp),
        "original_false_Kprime":list(kp_wrong),
        "K_norm_inverse_A":k_length(k),
        "false_Kp_norm_inverse_A":k_length(kp_wrong),
        "false_Kp_is_at_fraction_of_Gamma_M":2/3,
        "eigenvalue_TR_equality_rechecked_directly_against_original_EIGENVAL":False,
        "spin_conservation_counterexample":spin_conservation_counterexample(),
        "jobs":results,
        "scientific_status":"SOC_SAMPLED_GAP_PASS__SPIN_MIXING_NOT_YET_ESTABLISHED__OUTPUT_COMPLETION_WARNING",
        "new_PBS_jobs":0,"actual_material_Gamma_computed":False}

if __name__=="__main__":
    print(json.dumps(evaluate(),indent=2))
