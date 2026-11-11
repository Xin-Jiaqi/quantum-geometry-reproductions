#!/usr/bin/env python3
"""Stage 8k: ZERO-job nondegenerate-spin evidence ranking from EXISTING JSON.

This ranks falsifiable read-only WAVECAR-validation targets. It never claims
that pseudized spin proxies (PAW incomplete, parser dubious) are material
observables, and never calculates NMEE, Sx_cv or velocity vertices.
"""
from pathlib import Path
import json,math

HERE=Path(__file__).resolve().parent
JSON=HERE/"STAGE8I_SPIN_TEXTURE_20261022.json"
MIN_GAP_MEV=5.0
MIN_PERP=0.025

def normalized_label(original,xyz):
    u,v=xyz[:2]
    if abs(u)<1e-6 and abs(v)<1e-6:return "Gamma_TRIM__DEGENERACY_TRAP"
    if abs(u-1/3)<1e-5 and abs(v-1/3)<1e-5:return "TRUE_K"
    if abs(u+1/3)<1e-5 and abs(v+1/3)<1e-5:return "TRUE_K_PRIME"
    if abs(u+1/3)<1e-5 and abs(v-1/3)<1e-5:return "GAMMA_M_TWO_THIRDS__NOT_K_PRIME"
    if abs(v)<1e-6 and .4<abs(u)<.51:return "GAMMA_M_NEAR_M__NOT_EXACT_M"
    if abs(u-v)<1e-6 and .05<abs(u)<.30:return "GAMMA_K_INTERIOR__NOT_TRUE_BAND_EDGE"
    return "OTHER_CHECK_RECIPROCAL_COORDINATES"

def closest_visible_band_gap(bands,i):
    # This is the sampled LOCAL band isolation only within the archived subset.
    band=bands[i]
    neighbors=[abs(float(band["energy_eV"])-float(b["energy_eV"]))
      for b in bands if abs(int(b["band"])-int(band["band"]))==1]
    return min(neighbors)*1000 if neighbors else None

def classify(json_path=JSON,min_gap=MIN_GAP_MEV,min_spin=MIN_PERP):
    d=json.loads(Path(json_path).read_text())
    results=[]
    for job,source in d["jobs"].items():
        for old_name,point in source["per_k"].items():
            label=normalized_label(old_name,point["k"])
            bands=point["bands_near_gap"]
            for i,b in enumerate(bands):
                gap=closest_visible_band_gap(bands,i)
                perp=math.hypot(float(b["Sx"]),float(b["Sy"]))
                is_deg=(gap is None or gap<min_gap)
                item={"job":job,"old_label":old_name,"correct_location":label,
                      "fractional_k":[float(x) for x in point["k"]],
                      "band":int(b["band"]),"E_eV":float(b["energy_eV"]),
                      "S_perp_pseudized":round(perp,6),
                      "Sz_pseudized":float(b["Sz"]),
                      "spinor_norm_reported":float(b.get("norm",float("nan"))),
                      "min_adjacent_band_separation_meV":round(gap,5) if gap is not None else None,
                      "nondegenerate_within_visible_immediate_neighbors":not is_deg,
                      "candidate_from_archived_proxy":bool(not is_deg and perp>=min_spin),
                      "verified_physical_spin_matrix":False,
                      "verified_Sz_nonconservation":False}
                results.append(item)
    candidates=sorted(
        (x for x in results if x["candidate_from_archived_proxy"]),
        key=lambda x:(0 if x["job"]=="soc_AB_s0.00" else 1,
                      -x["min_adjacent_band_separation_meV"],
                      -x["S_perp_pseudized"]))
    discriminators=[]
    paired={}
    for item in results:
        paired[(item["correct_location"],item["band"],item["job"])]=item
    for name in ("GAMMA_M_NEAR_M__NOT_EXACT_M",
                 "GAMMA_K_INTERIOR__NOT_TRUE_BAND_EDGE",
                 "GAMMA_M_TWO_THIRDS__NOT_K_PRIME"):
        for band in (35,36):
            a=paired.get((name,band,"soc_AB_s0.00"))
            m=paired.get((name,band,"soc_GLIDE_s0.50"))
            if a and m:
                discriminators.append({
                  "location":name,"band_index":band,
                  "AB_perp_spin_proxy":a["S_perp_pseudized"],
                  "midpoint_perp_spin_proxy":m["S_perp_pseudized"],
                  "AB_min_adjacent_gap_meV":a["min_adjacent_band_separation_meV"],
                  "midpoint_min_adjacent_gap_meV":m["min_adjacent_band_separation_meV"],
                  "cross_stack_same_band_index_is_not_proven_same_orbital":True})
    return {
      "stage":"8k","evidence":"archived pseudized WAVECAR-derived spin JSON, NOT raw-wavefunction independent reparse",
      "source_file":str(Path(json_path).name),"max_additional_PBS_jobs":0,
      "primary_question":"Can actual nondegenerate spinors show robust in-plane spin, implying [H,Sz]!=0?",
      "prior_context":"Gamma Kramers pair is gauge ambiguous; false-Kprime is interior Gamma-M; same-Sz K bands 35/36 NOT opposite-spin splitting",
      "min_band_isolation_meV":min_gap,"min_proxy_transverse_spin":min_spin,
      "positive_proxy_candidates":candidates,
      "matched_k_band_index_comparison":discriminators,
      "sampled_states_count":len(results),
      "spin_mixing_gate":"PROXY_ONLY__PAW_AND_WAVECAR_PARSER_NOT_VERIFIED",
      "actual_spin_velocity_Kubo_elements":"NOT_COMPUTED",
      "static_material_NMEE":"NOT_COMPUTED",
      "required_remote_readonly_steps":[
        "Verify EIGENVAL eigenvalue arrays and WAVECAR wavefunction records jointly at identical k and band indexes; avoid pymatgen band_energy if reader disagrees.",
        "Check spinor coefficient lengths, G-vector ordering for two spin components, normalization and time reversal; independently compute pseudized Sx/Sy/Sz from BOTH spin components in a common plane-wave basis.",
        "Use original on-site PAW projection/augmentation with license and pseudopotential data staying on HPC; distinguish pseudized proxy from full physical spin.",
        "At nondegenerate AB Gamma-M interior and near-M plus Gamma-K interior, test transverse spin against band separation, change of phase, and true K/−K spin reversal.",
        "At glide midpoint treat false-Kprime/Gamma-M 2/3 degenerate bands as GAUGE AMBIGUOUS, not a robust counterexample. Compare AB and midpoint on identical k ONLY after band/orbital character tracking.",
        "Produce sanitized numerical QA and scientific judgment. No qsub, no new DFT, no Wannier, no raw/private file uploads."
      ],
      "global_qualification":"Nondegeneracy is screened only among immediate visible bands; pseudized spinor wavefunctions and full PAW operator have not been authenticated."
    }

if __name__=="__main__":
    print(json.dumps(classify(),indent=2,ensure_ascii=False))
