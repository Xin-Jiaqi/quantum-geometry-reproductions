#!/usr/bin/env python3
"""Stage8l: zero-job scientific integrity check of TWO archived SOC JSON ledgers.

The raw WAVECAR is NOT available inside GitHub; this does not replace onsite
byte-level analysis. Find and explicitly report off-by-one band index in
Stage8k key_results while preserving the unedited source artifacts.
"""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
NEW=ROOT/"STAGE8K_INDEPENDENT_SPIN_EVIDENCE_20261022.json"
OLD=ROOT/"STAGE8I_SPIN_TEXTURE_20261022.json"
LOCATIONS={"GAMMA_M_2_3__OLD_FALSE_Kp":"Kp",
           "GAMMA_M_NEAR_M":"M","GAMMA_K_INTERIOR":"Lambda","TRUE_K":"K"}

def verify(old_path=OLD,new_path=NEW):
    old=json.loads(Path(old_path).read_text())
    new=json.loads(Path(new_path).read_text())
    issues=[]
    rows=[]
    for job,locations in new["key_results"].items():
        for place,saved in locations.items():
            previous=old["jobs"][job]["per_k"][LOCATIONS[place]]["bands_near_gap"]
            old35=next(x for x in previous if x["band"]==35)
            old36=next(x for x in previous if x["band"]==36)
            saved_energy=float(saved["E_band35_eV"])
            saved_spin=saved["S_grouped"]
            distance35=abs(saved_energy-float(old35["energy_eV"]))
            distance36=abs(saved_energy-float(old36["energy_eV"]))
            spin36_max_diff=max(abs(float(a)-float(b))
                for a,b in zip(saved_spin,(old36["Sx"],old36["Sy"],old36["Sz"])))
            spin35_max_diff=max(abs(float(a)-float(b))
                for a,b in zip(saved_spin,(old35["Sx"],old35["Sy"],old35["Sz"])))
            # At degenerate energies only spin/orbital matching resolves identity.
            unique=abs(float(old35["energy_eV"])-float(old36["energy_eV"]))>5e-6
            if unique and not (distance36<=1.0e-5 and spin36_max_diff<=1.5e-4):
                raise RuntimeError(f"Unexpected mismatch for {job} {place}: "
                                   f"E36 diff={distance36}, S36 diff={spin36_max_diff}")
            label_ok=distance35<=1e-5 and spin35_max_diff<=1.5e-4
            if unique and not label_ok:
                issues.append({"job":job,"place":place,
                    "issue":"E_band35_eV and S_grouped store BAND36 data",
                    "saved_E_eV":saved_energy,
                    "real_band35_E_eV":old35["energy_eV"],
                    "real_band36_E_eV":old36["energy_eV"],
                    "independent_saved_S_grouped_corresponds_to_band":36,
                    "saved_S_perp":round((saved_spin[0]**2+saved_spin[1]**2)**.5,6)})
            rows.append({"job":job,"location":place,
                "legacy_E_band35_eV":old35["energy_eV"],
                "legacy_E_band36_eV":old36["energy_eV"],
                "new_archived_field_E_band35_eV":saved_energy,
                "new_field_matches_actual_band36":unique and distance36<=1e-5,
                "degenerate_energy_only_cannot_disambiguate":not unique,
                "new_spin_matches_legacy_band36":spin36_max_diff<=1.5e-4})
    return {
      "stage":"8l","scientific_status":"OFF_BY_ONE_EVIDENCE_LABEL_ERROR_CONFIRMED__RAW_RECHECK_TARGETED",
      "new_vasp_jobs":0,"raw_wavecar_examined_by_this_check":False,
      "archived_artifacts_unchanged":True,
      "unexpected_band35_label_issues":issues,
      "all_records":rows,
      "physics":"The new independently saved spinor evidence confirms nonzero pseudo transverse spin in actual BAND36 at the AB three nondegenerate candidates, whereas its prose review describes BAND35. Both candidates may be nondegenerate; reconcile against onsite raw WAVECAR using explicit zero-based indices before quoting band-resolved spin.",
      "PAW_physical_spin":"NOT_QUANTITATIVELY_VALIDATED",
      "Kubo_static_NMEE":"NOT_COMPUTED"}

if __name__=="__main__":
    print(json.dumps(verify(),indent=2,ensure_ascii=False))
