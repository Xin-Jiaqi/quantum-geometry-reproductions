# Publication and provenance boundary

- The public `stage6_harness/` holds **task contracts and synthetic unit tests only**. Never upload VASP's licensed `POTCAR`, proprietary binaries, personal cluster credentials or production structures whose redistribution rights are uncertain.
- Real DFT calculations and VASP output live in a **private** campaign root. Share a reviewed sanitized receipt, convergence summaries and SHA256 fingerprints as appropriate. An input SHA256 does not establish provenance; sources/permissions and magnetic mapping must be reviewed separately.
- Source data (published Zhu et al. AgF₂ BSF, reported domain pairing and magnetism), model evidence (earlier 8-band synthetic hopping study) and new *computed* first-principles evidence must be tagged differently.
- For a preliminary unsuccessful or blocked run, retain original failure reports and settings. No silent retry with modified U, magnetic order or k mesh to obtain desired switching parity.
- Physical response signs/normalizations must never be selected by fitting the published Xiang Fig. 2(d) plot (the absolute magnitude is an acknowledged open question).
- DFT output and independent scientific claim approval remain manual gates, not GitHub CI results.
