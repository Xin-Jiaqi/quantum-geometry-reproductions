# Stage 8d: official BiDB download handoff and WSe2 provenance gate

**Status: BLOCKED_ORIGINAL_RELAXED_STRUCTURE_NOT_RECEIVED.** No ab-initio relaxation, Wannier calculation, or real-material static/driven-DC NMEE performed. This is *independent* BSF × NMEE follow-up, not a published Jia et al. PRL result or change to its frozen reproduction.

## Where to obtain the actual raw data

The high-priority target is the **original 3R-WSe2 AB/BA six-atom relaxed bilayer** geometry in CIF/POSCAR or a first-party database containing it.

1. **BiDB original database docs**, https://cmr.fysik.dtu.dk/bidb/bidb.html . Web-indexed content says `Download data: bidb.db`. Current automated retrieval of this **legacy page returns 404**; it may also fail in the user's browser. Do not invent a working direct `.db` URL.
2. **BiDB archival DOI explicitly cited in the publication's Data Availability**, https://doi.org/10.11583/DTU.24849819 . Try this persistent front door if the older docs page fails; select original file with geometry. This DOI/landing page could not be retrieved automatically here either. Citation: Pakdel et al., *Nature Communications* 15, 932 (2024), DOI 10.1038/s41467-024-45003-w.
3. **BMDS alternate author dataset**, https://doi.org/10.6084/m9.figshare.21799475 . Barik and Woods, *Scientific Data* 10, 232 (2023), DOI 10.1038/s41597-023-02146-7; 760 relaxed bilayers. A DOS spreadsheet or generator code alone does not certify the target AB atomic motif.
4. The already-located [MoS2 AB-3R Wannier90/VASP archive](https://zenodo.org/records/21480888) explicitly names `AB-3R/POSCAR` but is **restricted**. This is not a public direct-download alternative unless the owner grants legitimate access.

**Requested upload:** unchanged original `bidb.db`, or the primary-source `WSe2` 3R bilayer `POSCAR`/`CIF`/`zip` along with its source page URL. If the site cannot be opened, a screenshot of the error is useful. Avoid downloading terabytes or a large irrelevant full dataset when a small author-coordinate file is provided.

## Critical methodological caveat from the actual BiDB paper

The extensive screening computes binding energies via a **z scan of the frozen monolayer atomic geometries**, supplemented by a **separate full PBE-D3 relaxation** of some bilayer structures. These are not synonyms; an extracted `W2Se4` geometry or recorded `distance` does **not** imply it was fully relaxed. For the downloaded target, retain its `bilayer_uid`, structural provenance, and the relaxation-method metadata. Distinguish also bulk 3R periodic structures from isolated bilayer slabs.

## Reproducible local extraction, fail closed

The separate study module [extract_bidb_targets.py](extract_bidb_targets.py) requires `pip install ase` **only for actual extraction**, plus an original ASE-SQLite database file:

```bash
python -m pip install ase
python extract_bidb_targets.py /path/to/bidb.db --out extracted_targets
```

It exports candidates with precise atom counts W2Se4/Mo2S4/W2S4, original database SHA-256 and ASE row IDs; every geometry remains explicitly labelled `UNVERIFIED_AB_BA_AND_RELAXATION`. It **never** assigns AA/AB/BA registry from the mere formula or proves a DFT relaxation. Unit tests cover exact 6-site formulas, reject unsupported database shapes, and verify checksum generation; they cannot substitute for checking actual original author bytes.

## Physics gate after raw data arrives

For each candidate confirm exact source, atomic species, two actual MX2 layers, lattice/vacuum, W-on-Se/Se-on-W stacking and author relaxation protocol. Then feed coordinates to the pinned group-theory toolkit for the *actual atom motif*. Stage 8b's ideal D3h reference and C3v AB/BA predictions are an analytical comparator, **not** preapproved labels for a potentially symmetry-lowered relaxed structure.

Before Wannier, inspect full-BZ gap and Gamma/Lambda/K spin/orbital contributions; stage-7 synthetic K/Kprime focus is not automatically transferable to 3R-WSe2. Spinor H, Sxyz and basis overlaps/covariant velocities are still required. Preserve the unresolved difference between insulating equilibrium/static geometry and driven-DC NMEE.

**Reproduction segregation:** only add files under `studies/bsf-nmee-connection/material-screening/` and update `studies/bsf-nmee-connection/RESEARCH_INDEX.md` plus the navigation-only `papers/.../research/README.md`. Do **not** alter the frozen original paper module or historic checksummed stage README.
