# Round 9: IGMC validated-method release candidate (LOCAL DRAFT)

## Meaning of this release

Independent reproduction of selected clean two-band model methods from Xiang et al., PRL **134**, 116301 (2025), DOI: `10.1103/PhysRevLett.134.116301`. This is **not** the authors' code, not a verbatim reproduction of all figures, and not an adjudication of an alleged paper error. Scientific acceptance remains `PARTIAL`.

### Accepted validated subsets

- Eq. (11) conduction and Eq. (12) displacement model kernels with stated small-frequency conditions and Pauli-spin normalization. The Rashba test, Dirac insulating response, and metallic Fermi-boundary integrations have independent analytic/numerical checks.
- Main Eq. (12): Gauss and midpoint BZ integrations, finite-field Wilson loop, direct Liouville, d-vector algebra, and wavepacket/lesser GF cross-check support the same low-frequency response in the tested models.
- Symmetry: 122 magnetic point groups; 244 full C/D tensor subspaces; 243 exact matches to printed SM table shapes. One Table II relation remains disputed.

### Unresolved / NOT accepted

1. Fig. 2(d) paper-plot absolute ordinate: model-converged insulating `+0.01775389 A/(m T)` and metallic `-0.06417845 A/(m T)` at ±1 eV, `g=2`, cyclic f=1 THz; the original plot's normalization/source data have not been independently recovered. No factor-2 rescaling permitted.
2. Main Eq. (14)'s single antiunitary factor vs channel-resolved C/D real-response symmetries: author's complex-frequency convention not resolved.
3. SM Eq. (34) intermediate factor and Eq. (46) Fourier coefficient: local source-formula issues documented but not substituted into accepted physics. SM Table II `-42'm'` equality source audit remains open.
4. Nondegenerate clean two-band results do not prove degenerate multiband, disorder-dressed response or real-material Wannier IGMC.
5. External `group-theory-operations-toolkit` native API integration is optional and was skipped in this offline runtime; numerical generators are pinned with their BSD-3-Clause notice.

## How to verify

Run from **this paper subdirectory**:

```bash
python -m pip install -e .
python scripts/verify_round9.py
```

Or from a repository checkout, follow `.github/workflows/igmc.yml`. This performs Python unit tests and frozen regression, but does **not** prove GitHub CI until published and checked on GitHub. Prefer a fresh Python 3.11 environment with dependencies resolved; isolated source re-extraction can be validated offline using pre-installed dependencies.

## Provenance / publication checklist

- Separate paper directory only; existing Jia NMEE freeze and workflow remain unchanged.
- Source code and synthetic model plots, numerical data, and own explanatory notes only; no article or SM publisher PDF, scanned table images, project proprietary archives, or user's long original notes.
- `release/release_contract.json` records exact gate, known missing data, and no empirical rescaling.
- `.github/workflows/igmc.yml` staged in the distribution root but **not yet installed remotely**.
- Public license for THIS IGMC code/module still **not selected**, so public-source sharing requires a license decision; toolkit generator data has a separate BSD-3-Clause attribution.
- Correspondence draft is not sent; author graph/source code still absent.
- The current distribution is a reviewable LOCAL DRAFT; do not call its GitHub deployment or CI PASS.
