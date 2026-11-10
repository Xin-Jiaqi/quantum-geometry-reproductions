# Figure provenance and QA

- Claim: a positive gate weakens both opposite local responses in A, strengthens
  both in B, while the net changes have the same sign.
- Quantitative grid, 183 x 86 mm. Python/matplotlib only; PDF/SVG editable text,
  600-dpi PNG. Visually inspected: no clipping or overlapping labels.
- 10 independent states (5 reused + 5 new), 3 response cutoffs, 2 magnetic-probe
  steps: all 60 response records retained in source_data.csv. The 24 positive
  gate differences are retained in positive_gate_changes.csv. The figure shows
  the prespecified highest cutoff/standard step, not a favorable subset.
- No fitted or synthetic points, statistical replicates, error bars, or p values.
  Lines guide the eye; net baselines are subtracted separately within each mesh.
- Geometry hashes, layer-sum closure, insulating gaps, opposite-field mirror
  relations and magnetic-step changes checked by summarize.py. No claim of
  fully converged net amplitude, complete experimental response, or safe
  nonswitching gate range.
- Source preflight warnings about TIFF omission (PNG supplied) and width parser
  (183/25.4 inches is exactly 183 mm) are benign for this research archive.
- No external illustration, private manuscript figure, or image manipulation.
