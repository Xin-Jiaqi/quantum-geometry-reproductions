# Round 9 handoff | 2026-10-20

**Scientific status:** `PARTIAL`; **engineering status:** `PASS_LOCAL_DRAFT`; **remote deployment:** NOT DONE.

Paths:

- `/mnt/data/igmc_reproduction_round9/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`
- Destination within `Xin-Jiaqi/quantum-geometry-reproductions`: `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`.
- Staged workflow outside paper module: `/mnt/data/igmc_reproduction_round9/.github/workflows/igmc.yml`.

Commands:

```bash
cd /mnt/data/igmc_reproduction_round9/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
PYTHONPATH=src python scripts/verify_round9.py
PYTHONPATH=src python -m unittest discover -s tests -v
```

Inputs: R1–R8 `parameters*.json` and all frozen numerical baselines unchanged. Round 9 contract `release/release_contract.json`; output text `PASS_LOCAL_DRAFT` and saved logs external to the paper package.

**Evidence:** 76 tests discovered (75 pass, 1 optional native-toolkit skip). Strict old freeze unchanged and all prior model/SM and 122 MPG acceptance checks rerun locally. No claim of full published figure or original-author data agreement. No author email. No GitHub commit or CI run.

**Important:** The actual public LICENSE of this paper module remains undecided; BSD-3-Clause applies to the attributed toolkit-derived magnetic generator dataset ONLY. The user should approve publication and license before an actual push/PR. Review `docs/RELEASE_CANDIDATE_ROUND9.md`.

**Next:** after approval, use GitHub repo's branch/PR interface to add exactly this new paper subtree and the new IGMC workflow, then inspect workflow logs and report exact SHA. Keep frozen Jia NMEE directory and nmee.yml unchanged. Author figure Fig.2(d) normalization and Eq.14 convention remain open.
