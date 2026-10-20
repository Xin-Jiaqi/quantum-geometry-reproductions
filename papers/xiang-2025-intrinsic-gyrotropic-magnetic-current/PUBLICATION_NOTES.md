# GitHub release scope (validated-method subset)

This is an independently reproduced subset of Xiang et al., *Phys. Rev. Lett.* **134**, 116301 (2025); it is not the authors’ original implementation. The scientific status remains **PARTIAL**. In particular, the published Fig. 2(d) absolute ordinate, the Eq. (14) antiunitary convention, and a Table II tensor relation remain unresolved. The code MUST NOT add a numerical factor to match the original figure.

This published module is separate from the frozen Jia et al. NMEE module. The files under `archive/round7_integrity_files/` preserve internal reproducibility hashes of this specific paper’s earlier public scientific text; they do not include the publisher PDFs, scanned tables, unpublished material data, or user-provided private ZIPs. All generated figures are independent computational outputs.

The `release/round9_receipt.json` and `release/release_contract.json` are immutable historical **local draft** snapshots, included to make the Round 9 verifier reproducible. Their `LOCAL_DRAFT` fields describe that frozen checkpoint, not the later GitHub publication state. `python scripts/verify_round9.py` tests that historical physical baseline; it is not a GitHub deployment check.

## Workflow integration

A signed W90Flow top-K package contains a model Hamiltonian, not necessarily the same-gauge spin and full position/Berry operators required for a physical IGMC. Real-material IGMC is therefore **BLOCKED_OPERATOR_INPUT**, regardless of whether a package’s SHA256 checks pass. The existing NMEE Stage-1 handoff demonstrates a read-only contract, but does not establish full IGMC support.

## License

The independently developed Python implementation is released under [BSD-3-Clause](LICENSE). Reused magnetic-point-group generator data from `group-theory-operations-toolkit` retain their own [BSD-3-Clause attribution](data/TOOLKIT_BSD3_LICENSE.txt). Original article text and publisher figures are not relicensed here.
