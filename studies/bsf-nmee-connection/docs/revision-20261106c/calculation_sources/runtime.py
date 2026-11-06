"""Locate the repository for the preserved calculation sources."""
from pathlib import Path
import os

HERE = Path(__file__).resolve()
if os.environ.get("NMEE_REPO_ROOT"):
    REPO = Path(os.environ["NMEE_REPO_ROOT"]).expanduser().resolve()
else:
    REPO = next((p for p in HERE.parents if (p / "studies/bsf-nmee-connection/bsf_nmee").is_dir()), None)
if REPO is None or not (REPO / "studies/bsf-nmee-connection/bsf_nmee").is_dir():
    raise RuntimeError("Run within the repository or set NMEE_REPO_ROOT to its root")
STUDY = REPO / "studies/bsf-nmee-connection"
OUTPUT = Path(os.environ.get("NMEE_AUDIT_OUTPUT", str(HERE.parent.parent))).expanduser().resolve()
OUTPUT.mkdir(parents=True, exist_ok=True)
