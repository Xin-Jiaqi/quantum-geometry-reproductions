"""Reproduce PRL 135 146303 model branches as machine-readable JSON."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from paper_models import reproduce

if __name__=='__main__':
    print(json.dumps(reproduce(),sort_keys=True,indent=2))
