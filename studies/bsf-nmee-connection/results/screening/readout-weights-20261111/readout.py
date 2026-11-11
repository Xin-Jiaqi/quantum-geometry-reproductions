"""Readout of a fixed weighted sum of layer spin responses, from archived SCFs."""
import json
from math import isclose
from pathlib import Path

folder = Path(__file__).resolve().parent
review = folder.parent / "review7268-20261111"
anchors = json.loads((review / "matched_kernel_interfaces.json").read_text())
pilot = json.loads((review / "prediction_review.json").read_text())
states = [row for row in anchors if row["interface"] == "S/S"]
states.extend(dict(row, interface="S/Se") for row in pilot["mesh_metrics"])

rows = []
for state in states:
    net, lower, upper = state["G"]
    assert isclose(lower + upper, net, abs_tol=1e-8)
    assert upper < 0 < lower
    rows.append({
        "interface": state["interface"],
        "mesh": state["mesh"],
        "bands": state["bands"],
        "G_lower": lower,
        "G_upper": upper,
        "G_equal_weights": net,
        "r_zero": -lower / upper,
        "G_at_r_0.25": lower + 0.25 * upper,
    })

lower_bound = max(row["r_zero"] for row in rows if row["interface"] == "S/S")
upper_bound = min(row["r_zero"] for row in rows if row["interface"] == "S/Se")
assert lower_bound < 0.25 < upper_bound
result = {
    "sources": ["review7268-20261111/matched_kernel_interfaces.json",
                "review7268-20261111/prediction_review.json"],
    "definition": "M = C_lower + r*C_upper; r = w_upper/w_lower, positive constant weights, w_lower normalized to 1",
    "G_units": "1e-9 (V^2 T)^-1",
    "opposite_slope_open_interval_on_reported_settings": [lower_bound, upper_bound],
    "rows": rows,
}
(folder / "weights.json").write_text(json.dumps(result, indent=2) + "\n")
