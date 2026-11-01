"""Small shared export helpers for independently recomputed paper figures."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def setup():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                        'axes.spines.top': False, 'axes.spines.right': False,
                        'legend.frameon': False, 'svg.fonttype': 'none',
                        'pdf.fonttype': 42, 'figure.dpi': 120})


def export(fig, folder, name, arrays):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / (name + '.png'), dpi=180, bbox_inches='tight')
    with plt.rc_context({'svg.fonttype': 'none'}):
        fig.savefig(folder / (name + '.svg'), bbox_inches='tight')
    np.savez_compressed(folder / (name + '.npz'), **arrays)
    plt.close(fig)


def metadata(folder, values):
    Path(folder).mkdir(parents=True, exist_ok=True)
    Path(folder, 'validation.json').write_text(
        json.dumps(values, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
