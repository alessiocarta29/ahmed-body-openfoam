#!/usr/bin/env python3
"""
Convergence history of Cd and Cl, and the statistics quoted in the README.

Reads every results/data/coefficient*.dat (one file per run segment, e.g.
coefficient_0000.dat and coefficient_0500.dat after a restart), merges them in
iteration order and averages over the last N iterations.

Usage:  python3 scripts/plot_coefficients.py [N]      (default N = 250)
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "data"
FIGURE = ROOT / "results" / "figures" / "convergence.png"

N_AVG = int(sys.argv[1]) if len(sys.argv) > 1 else 250

# Experiment (Meile et al., 2011), stilt loads subtracted
CD_EXP, CL_EXP = 0.299, 0.345
# Pressure drag of the stilts, integrated in ParaView (see README)
CD_STILTS = 0.034


def read_coefficients(path):
    """Return {column name: array} from an OpenFOAM coefficient.dat file."""
    names = None
    with open(path) as f:
        for line in f:
            if line.startswith("# Time"):
                names = [h.strip() for h in line[1:].split("\t") if h.strip()]
    values = np.loadtxt(path, comments="#")
    return {name: values[:, i] for i, name in enumerate(names)}


files = sorted(DATA.glob("coefficient*.dat"))
if not files:
    sys.exit(f"No coefficient*.dat files in {DATA}")
segments = [read_coefficients(f) for f in files]

it = np.concatenate([s["Time"] for s in segments])
coeff = {name: np.concatenate([s[name] for s in segments]) for name in ("Cd", "Cl")}

# A restart can repeat an iteration already present: keep its last occurrence
_, first_in_reversed = np.unique(it[::-1], return_index=True)
keep = len(it) - 1 - first_in_reversed
it = it[keep]
coeff = {name: y[keep] for name, y in coeff.items()}

window = it > it[-1] - N_AVG
print(f"Statistics over iterations {int(it[window][0])}-{int(it[-1])}")

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
for ax, name, exp in zip(axes, ("Cd", "Cl"), (CD_EXP, CL_EXP)):
    y = coeff[name]
    mean, std = y[window].mean(), y[window].std()
    print(f"  {name}: {mean:.3f} +/- {std:.3f}   (experiment {exp:.3f}, {100 * (mean - exp) / exp:+.0f} %)")

    ax.plot(it, y, lw=0.7, color="tab:blue", label=name)
    ax.fill_between(it[window], mean - std, mean + std, color="tab:orange", alpha=0.3,
                    label=f"mean $\\pm$ std, last {N_AVG} it.")
    ax.axhline(exp, color="k", ls="--", lw=1, label="experiment")
    if name == "Cd":
        ax.axhline(mean - CD_STILTS, color="tab:green", ls=":", lw=1.5,
                   label="mean, stilts removed")
        print(f"  Cd without stilts: {mean - CD_STILTS:.3f}   ({100 * (mean - CD_STILTS - exp) / exp:+.0f} %)")

    # The first iterations after the potential-flow start are far off scale
    low, high = np.percentile(y[it > 50], [1, 99])
    low, high = min(low, exp), max(high, exp)
    pad = 0.1 * (high - low)
    ax.set_ylim(low - pad, high + pad)
    ax.set_ylabel(name)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)

axes[1].set_xlabel("SIMPLE iteration")
fig.tight_layout()
FIGURE.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(FIGURE, dpi=150)
print(f"Figure written to {FIGURE.relative_to(ROOT)}")
