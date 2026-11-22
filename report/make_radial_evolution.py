"""Time evolution of the 1D axisymmetric radial dam break."""
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
CASE = ROOT.parent/"SWEshockFoam"/"radial1D"


def read_scalar(path):
    txt = Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;", txt, re.S)
    return np.array([float(v) for v in m.group(2).split()])


def read_vector_x(path):
    txt = Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;", txt, re.S)
    return np.array([float(v.split()[0]) for v in re.findall(r"\(([^()]*)\)", m.group(2))])


fig, ax = plt.subplots(1, 2, figsize=(12.6, 4.4))
for t, name, c in ((0.5, "0.5", "#4c78a8"), (1.0, "1", "#f58518"), (2.0, "2", "#54a24b")):
    h = read_scalar(CASE/name/"h")
    u = read_vector_x(CASE/name/"hU")/np.maximum(h, 1e-12)
    r = (np.arange(h.size) + 0.5)*(10.0/h.size)
    ax[0].plot(r, h, color=c, lw=1.6, label=f"t={t:g} s")
    ax[1].plot(r, u, color=c, lw=1.6, label=f"t={t:g} s")

ax[0].set_xlim(0, 7)
ax[0].set_ylim(0, 1.05)
ax[0].set_xlabel(r"$r$ [m]")
ax[0].set_ylabel(r"$h$ [m]")
ax[0].set_title("Depth")
ax[0].grid(alpha=0.3)
ax[0].legend(fontsize=9)

ax[1].set_xlim(0, 7)
ax[1].set_xlabel(r"$r$ [m]")
ax[1].set_ylabel(r"$u_r$ [m/s]")
ax[1].set_title("Radial velocity")
ax[1].grid(alpha=0.3)
ax[1].legend(fontsize=9)

fig.tight_layout()
fig.savefig(ROOT/"figures"/"radial_evolution.png", dpi=170)
print("saved", ROOT/"figures"/"radial_evolution.png")
