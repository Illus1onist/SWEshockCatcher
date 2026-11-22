"""Compare the 1D axisymmetric and 2D radial dam-break solutions.

Reads
    radial1D/2/h, radial1D/2/hU   (1D r-weighted OpenFOAM case)
    radial2D/2/h, radial2D/2/hU   (2D Cartesian OpenFOAM case)
and writes
    radial2d_summary.json
    radial_1d_vs_2d.png
"""
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
NX = NY = 400
DX = 20.0/NX


def read_scalar(path):
    txt = Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;", txt, re.S)
    if m:
        return np.array([float(v) for v in m.group(2).split()])
    m = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", txt)
    return np.full(1, float(m.group(1)))


def read_vector(path):
    txt = Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;", txt, re.S)
    if m:
        return np.array([[float(x) for x in v.split()]
                         for v in re.findall(r"\(([^()]*)\)", m.group(2))])
    m = re.search(r"internalField\s+uniform\s+\(([^()]*)\)\s*;", txt)
    v = [float(x) for x in m.group(1).split()]
    return np.array([v])


def front(x, h):
    g = np.abs(np.gradient(h, x))
    i = int(np.argmax(g))
    return float(x[i]), float(h[i]), float(g[i])


# 1D axisymmetric
h1 = read_scalar(ROOT/"radial1D"/"2"/"h")
u1 = read_vector(ROOT/"radial1D"/"2"/"hU")[:, 0]/np.maximum(h1, 1e-12)
r1 = (np.arange(h1.size) + 0.5)*(10.0/h1.size)
r1f, h1f, _ = front(r1, h1)

# 2D Cartesian
h2 = read_scalar(ROOT/"radial2D"/"2"/"h").reshape(NY, NX)   # [j, i]
v2 = read_vector(ROOT/"radial2D"/"2"/"hU")
ux = v2[:, 0].reshape(NY, NX)/np.maximum(h2, 1e-12)
uy = v2[:, 1].reshape(NY, NX)/np.maximum(h2, 1e-12)
x = -10.0 + (np.arange(NX) + 0.5)*DX
y = -10.0 + (np.arange(NY) + 0.5)*DX
X, Y = np.meshgrid(x, y)
R = np.hypot(X, Y)
TH = np.arctan2(Y, X)
ur = ux*np.cos(TH) + uy*np.sin(TH)

h0 = read_scalar(ROOT/"radial2D"/"0"/"h").reshape(NY, NX)
mass0 = h0.sum()*DX*DX
mass = h2.sum()*DX*DX

# radial average of the 2D field
bins = np.arange(0.0, 10.0 + 0.5*DX, DX)
ib = np.digitize(R.ravel(), bins) - 1
hr = np.full(bins.size - 1, np.nan)
urAvg = np.full(bins.size - 1, np.nan)
for b in range(bins.size - 1):
    m = ib == b
    if m.any():
        hr[b] = h2.ravel()[m].mean()
        urAvg[b] = ur.ravel()[m].mean()
rc = 0.5*(bins[1:] + bins[:-1])
ok = np.isfinite(hr)
rcf, h2f, _ = front(rc[ok], hr[ok])
u2max = float(np.nanmax(np.abs(urAvg)))
h2max = float(np.nanmax(hr))
r2max = float(rc[np.nanargmax(hr)])

summary = {
    "t": 2.0,
    "1d": {
        "h_center": float(h1[0]), "h_max": float(h1.max()),
        "r_front": r1f, "h_front": h1f, "u_max": float(np.abs(u1).max()),
        "h_min": float(h1.min()),
    },
    "2d": {
        "h_center": float(h2[199:201, 199:201].mean()),
        "h_max": h2max, "r_at_h_max": r2max,
        "r_front": rcf, "h_front": h2f, "u_max": u2max,
        "h_min": float(h2.min()),
        "mass": float(mass), "mass0": float(mass0),
        "mass_rel": float(mass/mass0 - 1.0),
        "symmetry_x": float(np.abs(h2 - h2[:, ::-1]).max()),
        "symmetry_y": float(np.abs(h2 - h2[::-1, :]).max()),
        "symmetry_xy": float(np.abs(h2 - h2.T).max()),
    },
}
(ROOT/"radial2d_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

fig, ax = plt.subplots(2, 2, figsize=(14, 11))

ax[0, 0].plot(r1, h1, lw=1.6, label="1D axisymmetric (radial1D)")
ax[0, 0].plot(rc[ok], hr[ok], "--", lw=1.6, label="2D radial average (radial2D)")
ax[0, 0].set_xlim(0, 7)
ax[0, 0].set_ylim(0, 1.05)
ax[0, 0].set_xlabel("r [m]")
ax[0, 0].set_ylabel("h [m]")
ax[0, 0].set_title("Depth at t = 2 s")
ax[0, 0].grid(alpha=0.3)
ax[0, 0].legend(fontsize=8)

ax[0, 1].plot(r1, u1, lw=1.6, label="1D axisymmetric")
ax[0, 1].plot(rc[ok], urAvg[ok], "--", lw=1.6, label="2D radial average")
ax[0, 1].set_xlim(0, 7)
ax[0, 1].set_xlabel("r [m]")
ax[0, 1].set_ylabel("u_r [m/s]")
ax[0, 1].set_title("Radial velocity at t = 2 s")
ax[0, 1].grid(alpha=0.3)
ax[0, 1].legend(fontsize=8)

pcm = ax[1, 0].pcolormesh(X, Y, h2, shading="auto", cmap="viridis", vmin=0.0, vmax=0.25)
ax[1, 0].set_xlim(-7, 7)
ax[1, 0].set_ylim(-7, 7)
ax[1, 0].set_aspect("equal")
ax[1, 0].set_xlabel("x [m]")
ax[1, 0].set_ylabel("y [m]")
ax[1, 0].set_title("2D depth field at t = 2 s")
fig.colorbar(pcm, ax=ax[1, 0], shrink=0.85, label="h [m]")

# symmetry cuts
j0 = int(np.argmin(np.abs(y)))
i0 = int(np.argmin(np.abs(x)))
ax[1, 1].plot(x, h2[j0, :], lw=1.2, label="h(x, y=0.025)")
ax[1, 1].plot(y, h2[:, i0], "--", lw=1.2, label="h(x=0.025, y)")
ax[1, 1].plot(r1, h1, ":", lw=1.4, label="1D axisymmetric")
ax[1, 1].set_xlim(0, 7)
ax[1, 1].set_ylim(0, 1.05)
ax[1, 1].set_xlabel("x or y [m]")
ax[1, 1].set_ylabel("h [m]")
ax[1, 1].set_title("Axis cuts and 1D reference")
ax[1, 1].grid(alpha=0.3)
ax[1, 1].legend(fontsize=8)

fig.tight_layout()
fig.savefig(ROOT/"radial_1d_vs_2d.png", dpi=150)
print("saved", ROOT/"radial_1d_vs_2d.png")
