"""Analyse all SWEshockFoam OpenFOAM cases."""
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
G = 9.81


def read_scalar(path):
    txt = Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;", txt, re.S)
    return np.array([float(v) for v in m.group(2).split()])


def read_vector_x(path):
    txt = Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;", txt, re.S)
    return np.array([float(v.split()[0]) for v in re.findall(r"\(([^()]*)\)", m.group(2))])


def stoker(hL, hR, g=G):
    cL = np.sqrt(g*hL)
    def f(h2):
        c2 = np.sqrt(g*h2)
        return 2.0*(cL - c2) - (h2 - hR)*np.sqrt(0.5*g*(h2 + hR)/(h2*hR))
    a, b = hR, hL
    for _ in range(200):
        m = 0.5*(a + b)
        if f(a)*f(m) <= 0: b = m
        else: a = m
    h2 = 0.5*(a + b)
    c2 = np.sqrt(g*h2)
    u2 = 2.0*(cL - c2)
    return h2, u2, h2*u2/(h2 - hR)


def stoker_profile(x, t, hL, hR, g=G):
    h2, u2, cs = stoker(hL, hR, g)
    cL = np.sqrt(g*hL)
    c2 = np.sqrt(g*h2)
    h = np.full_like(x, hR)
    u = np.zeros_like(x)
    h[x < -cL*t] = hL
    fan = (x >= -cL*t) & (x <= (u2 - c2)*t)
    c = (2.0*cL - x[fan]/t)/3.0
    h[fan] = c*c/g
    u[fan] = (2.0/3.0)*(cL + x[fan]/t)
    pl = (x > (u2 - c2)*t) & (x < cs*t)
    h[pl] = h2
    u[pl] = u2
    return h, u, h2, u2, cs


def ritter(x, t, hL, g=G):
    c0 = np.sqrt(g*hL)
    h = np.where(x < -c0*t, hL, np.where(x < 2.0*c0*t, (2.0*c0 - x/t)**2/(9.0*g), 0.0))
    return h, c0


def cross(x, f, level):
    for i in range(len(x)-1, 0, -1):
        if f[i] < level <= f[i-1]:
            return x[i-1] + (level - f[i-1])*(x[i] - x[i-1])/(f[i] - f[i-1])
    return None


def grid(xmin, xmax, n):
    dx = (xmax - xmin)/n
    return np.linspace(xmin + 0.5*dx, xmax - 0.5*dx, n), dx


summary = {}

# wet
n = 2000
x, dx = grid(-10, 10, n)
h = read_scalar(ROOT / "pipe1D/1/h")
u = read_vector_x(ROOT / "pipe1D/1/hU")/np.maximum(h, 1e-12)
hex, uex, h2, u2, cs = stoker_profile(x, 1.0, 1.0, 0.1)
xf = cross(x, h, 0.5*(h2 + 0.1))
mask = ((x > (u2 - np.sqrt(G*h2)) + 0.3) & (x < xf - 0.3)) if xf else np.zeros(n, bool)
summary["wet"] = dict(h2_exact=h2, u2_exact=u2, cs_exact=cs,
                     h2_num=float(np.median(h[mask])) if mask.any() else None,
                     u2_num=float(np.median(u[mask])) if mask.any() else None,
                     cs_num=float(xf) if xf else None, front=float(xf) if xf else None,
                     h_min=float(h.min()), h_max=float(h.max()))

# dry
h_dry = read_scalar(ROOT / "pipe1D_dry/1/h")
hR_dry, c0 = ritter(x, 1.0, 1.0)
idx = np.where(h_dry > 1e-6)[0]
xf_dry = x[idx[-1]] if idx.size else None
summary["dry"] = dict(c0=c0, front=float(xf_dry) if xf_dry is not None else None,
                      h_min=float(h_dry.min()), h_max=float(h_dry.max()))

# nearly dry
n2 = 2000
x2, dx2 = grid(-160, 160, n2)
h_nd = read_scalar(ROOT / "pipe1D_nearly/16.5/h")
u_nd = read_vector_x(ROOT / "pipe1D_nearly/16.5/hU")/np.maximum(h_nd, 1e-12)
hndex, undex, h2n, u2n, csn = stoker_profile(x2, 16.5, 1.0, 1e-3)
xf_nd = cross(x2, h_nd, 0.5*(h2n + 1e-3))
mask_nd = ((x2 > (u2n - np.sqrt(G*h2n))*16.5 + 0.5) & (x2 < xf_nd - 0.5)) if xf_nd else np.zeros(n2, bool)
summary["nearly_dry"] = dict(h2_exact=h2n, u2_exact=u2n, cs_exact=csn,
                             h2_num=float(np.median(h_nd[mask_nd])) if mask_nd.any() else None,
                             u2_num=float(np.median(u_nd[mask_nd])) if mask_nd.any() else None,
                             cs_num=float(xf_nd/16.5) if xf_nd else None,
                             front=float(xf_nd) if xf_nd else None,
                             h_min=float(h_nd.min()), h_max=float(h_nd.max()))

# radial
n3 = 2000
xr, dr = grid(0, 10, n3)
h_r = read_scalar(ROOT / "radial1D/2/h")
u_r = read_vector_x(ROOT / "radial1D/2/hU")/np.maximum(h_r, 1e-12)
grad = np.abs(np.gradient(h_r, dr))
ifr = int(np.argmax(grad))
summary["radial"] = dict(h_center=float(h_r[0]), h_max=float(h_r.max()),
                         r_front=float(xr[ifr]), h_front=float(h_r[ifr]),
                         u_max=float(np.max(np.abs(u_r))), h_min=float(h_r.min()))

(ROOT / "all_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

fig, ax = plt.subplots(2, 2, figsize=(14, 10))

ax[0,0].plot(x, h, color="#1f77b4", lw=1.3, label="SWEshockFoam")
ax[0,0].plot(x, hex, "--", color="#d62728", lw=1.0, label="Stoker")
ax[0,0].set_xlim(-5, 8); ax[0,0].set_ylim(0, 1.05)
ax[0,0].set_title("Wet pipe, t=1 s"); ax[0,0].set_xlabel("x [m]"); ax[0,0].set_ylabel("depth [m]")
ax[0,0].grid(True, alpha=0.3); ax[0,0].legend(fontsize=8)

ax[0,1].plot(x, h_dry, color="#2ca02c", lw=1.3, label="SWEshockFoam")
ax[0,1].plot(x, hR_dry, "--", color="#d62728", lw=1.0, label="Ritter")
ax[0,1].set_xlim(-5, 8); ax[0,1].set_ylim(0, 1.05)
ax[0,1].set_title("Dry bed, t=1 s"); ax[0,1].set_xlabel("x [m]"); ax[0,1].set_ylabel("depth [m]")
ax[0,1].grid(True, alpha=0.3); ax[0,1].legend(fontsize=8)

ax[1,0].plot(x2, h_nd, color="#9467bd", lw=1.2, label="SWEshockFoam")
ax[1,0].plot(x2, hndex, "--", color="#d62728", lw=1.0, label="Stoker")
ax[1,0].set_xlim(-10, 100); ax[1,0].set_ylim(0, 1.05)
ax[1,0].set_title("Nearly dry (1 mm), t=16.5 s"); ax[1,0].set_xlabel("x [m]"); ax[1,0].set_ylabel("depth [m]")
ax[1,0].grid(True, alpha=0.3); ax[1,0].legend(fontsize=8)

ax[1,1].plot(xr, h_r, color="#ff7f0e", lw=1.3, label="SWEshockFoam")
ax[1,1].set_xlim(0, 8); ax[1,1].set_ylim(0, 1.05)
ax[1,1].set_title("1D radial, t=2 s"); ax[1,1].set_xlabel("r [m]"); ax[1,1].set_ylabel("depth [m]")
ax[1,1].grid(True, alpha=0.3); ax[1,1].legend(fontsize=8)

fig.tight_layout()
fig.savefig(ROOT / "all_cases.png", dpi=150)
print("saved", ROOT / "all_cases.png")
