"""Validate SWE1D: Stoker (wet), Ritter (dry), nearly dry, radial."""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from swe1d import SWE1D

G = 9.81
OUT = Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)


def stoker(hL, hR, g=G):
    cL = np.sqrt(g * hL)

    def f(h2):
        c2 = np.sqrt(g * h2)
        u2 = 2.0 * (cL - c2)
        rhs = (h2 - hR) * np.sqrt(0.5 * g * (h2 + hR) / (h2 * hR))
        return u2 - rhs

    a, b = hR, hL
    for _ in range(200):
        m = 0.5 * (a + b)
        if f(a) * f(m) <= 0.0:
            b = m
        else:
            a = m
    h2 = 0.5 * (a + b)
    c2 = np.sqrt(g * h2)
    u2 = 2.0 * (cL - c2)
    cs = h2 * u2 / (h2 - hR)
    return h2, u2, cs


def stoker_profile(x, t, hL, hR, g=G):
    h2, u2, cs = stoker(hL, hR, g)
    cL = np.sqrt(g * hL)
    c2 = np.sqrt(g * h2)
    h = np.full_like(x, hR)
    u = np.zeros_like(x)
    h[x < -cL * t] = hL
    fan = (x >= -cL * t) & (x <= (u2 - c2) * t)
    c = (2.0 * cL - x[fan] / t) / 3.0
    h[fan] = c * c / g
    u[fan] = (2.0 / 3.0) * (cL + x[fan] / t)
    plateau = (x > (u2 - c2) * t) & (x < cs * t)
    h[plateau] = h2
    u[plateau] = u2
    return h, u, h2, u2, cs


def ritter(x, t, hL, g=G):
    c0 = np.sqrt(g * hL)
    h = np.where(x < -c0 * t, hL,
        np.where(x < 2.0 * c0 * t, (1.0 / (9.0 * g)) * (2.0 * c0 - x / t) ** 2, 0.0))
    u = np.where(x < -c0 * t, 0.0,
        np.where(x < 2.0 * c0 * t, (2.0 / 3.0) * (c0 + x / t), 0.0))
    return h, u, c0


def cross(x, h, level):
    for i in range(len(x) - 1, 0, -1):
        if h[i] < level <= h[i - 1]:
            x0, x1 = x[i - 1], x[i]
            h0, h1 = h[i - 1], h[i]
            return x0 + (level - h0) * (x1 - x0) / (h1 - h0)
    return None


def run_pipe(n, xmin, xmax, hL, hR, t_end):
    dx = (xmax - xmin) / n
    x = np.linspace(xmin + 0.5 * dx, xmax - 0.5 * dx, n)
    h = np.where(x < 0.0, hL, hR)
    u = np.zeros(n)
    s = SWE1D(x, g=G, geometry="cartesian", riemann="hll", limiter="mc",
              cfl=0.3, dry_tol=1e-8, h_min=1e-10)
    s.set_state(h, u)
    s.run(t_end, ("transmissive", "transmissive"))
    return s


def run_radial(n, rmax, hL, hR, t_end):
    dr = rmax / n
    r = np.linspace(0.5 * dr, rmax - 0.5 * dr, n)
    h = np.where(r < 1.0, hL, hR)
    u = np.zeros(n)
    s = SWE1D(r, g=G, geometry="radial", riemann="hll", limiter="mc",
              cfl=0.3, dry_tol=1e-8, h_min=1e-10)
    s.set_state(h, u)
    s.run(t_end, ("symmetry", "transmissive"))
    return s


summary = {}

# 1) wet pipe
s = run_pipe(2000, -10, 10, 1.0, 0.1, 1.0)
x = s.x
h, u, h2, u2, cs = stoker_profile(x, 1.0, 1.0, 0.1)
xf = cross(x, s.h, 0.5 * (h2 + 0.1))
mask = ((x > (u2 - np.sqrt(G * h2)) * 1.0 + 0.2) & (x < xf - 0.2)) if xf else np.zeros_like(x, bool)
summary["wet"] = dict(h2_exact=h2, u2_exact=u2, cs_exact=cs,
                      h2_num=float(np.median(s.h[mask])) if mask.any() else None,
                      u2_num=float(np.median(s.u[mask])) if mask.any() else None,
                      cs_num=float(xf / 1.0) if xf else None,
                      front=float(xf) if xf else None,
                      h_min=float(s.h.min()), h_max=float(s.h.max()),
                      steps=s.steps)

# 2) dry pipe
s_dry = run_pipe(2000, -10, 10, 1.0, 0.0, 1.0)
x_dry = s_dry.x
h_dry, u_dry, c0 = ritter(x_dry, 1.0, 1.0)
idx_dry = np.where(s_dry.h > 1e-6)[0]
xf_dry = x_dry[idx_dry[-1]] if idx_dry.size else None
summary["dry"] = dict(c0=c0, front=float(xf_dry) if xf_dry is not None else None,
                      cs_num=float(xf_dry / 1.0) if xf_dry is not None else None,
                      h_min=float(s_dry.h.min()), h_max=float(s_dry.h.max()),
                      steps=s_dry.steps)

# 3) nearly dry pipe
s_nd = run_pipe(2000, -160, 160, 1.0, 1e-3, 16.5)
x_nd = s_nd.x
h_nd, u_nd, h2_nd, u2_nd, cs_nd = stoker_profile(x_nd, 16.5, 1.0, 1e-3)
xf_nd = cross(x_nd, s_nd.h, 0.5 * (h2_nd + 1e-3))
mask_nd = ((x_nd > (u2_nd - np.sqrt(G * h2_nd)) * 16.5 + 0.5) & (x_nd < xf_nd - 0.5)) if xf_nd else np.zeros_like(x_nd, bool)
summary["nearly_dry"] = dict(h2_exact=h2_nd, u2_exact=u2_nd, cs_exact=cs_nd,
                             h2_num=float(np.median(s_nd.h[mask_nd])) if mask_nd.any() else None,
                             u2_num=float(np.median(s_nd.u[mask_nd])) if mask_nd.any() else None,
                             cs_num=float(xf_nd / 16.5) if xf_nd else None,
                             front=float(xf_nd) if xf_nd else None,
                             h_min=float(s_nd.h.min()), h_max=float(s_nd.h.max()),
                             steps=s_nd.steps)

# 4) radial
s_r = run_radial(2000, 10.0, 1.0, 0.1, 2.0)
xr = s_r.x
grad = np.abs(np.gradient(s_r.h, s_r.dx))
ifront = int(np.argmax(grad))
summary["radial"] = dict(h_center=float(s_r.h[0]), h_max=float(s_r.h.max()),
                         r_front=float(xr[ifront]), h_front=float(s_r.h[ifront]),
                         u_max=float(np.max(np.abs(s_r.u))),
                         h_min=float(s_r.h.min()), steps=s_r.steps)

(OUT / "validation_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

# ------------------------------------------------------------------ figure
fig, axes = plt.subplots(2, 3, figsize=(18, 9))

ax = axes[0, 0]
ax.plot(x, s.h, color="#1f77b4", lw=1.4, label="num")
ax.plot(x, h, "--", color="#d62728", lw=1.0, label="Stoker")
ax.set_xlim(-5, 8)
ax.set_ylim(0, 1.05)
ax.set_xlabel("x [m]")
ax.set_ylabel("h [m]")
ax.set_title("Wet pipe: h(x), t=2 s")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

ax = axes[0, 1]
ax.plot(x, s.u, color="#1f77b4", lw=1.2, label="num")
ax.plot(x, u, "--", color="#d62728", lw=1.0, label="Stoker")
ax.set_xlim(-5, 8)
ax.set_xlabel("x [m]")
ax.set_ylabel("u [m/s]")
ax.set_title("Wet pipe: u(x), t=2 s")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

ax = axes[0, 2]
ax.plot(x_dry, s_dry.h, color="#2ca02c", lw=1.4, label="num")
ax.plot(x_dry, h_dry, "--", color="#d62728", lw=1.0, label="Ritter")
ax.set_xlim(-5, 8)
ax.set_ylim(0, 1.05)
ax.set_xlabel("x [m]")
ax.set_ylabel("h [m]")
ax.set_title("Dry pipe: h(x), t=2 s")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

ax = axes[1, 0]
ax.plot(x_nd, s_nd.h, color="#9467bd", lw=1.2, label="num")
ax.plot(x_nd, h_nd, "--", color="#d62728", lw=1.0, label="Stoker")
ax.set_xlim(-10, 100)
ax.set_ylim(0, 1.05)
ax.set_xlabel("x [m]")
ax.set_ylabel("h [m]")
ax.set_title("Nearly dry (hR=1 mm): h(x), t=16.5 s")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

ax = axes[1, 1]
ax.plot(xr, s_r.h, color="#ff7f0e", lw=1.4, label="num")
ax.set_xlim(0, 8)
ax.set_ylim(0, 1.05)
ax.set_xlabel("r [m]")
ax.set_ylabel("h [m]")
ax.set_title("1D radial: h(r), t=2 s")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

ax = axes[1, 2]
ax.plot(xr, s_r.u, color="#ff7f0e", lw=1.4, label="num")
ax.set_xlim(0, 8)
ax.set_xlabel("r [m]")
ax.set_ylabel("u [m/s]")
ax.set_title("1D radial: u(r), t=2 s")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

fig.tight_layout()
fig.savefig(OUT / "validation.png", dpi=150)
print("saved", OUT / "validation.png")
