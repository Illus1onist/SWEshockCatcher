"""Grid-convergence study for the wet 1D dam break at t = 1 s.

Reads convergence/results/h_N#### and hU_N#### (copies of the OpenFOAM
ascii fields at t = 1), compares with the Stoker solution on the cell
centres, and writes convergence.json / convergence.png.
"""
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
G = 9.81
HL, HR = 1.0, 0.1
TEND = 1.0
XMIN, XMAX = -10.0, 10.0


def read_scalar(path):
    txt = Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;", txt, re.S)
    if m:
        return np.array([float(v) for v in m.group(2).split()])
    m = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", txt)
    return np.full(1, float(m.group(1)))


def read_vector_x(path):
    txt = Path(path).read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(\s*(.*?)\s*\)\s*;", txt, re.S)
    if m:
        return np.array([float(v.split()[0]) for v in re.findall(r"\(([^()]*)\)", m.group(2))])
    m = re.search(r"internalField\s+uniform\s+\(([^()]*)\)\s*;", txt)
    return np.full(1, float(m.group(1).split()[0]))


def stoker(hL, hR, g=G):
    cL = np.sqrt(g*hL)

    def f(h):
        c = np.sqrt(g*h)
        return 2.0*(cL - c) - (h - hR)*np.sqrt(0.5*g*(h + hR)/(h*hR))

    a, b = hR, hL
    for _ in range(200):
        m = 0.5*(a + b)
        if f(a)*f(m) <= 0:
            b = m
        else:
            a = m
    h2 = 0.5*(a + b)
    c2 = np.sqrt(g*h2)
    u2 = 2.0*(cL - c2)
    cs = h2*u2/(h2 - hR)
    return h2, u2, cs


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


Ns = [128, 256, 512, 1024, 2048]
rows = []
for N in Ns:
    h = read_scalar(ROOT/"results"/f"h_N{N}")
    hu = read_vector_x(ROOT/"results"/f"hU_N{N}")
    dx = (XMAX - XMIN)/N
    x = XMIN + (np.arange(N) + 0.5)*dx
    he, ue, h2, u2, cs = stoker_profile(x, TEND, HL, HR)
    hue = he*ue
    e_h = np.mean(np.abs(h - he))
    e_hu = np.mean(np.abs(hu - hue))
    e_u = np.mean(np.abs(hu/np.maximum(h, 1e-12) - ue))
    rows.append(dict(N=N, dx=dx, L1_h=float(e_h), L1_hu=float(e_hu),
                     L1_u=float(e_u), Linf_h=float(np.max(np.abs(h - he))),
                     h_min=float(h.min()), h_max=float(h.max())))

for i in range(len(rows) - 1):
    rows[i]["p_h"] = float(np.log(rows[i]["L1_h"]/rows[i+1]["L1_h"])/np.log(2.0))
    rows[i]["p_hu"] = float(np.log(rows[i]["L1_hu"]/rows[i+1]["L1_hu"])/np.log(2.0))
rows[-1]["p_h"] = None
rows[-1]["p_hu"] = None

summary = dict(Ns=Ns, hL=HL, hR=HR, t=TEND, h2=h2, u2=u2, cs=cs, rows=rows)
(ROOT/"convergence.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

fig, ax = plt.subplots(figsize=(7.2, 5.2))
dx = np.array([r["dx"] for r in rows])
eh = np.array([r["L1_h"] for r in rows])
ehu = np.array([r["L1_hu"] for r in rows])
ax.loglog(dx, eh, "o-", lw=1.5, label=r"$L_1(h)$")
ax.loglog(dx, ehu, "s--", lw=1.5, label=r"$L_1(hu)$")
# reference slopes
c1 = eh[0]/dx[0]
c2 = ehu[0]/dx[0]
ax.loglog(dx, c1*dx, ":", color="gray", lw=1.0, label=r"$\propto \Delta x$")
ax.loglog(dx, c2*dx, ":", color="gray", lw=1.0)
ax.set_xlabel(r"$\Delta x$ [m]")
ax.set_ylabel(r"$L_1$ error")
ax.set_title("Wet dam break, t = 1 s: grid convergence")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(ROOT/"convergence.png", dpi=160)
print("saved", ROOT/"convergence.png")
