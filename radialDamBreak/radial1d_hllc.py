"""Standalone 1D axisymmetric shallow-water solver for the radial dam break.

Equations (r-weighted conservative form):

    d_t (r h)     + d_r (r h u)                 = 0
    d_t (r h u)   + d_r (r h u^2 + (g/2) r h^2) = (g/2) h^2

with h = h(r,t), u = u(r,t), r >= 0.  The r-weighting removes the 1/r
singularity; the geometric source is S = (g/2) h^2 in the momentum.

Numerics: finite volume, HLLC flux, MUSCL reconstruction with the MC
limiter, positivity floor, symmetry at r=0, transmissive outer boundary.
This is the same family as the planned custom OpenFOAM solver.
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G = 9.81
HL = 1.0
HR = 0.1
R0 = 1.0
TEND = 2.0
RMAX = 10.0
N = 2000
CFL = 0.4

OUT = Path(__file__).resolve().parent
dr = RMAX / N
r = (np.arange(N) + 0.5) * dr          # cell centres
rf = np.arange(N + 1) * dr             # faces, rf[0] = 0

h = np.where(r < R0, HL, HR)
u = np.zeros(N)
q0 = r * h
q1 = r * h * u


def minmod(a, b, c):
    return np.where((a > 0) & (b > 0) & (c > 0), np.minimum(np.minimum(a, b), c),
           np.where((a < 0) & (b < 0) & (c < 0), np.maximum(np.maximum(a, b), c), 0.0))


def slopes(q):
    dqL = q - np.roll(q, 1)
    dqR = np.roll(q, -1) - q
    return minmod(2.0 * dqL, 0.5 * (np.roll(q, -1) - np.roll(q, 1)), 2.0 * dqR)


def hllc_flux(hL, uL, hR, uR, rface):
    cL = np.sqrt(G * hL)
    cR = np.sqrt(G * hR)
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)
    pL = 0.5 * G * hL * hL
    pR = 0.5 * G * hR * hR
    num = pR - pL + hL * uL * (SL - uL) - hR * uR * (SR - uR)
    den = hL * (SL - uL) - hR * (SR - uR)
    SM = np.where(np.abs(den) > 1e-14, num / np.where(np.abs(den) > 1e-14, den, 1.0), 0.5 * (uL + uR))

    qL = np.stack([hL, hL * uL], axis=-1)
    qR = np.stack([hR, hR * uR], axis=-1)
    FL = np.stack([hL * uL, hL * uL * uL + 0.5 * G * hL * hL], axis=-1)
    FR = np.stack([hR * uR, hR * uR * uR + 0.5 * G * hR * hR], axis=-1)

    hsL = hL * (SL - uL) / np.where(np.abs(SL - SM) > 1e-14, SL - SM, 1.0)
    hsR = hR * (SR - uR) / np.where(np.abs(SR - SM) > 1e-14, SR - SM, 1.0)
    hsL = np.maximum(hsL, 1e-12)
    hsR = np.maximum(hsR, 1e-12)
    qsL = np.stack([hsL, hsL * SM], axis=-1)
    qsR = np.stack([hsR, hsR * SM], axis=-1)
    FsL = FL + SL[..., None] * (qsL - qL)
    FsR = FR + SR[..., None] * (qsR - qR)

    F = np.where((SL >= 0)[..., None], FL,
        np.where((SM >= 0)[..., None], FsL,
        np.where((SR > 0)[..., None], FsR, FR)))
    return F * rface[..., None]


t = 0.0
history = {0.0: (h.copy(), u.copy())}
step = 0
while t < TEND:
    h = np.maximum(q0 / r, 1e-12)
    u = q1 / np.maximum(q0, 1e-30)
    u[0] = 0.0
    q1[0] = 0.0
    c = np.sqrt(G * h)
    dt = CFL * dr / np.max(np.abs(u) + c)
    if t + dt > TEND:
        dt = TEND - t

    s0 = slopes(q0)
    s1 = slopes(q1)
    # symmetry at r = 0: q0 = r h is odd, q1 = r h u is even
    s0[0] = minmod(np.array([4.0 * q0[0]]),
                   np.array([0.5 * (q0[1] + q0[0])]),
                   np.array([2.0 * (q0[1] - q0[0])]))[0]
    s1[0] = minmod(np.array([0.0]),
                   np.array([0.5 * (q1[1] - q1[0])]),
                   np.array([2.0 * (q1[1] - q1[0])]))[0]
    q0L = q0 + 0.5 * s0
    q0R = np.roll(q0, -1) - 0.5 * np.roll(s0, -1)
    q1L = q1 + 0.5 * s1
    q1R = np.roll(q1, -1) - 0.5 * np.roll(s1, -1)

    rface = rf[1:-1]
    hL = np.maximum(q0L[:-1] / rface, 1e-12)
    hR = np.maximum(q0R[:-1] / rface, 1e-12)
    uL = q1L[:-1] / np.maximum(q0L[:-1], 1e-30)
    uR = q1R[:-1] / np.maximum(q0R[:-1], 1e-30)

    F = np.zeros((N + 1, 2))
    F[0] = 0.0
    F[1:-1] = hllc_flux(hL, uL, hR, uR, rface)
    # transmissive outer: same state as the last cell
    hOut = np.maximum(q0[-1] / r[-1], 1e-12)
    uOut = q1[-1] / max(q0[-1], 1e-30)
    F[-1] = hllc_flux(np.array([hOut]), np.array([uOut]), np.array([hOut]), np.array([uOut]),
                      np.array([RMAX]))[0]

    S = np.zeros((N, 2))
    S[:, 1] = 0.5 * G * h * h

    q0 = q0 - dt / dr * (F[1:, 0] - F[:-1, 0])
    q1 = q1 - dt / dr * (F[1:, 1] - F[:-1, 1]) + dt * S[:, 1]

    q0 = np.maximum(q0, 1e-12 * r)
    q1[0] = 0.0
    t += dt
    step += 1
    for tt in (0.5, 1.0, 2.0):
        if abs(t - tt) < 0.5 * dt and tt not in history:
            h = np.maximum(q0 / r, 1e-12)
            u = q1 / np.maximum(q0, 1e-30)
            u[0] = 0.0
            history[tt] = (h.copy(), u.copy())

h = np.maximum(q0 / r, 1e-12)
u = q1 / np.maximum(q0, 1e-30)
u[0] = 0.0
history[TEND] = (h.copy(), u.copy())

# front: max |dh/dr|
grad = np.abs(np.gradient(h, dr))
ifront = int(np.argmax(grad))
r_front = r[ifront]
h_front = h[ifront]
u_front = u[ifront]
h_center = h[0]
u_max = float(np.max(np.abs(u)))

summary = {
    "N": N, "dr": dr, "t_end": TEND, "R0": R0, "hL": HL, "hR": HR,
    "h_center": float(h_center), "h_max": float(np.max(h)),
    "r_front": float(r_front), "h_front": float(h_front), "u_front": float(u_front),
    "u_max": u_max, "steps": step,
    "h_at_r1": float(np.interp(1.0, r, h)), "u_at_r1": float(np.interp(1.0, r, u)),
}
OUT.mkdir(exist_ok=True)
(OUT / "radial1d_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

fig, ax = plt.subplots(1, 2, figsize=(15, 5.5))
for tt, c in ((0.0, "#7f7f7f"), (0.5, "#1f77b4"), (1.0, "#2ca02c"), (2.0, "#d62728")):
    hh, uu = history[tt]
    ax[0].plot(r, hh, color=c, lw=1.6, label=f"t={tt:g}")
    ax[1].plot(r, uu, color=c, lw=1.6, label=f"t={tt:g}")
ax[0].set_xlim(0, 8)
ax[0].set_ylim(0, 1.05)
ax[0].set_xlabel("r [m]")
ax[0].set_ylabel("h [m]")
ax[0].set_title("1D radial SWE: depth")
ax[0].grid(True, alpha=0.3)
ax[0].legend(fontsize=8)
ax[1].set_xlim(0, 8)
ax[1].set_xlabel("r [m]")
ax[1].set_ylabel("u [m/s]")
ax[1].set_title("1D radial SWE: radial velocity")
ax[1].grid(True, alpha=0.3)
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "radial1d_solution.png", dpi=150)
print("saved", OUT / "radial1d_solution.png")
