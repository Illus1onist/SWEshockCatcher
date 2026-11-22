"""Schematic of the 1D pipe and 2D radial dam-break setups."""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

ROOT = Path(__file__).resolve().parent
FIG = ROOT/"figures"
FIG.mkdir(exist_ok=True)

fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.4))

# ---------------- (a) 1D pipe ----------------
x = np.linspace(-10, 10, 800)
ax[0].fill_between(x, 0, np.where(x <= 0, 1.0, 0.1),
                   color="#4c78a8", alpha=0.65, lw=0)
ax[0].axhline(0, color="k", lw=2.0)
ax[0].axvline(0, color="k", ls="--", lw=1.0)
ax[0].annotate(r"$h_L=1$", xy=(-5.0, 0.62), ha="center", fontsize=12)
ax[0].annotate(r"$h_R=0.1$", xy=(5.0, 0.34), ha="center", fontsize=12)
ax[0].annotate("dam", xy=(0.0, 1.10), xytext=(0.8, 1.22),
               arrowprops=dict(arrowstyle="->", lw=0.9), fontsize=10)
ax[0].text(-9.6, -0.12, "transmissive", fontsize=9)
ax[0].text(8.1, -0.12, "transmissive", fontsize=9)
ax[0].set_xlim(-10, 10)
ax[0].set_ylim(-0.2, 1.35)
ax[0].set_xlabel(r"$x$ [m]")
ax[0].set_ylabel(r"$h$ [m]")
ax[0].set_title("(a) 1D pipe: wet dam break")
ax[0].grid(alpha=0.25)

# ---------------- (b) 2D radial tank ----------------
ax[1].add_patch(Rectangle((-10, -10), 20, 20, facecolor="#4c78a8",
                          alpha=0.15, edgecolor="k", lw=2.0))
ax[1].add_patch(Circle((0, 0), 1.0, facecolor="#4c78a8",
                       alpha=0.75, edgecolor="k", lw=0.8))
ax[1].annotate("", xy=(0, 0), xytext=(4.0, 0),
               arrowprops=dict(arrowstyle="->", lw=1.2))
ax[1].text(2.1, 0.55, r"$r$", fontsize=13)
ax[1].text(-0.2, 1.65, r"$h_L$", ha="center", fontsize=13)
ax[1].text(5.3, 2.5, r"$h_R$", fontsize=13)
ax[1].annotate(r"$R_0=1$", xy=(0.7, 0.7), xytext=(2.0, 4.0),
               arrowprops=dict(arrowstyle="->", lw=0.9), fontsize=10)
ax[1].text(-9.4, 9.2, "closed tank", fontsize=9)
ax[1].set_xlim(-10, 10)
ax[1].set_ylim(-10, 10)
ax[1].set_aspect("equal")
ax[1].set_xlabel(r"$x$ [m]")
ax[1].set_ylabel(r"$y$ [m]")
ax[1].set_title("(b) 2D radial tank, $t=0$")
ax[1].grid(alpha=0.25)

fig.tight_layout()
fig.savefig(FIG/"setup.png", dpi=180)
print("saved", FIG/"setup.png")
