# customSolver — 1D shallow-water solver

A compact 1D SWE solver: **HLL (2-wave) + MUSCL/MC + forward Euler +
positivity/dry tolerance**, Cartesian and axisymmetric (`r`-weighted).

## Files

- `swe1d.py` — the solver (`SWE1D` class).
- `validate_swe1d.py` — Stoker/Ritter/nearly-dry/radial validation.
- `validation_summary.json`, `validation.png` — results.

## Equations

Cartesian:

```
d_t h     + d_x (h u)                 = 0
d_t (h u) + d_x (h u^2 + (g/2) h^2)   = -g h d_x h0
```

Axisymmetric (`r`-weighted):

```
d_t (r h)   + d_r (r h u)                   = 0
d_t (r h u) + d_r (r h u^2 + (g/2) r h^2)   = (g/2) h^2 - g r h d_r h0
```

`u(0) = 0`; the first face at `r=0` has zero flux.

## Design choices

1. **Python/NumPy prototype** — fast iteration; the 1D reference for the
   eventual OpenFOAM port.
2. **Cell-centered FV, uniform grid, one ghost cell per side** — simple and
   conservative; BCs are imposed in the ghosts.
3. **HLL for 1D** — the 1D SWE has two waves `u±c`; HLL is natural,
   positive and robust. HLLC is reserved for 2D, where the contact/shear
   wave exists.
4. **MUSCL + MC limiter in primitive (`h`, `u`)** — second order, TVD, less
   diffusive than minmod; primitive variables are safer for `h -> 0`.
5. **Davis wave speeds** — `S_L = min(uL-cL, uR-cR)`,
   `S_R = max(uL+cL, uR+cR)`; cheap and positivity-friendly.
6. **Positivity / dry tolerance** — `h_min`, `dry_tol`; `u = 0` on dry
   cells; `h` floor at interfaces.
7. **Forward Euler, CFL 0.3** — first order in time, second in space; enough
   for the dam break; generous CFL margin.
8. **Flat bed `h0 = 0` for the dam break** — the `h0` source is implemented
   but not exercised; a well-balanced bump case is a later step.
9. **BCs** — transmissive, wall, symmetry; symmetry sets `u(0)=0`.
10. **Radial `r`-weighted form** — removes the `1/r` singularity and uses
    the same HLL/MUSCL machinery.

## Validation

| Case | Reference | Result |
|---|---|---|
| wet 10:1, t=1 s | Stoker `h2=0.39617`, `u2=2.32135`, `cs=3.10513` | `h2=0.39622` (+0.01 %), `u2=2.32180`, `cs=3.10647` (+0.04 %) |
| dry, t=1 s | Ritter front `2c0=6.264` | front `~7.4` by `h>1e-6` (smeared ahead); `h_min=1e-10`, no negative |
| nearly dry `hR=1e-3`, t=16.5 s | Stoker `h2=0.06683`, `u2=4.64480`, `cs=4.71536` | `h2=0.06339` (−5.2 %), `u2=4.64701`, `cs=4.72904` (+0.29 %), front `78.0 m`; `h_min=0.001`, no negative |
| 1D radial, t=2 s | no closed form | `h_center=0.0138`, `h_max=0.2766`, `r_front=5.64`, `h_front=0.1734`, `u_max=1.407`; **oscillatory** (superseded by `SWEshockFoam/radial1D`) |

The wet and nearly-dry results are good. The dry case no longer crashes
(`shallowWaterFoam` did) and stays positive, although the front is
smeared by the HLL diffusion; the `h>1e-6` front metric overestimates it.
The radial case is still **oscillatory** — the `r`-weighted source /
reconstruction needs Step-3 work (primitive reconstruction, well-balanced
radial source, possibly HLLC).  The OpenFOAM version
(`SWEshockFoam/radial1D`) supersedes this prototype and is much smoother.

## Usage

```python
from swe1d import SWE1D
s = SWE1D(x, g=9.81, geometry="cartesian", riemann="hll", limiter="mc")
s.set_state(h, u)
s.run(2.0, ("transmissive", "transmissive"))
```
