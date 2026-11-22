# pipeDamBreak — 2D Cartesian channel dam break (shallow water)

Baseline OpenFOAM case for the **pipe** (Cartesian) dam break. It is the
1D Ritter/Stoker problem embedded in a 2D channel: the flow is uniform in
`y`, so the exact 1D solution can be compared along any line `y = const`.

The case is written for the Foundation solver `shallowWaterFoam`
(PISO/PIMPLE for the shallow-water equations, **not** an HLLC
shock-capturing solver). The same fields and geometry are intended to be
reused by the custom HLLC/MUSCL solver.

## Geometry and mesh

- Channel: `x in [-5, 5] m`, `y in [0, 1] m`, one cell in `z` (0 to 0.1 m).
- `blockMesh`: `400 x 10 x 1` cells.
- Cell size: `dx = 0.025 m`, `dy = 0.1 m`, `dz = 0.1 m`.
- Total cells: 4,000 (very cheap).
- `frontAndBack` patches: `empty` (2D).
- Side walls: `bottom`, `top` (`y = 0, 1`).
- Ends: `left`, `right` (`x = -5, 5`).

The dam is at `x = 0`. The initial step is aligned with a cell face
(`x = 0` is a face of the 400-cell grid).

## Initial condition

- Left state: `h_L = 1.0 m` for `x < 0`.
- Right state: `h_R = 0.1 m` (wet bed) for `x > 0`.
- `hU = (0, 0, 0)` everywhere.
- Flat bed: `h0 = 0` (`constant/h0`).
- Gravity: `g = (0, 0, -9.81) m/s^2`.

The left state is set with `setFields` (`system/setFieldsDict`,
`boxToCell` from `x=-5` to `x=0`). The `0/h` and `0/hU` contain the
ambient right state; `setFields` overwrites the left half.

For the **dry** dam break (classic Ritter), change `h_R` from `0.1` to
`0` in both `system/setFieldsDict` and `0/h`. This is the stress case for
`shallowWaterFoam`; the custom positivity-preserving HLLC solver is the
target for it.

## Boundary conditions

- `h`: `zeroGradient` on `left`, `right`, `bottom`, `top`; `empty` on
  `frontAndBack`.
- `hU`: `zeroGradient` on `left`/`right` (transmissive ends); `slip` on
  `bottom`/`top` (side walls); `empty` on `frontAndBack`.
- `h0`: `zeroGradient` on all outer patches; `empty` on `frontAndBack`.

## Exact solution

### Dry bed: Ritter (1892)

With `h_L = h0`, `h_R = 0`, `u = 0`, `c0 = sqrt(g h0)`:

```
h = h0                                           x < -c0 t
h = (1/(9g)) (2 c0 - x/t)^2                      -c0 t < x < 2 c0 t
h = 0                                            x > 2 c0 t

u = 0                                            x < -c0 t
u = (2/3) (c0 + x/t)                             -c0 t < x < 2 c0 t
u = 0                                            x > 2 c0 t
```

Front speed `2 c0`, tail `-c0`; at `x = 0`: `h = 4/9 h0`,
`u = 2/3 c0`.

### Wet bed: Stoker (1957)

For `h_R > 0` the solution is a rarefaction (left) plus a bore
(right); the shock speed and the state behind the shock follow from the
Rankine-Hugoniot conditions. The formulas and code are in
[SWASHES](https://arxiv.org/abs/1110.0288), section 4.1.

The 1D solution is exact for the **pipe**; the 2D channel is uniform in
`y`, so comparing along `y = 0.5` is meaningful.

## Solver and numerical settings

- Application: `shallowWaterFoam` (Foundation).
- `PIMPLE`: `nOuterCorrectors 3`, `nCorrectors 1`,
  `nNonOrthogonalCorrectors 0`, `momentumPredictor yes`.
- `ddt`: `Euler`.
- `div(phiv,hU)`: `Gauss limitedLinearV 1`.
- `laplacian`: `Gauss linear corrected`.
- `grad`: `Gauss linear`.
- Linear solvers: `h` — PCG/DIC, `tol 1e-10`, `relTol 0.01`;
  `hU` — smoothSolver/GaussSeidel, `tol 1e-10`, `relTol 0.1`; final
  correctors use `relTol 0`.
- Time: `endTime 1 s`, `deltaT 0.0005 s`, `adjustTimeStep yes`,
  `maxCo 0.5`, `maxDeltaT 0.002 s`, write every `0.02 s`.

## Run

```sh
blockMesh
setFields
shallowWaterFoam
```

or `./Allrun`.

## What to look at

- `h(x,t)`, `u(x,t)` along `y = 0.5`: rarefaction fan, shock/dry front,
  the exact Ritter/Stoker profiles.
- Front position and speed; the depth at `x = 0`.
- `Fr = |u|/sqrt(gh)` and `lambda +/- = u +/- 2 sqrt(gh)`.
- Probes at `x = -2, 2, 4` (`y = 0.5`).
- 2D uniformity: profiles along different `y` should coincide (the
  channel is uniform); any `y` dependence is a 2D/numerical artifact.
- Grid convergence and limiters (`limitedLinearV` vs `vanLeerV` vs
  `upwind`); the front width.

## Caveats

- `shallowWaterFoam` is a PISO/PIMPLE pressure-based solver, not a
  Riemann/HLLC solver. It is stable for smooth/subcritical flows
  (`squareBump` tutorial) but is not designed for strong shocks or dry
  beds.
- The `h` equation is elliptic (pressure-like), not a conservative flux;
  positivity of `h` is not guaranteed.
- The case is written for Foundation OpenFOAM. In ESI OpenFOAM 2406
  `shallowWaterFoam` is absent; use Foundation or the custom solver.
  `blockMesh`/`checkMesh` work in 2406.
- The ends are transmissive (`zeroGradient`); the domain is 5 m. For the
  dry case the front reaches `x = 5 m` at about `0.8 s`; `endTime 1 s`
  includes the arrival/reflection at the right end. Reduce `endTime` or
  enlarge the domain for the pre-reflection window.
- The 1D exact solution is exact for the pipe; the 2D channel is uniform
  in `y`, so the comparison is meaningful. For the radial case use the
  1D radial reference instead.
