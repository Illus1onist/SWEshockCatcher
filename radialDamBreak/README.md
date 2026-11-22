# radialDamBreak — 2D radial dam break (shallow water)

Baseline OpenFOAM case for the radial dam-break verification problem.
It is written for the Foundation solver `shallowWaterFoam` (PISO/PIMPLE
for the shallow-water equations, **not** an HLLC shock-capturing solver).
The same `0/`, `constant/` and geometry are intended to be reused by the
custom HLLC/MUSCL solver that will be developed in this project.  The
custom-solver version of this case is
[`../SWEshockFoam/radial2D`](../SWEshockFoam/radial2D); it uses the same
400 x 400 mesh and the same `cylinderToCell` initial condition.

## Geometry and mesh

- 2D square domain: `x, y in [-10, 10] m`, one cell in `z` (0 to 0.1 m).
- `blockMesh`: `400 x 400 x 1` cells.
- Cell size: `0.05 x 0.05 x 0.1 m`.
- Total cells: 160,000.
- `frontAndBack` patches: `empty` (2D).
- Outer patches: `left`, `right`, `bottom`, `top` (all `patch`).

The initial water column has radius `R0 = 1 m`, i.e. 20 cells across the
radius at this resolution. The outer boundary is 10 radii away, so the
radial bore is not affected by reflections during the first ~1.5 s.

## Initial condition

- Circular column: `h = 1.0 m` for `r <= 1 m`.
- Ambient depth: `h = 0.1 m`.
- `hU = (0, 0, 0)` everywhere.
- Flat bed: `h0 = 0` (`constant/h0`).
- Gravity: `g = (0, 0, -9.81) m/s^2`.

The column is set with `setFields` (`system/setFieldsDict`,
`cylinderToCell`, radius 1 m). The `0/h` and `0/hU` files contain the
ambient state; `setFields` overwrites the column.

For the **dry** radial dam break, change the ambient depth in
`system/setFieldsDict` from `0.1` to `0` and in `0/h` from `0.1` to `0`.
This is the stress case for `shallowWaterFoam`: expect negative/oscillatory
`h` near the front; the custom positivity-preserving HLLC solver is the
target for this variant.

## Boundary conditions

- `h`: `zeroGradient` on `left`, `right`, `bottom`, `top`; `empty` on
  `frontAndBack`.
- `hU`: `slip` on the outer patches; `empty` on `frontAndBack`.
- `h0`: `zeroGradient` on the outer patches; `empty` on `frontAndBack`.

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
- Time: `endTime 2 s`, `deltaT 0.001 s`, `adjustTimeStep yes`,
  `maxCo 0.5`, `maxDeltaT 0.01 s`, write every `0.05 s`.

## Run

```sh
blockMesh
setFields
shallowWaterFoam
```

or, on a machine with the OpenFOAM environment:

```sh
./Allrun
```

The case is not run automatically.

## What to look at

- `h(x,y,t)` and `|hU|/h`: the circular bore front, its radius and speed.
- Radial symmetry: profiles along `x` and `y` should coincide.
- `h` at the centre and at the front; the front radius `R(t)`.
- 2D maps in ParaView (`h`, `U`, `hTotal = h + h0`).
- The `shallowWaterFoam` baseline is expected to smear/oscillate at the
  bore; this is the motivation for the custom HLLC/MUSCL solver.

## Reference

The radial dam break has a self-similar axisymmetric shallow-water
solution; it is the radial counterpart of the Ritter solution. Useful
references:

- Ritter (1892), dry-bed dam break (Cartesian).
- Stoker (1957), wet-bed dam break.
- Delestre et al., SWASHES: a library of shallow-water analytic solutions
  (arXiv:1110.0288) — Cartesian dam-break formulas and code.
- The radial dam-break similarity solution: compare the front radius and
  the `h(r,t)` profile against a 1D radial solver or the literature.

The 1D Cartesian Ritter solution is **not** the radial solution; use it
only as a sanity check along a radial line at early times, or use a 1D
radial reference.

## Caveats

- `shallowWaterFoam` is a PISO/PIMPLE pressure-based solver, not a
  Riemann/HLLC solver. It is stable for smooth/subcritical flows
  (`squareBump` tutorial) but is not designed for strong bores or dry
  beds.
- The `h` equation is elliptic (pressure-like), not a conservative flux;
  positivity of `h` is not guaranteed.
- The case is written for Foundation OpenFOAM. In ESI OpenFOAM 14 the
  solver may be absent; check `foamRun -list` or use the custom solver.
- The custom HLLC/MUSCL solver will reuse the same `0/`, `constant/` and
  geometry; only the `application` in `system/controlDict` and the
  numerical settings will change.
