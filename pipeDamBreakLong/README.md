# pipeDamBreakLong — 1-cell-wide, 200 m long pipe dam break

Long-channel version of `pipeDamBreak` for tracking the **shock shape over
distance**. The channel is one cell wide in `y` (a true 1D pipe embedded in
a 2D OpenFOAM mesh), 200 m long, and runs to `t = 50 s`. The shock travels
~155 m, so the numerical front has time to develop/smear and can be
compared with the exact Stoker solution at several stations.

## Geometry and mesh

- Channel: `x in [-20, 180] m`, `y in [0, 1] m` (one cell), `z in [0, 0.1] m`.
- `blockMesh`: `2000 x 1 x 1` cells.
- Cell size: `dx = 0.1 m`, `dy = 1 m`, `dz = 0.1 m`.
- Total cells: 2,000 (cheap).
- `frontAndBack` (`z`) = `empty`; the mesh is 2D in `x-y`.
- Side walls: `bottom`, `top` (`y = 0, 1`), `slip`.
- Ends: `left`, `right` (`x = -20, 180`), `zeroGradient` (transmissive).

The dam is at `x = 0`; the initial step is aligned with a cell face.

## Initial condition

- Left state: `h_L = 1.0 m` for `x < 0` (reservoir from `x = -20` to `0`).
- Right state: `h_R = 0.1 m` (wet bed) for `x > 0`.
- `hU = (0, 0, 0)`; flat bed `h0 = 0`; `g = (0, 0, -9.81) m/s^2`.

Set with `setFields` (`boxToCell` from `x = -20` to `0`).

## Exact solution (Stoker, wet bed)

For `h_L = 1`, `h_R = 0.1`, `g = 9.81`:

| Quantity | Value |
|---|---|
| `h_2` (behind shock) | `0.396175 m` |
| `u_2` (behind shock) | `2.321355 m/s` |
| `c_s` (shock speed) | `3.105134 m/s` |
| shock at `t = 50 s` | `x = 155.26 m` |
| `x = 150 m` arrival | `48.31 s` |
| `x = 100 m` arrival | `32.20 s` |
| `x = 50 m` arrival | `16.10 s` |
| `x = 20 m` arrival | `6.44 s` |

The left rarefaction head moves at `-c_L = -3.132 m/s`; it leaves the left
boundary early, and the transmissive end lets it out. The right boundary is
at 180 m; the shock does not reach it within 50 s, so there is no
reflection.

## Boundary conditions

- `h`: `zeroGradient` on `left/right/bottom/top`; `empty` on `frontAndBack`.
- `hU`: `zeroGradient` on `left/right`; `slip` on `bottom/top`; `empty` on
  `frontAndBack`.
- `h0`: `zeroGradient`; `empty` on `frontAndBack`.

## Solver and numerical settings

- Application: `shallowWaterFoam` (OpenFOAM 2406, native).
- `PIMPLE`: `nOuterCorrectors 3`, `nCorrectors 1`,
  `nNonOrthogonalCorrectors 0`, `momentumPredictor yes`.
- `ddt`: `Euler`; `div(phiv,hU)`: `Gauss limitedLinearV 1`;
  `laplacian`: `Gauss linear corrected`.
- `h` PCG/DIC, `hU` smoothSolver/GaussSeidel; final `relTol 0`.
- Time: `endTime 50 s`, `deltaT 0.001 s`, `adjustTimeStep yes`,
  `maxCo 0.5`, `maxDeltaT 0.02 s`, write every `0.5 s`.
- Probes at `x = 20, 50, 100, 150 m` (`y = 0.5`, `z = 0.05`).

## Run

```sh
blockMesh
setFields
shallowWaterFoam
```

or `./Allrun`.

## What to look at

- `h(x,t)` at `t = 10, 30, 50 s`: shock at `31, 93, 155 m`; the plateau
  `h_2`, the rarefaction, the front width.
- The front position `x_f(t)` and its slope; compare with `c_s = 3.105 m/s`.
- The front width vs distance/time (numerical diffusion): does the profile
  reach a steady shape, or keep spreading?
- The shock overshoot/undershoot; `h_min`, `h_max`; any negative `h`.
- Probe arrival times: `20 m` at ~6.4 s, `50 m` at ~16.1 s, `100 m` at
  ~32.2 s, `150 m` at ~48.3 s.
- The plateau `h_2` and `u_2` at late times; compare with the exact values.

## Caveats

- `shallowWaterFoam` is PISO/PIMPLE, not HLLC; the shock is smeared over a
  few cells; the long run shows how the numerical width/overshoot evolves.
- The wet bed is the stable regime; the dry bed (`h_R = 0`) is the stress
  test for positivity and is left for the custom HLLC solver.
- The left boundary is transmissive; the rarefaction leaves the domain.
  The exact Stoker solution assumes an infinite left reservoir, so compare
  the shock/plateau, not the far-left tail, after the rarefaction exits.
- `blockMesh` warns about `convertToMeters` (deprecated in 2406); the mesh
  is correct.

## Observed at t = 50 s (finite-reservoir effect)

The plateau behind the shock decays:

| t [s] | h2 (numerical) | h2 (Stoker) |
|---|---|---|
| 10 | 0.3962 | 0.3962 |
| 20 | 0.3912 | 0.3962 |
| 30 | 0.3028 | 0.3962 |
| 40 | 0.2612 | 0.3962 |
| 50 | 0.2370 | 0.3962 |

`h_max` at `t = 50` is `0.388`, i.e. the initial `h_L = 1` reservoir is
gone. This is **not** a failure of the Stoker comparison: the 20 m left
reservoir is drained by the rarefaction. The exact Stoker solution assumes
an infinite reservoir. The front speed remains close to the exact value
(`3.094` vs `3.105 m/s`), and the numerical front width decreases from
2.5 to 1.4 cells as the shock weakens.

For a long-time Stoker comparison, fix `h = 1` at the left boundary
(infinite reservoir) or extend the left reservoir substantially. With a
fixed left `h`, the plateau stays at `h2 = 0.396` and the front
speed/width can be studied at constant shock strength.

