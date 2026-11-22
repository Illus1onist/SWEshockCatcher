# pipeDamBreakDry — dry-bed pipe dam break (positivity stress test)

Dry-bed version of the pipe dam break: the right half is **empty**
(`h_R = 0`), the left half is a 1 m deep reservoir. This is the classic
Ritter problem and the most direct stress test for `shallowWaterFoam`:
the dry front has `h -> 0`, and the PISO/height equation has no
positivity limiter. Expect negative/oscillatory `h` near the front, or a
crash; this is exactly the motivation for the custom HLLC/positivity
solver.

## Geometry and mesh

- Channel: `x in [-160, 160] m`, `y in [0, 1] m` (one cell), `z in [0, 0.1] m`.
- `blockMesh`: `1600 x 1 x 1` cells.
- Cell size: `dx = 0.2 m`, `dy = 1 m`, `dz = 0.1 m`.
- Total cells: 1,600.
- `frontAndBack` (`z`) = `empty`; side walls `bottom/top` = `slip`;
  ends `left/right` = `zeroGradient` (transmissive).

The reservoir is the left half (`x in [-160, 0]`), i.e. half the domain.
This is enough for the dry Ritter wave: the front moves at `2 c0` and the
rarefaction head at `-c0`, so the rarefaction is twice slower and does not
drain the reservoir before the front reaches the right end.

## Initial condition

- `h_L = 1.0 m` for `x < 0`; `h_R = 0` for `x > 0`.
- `hU = (0, 0, 0)`; flat bed `h0 = 0`; `g = (0, 0, -9.81) m/s^2`.
- Set with `setFields` (`boxToCell` from `x = -160` to `0`; the default is
  `h = 0`).

## Exact solution (Ritter, dry bed)

With `c0 = sqrt(g h0) = 3.132 m/s`:

```
h = h0                                           x < -c0 t
h = (1/(9g)) (2 c0 - x/t)^2                      -c0 t < x < 2 c0 t
h = 0                                            x > 2 c0 t

u = 0                                            x < -c0 t
u = (2/3) (c0 + x/t)                             -c0 t < x < 2 c0 t
u = 0                                            x > 2 c0 t
```

Front speed `2 c0 = 6.264 m/s`; tail `-c0 = -3.132 m/s`; at `x = 0`:
`h = 4/9 h0 = 0.444 m`, `u = 2/3 c0 = 2.088 m/s`. The front reaches
`x = 160 m` at `t = 25.55 s`; `endTime = 25 s` keeps it just inside the
domain.

## Numerical settings

- `shallowWaterFoam` (2406): `PIMPLE 3/1/0`, `momentumPredictor yes`.
- `Euler`; `div(phiv,hU) = limitedLinearV 1`; `laplacian corrected`.
- `h` PCG/DIC, `hU` smoothSolver/GaussSeidel; final `relTol 0`.
- `endTime 25 s`, `deltaT 0.001`, `maxCo 0.5`, `maxDeltaT 0.02`,
  write every `0.25 s`.
- Probes at `x = 20, 50, 100, 150 m`.

## Run

```sh
blockMesh
setFields
shallowWaterFoam
```

or `./Allrun`.

## What to look for

- `h(x,t)` along the channel: the rarefaction fan, the dry front, the
  exact Ritter profile.
- `h` at the front: does it go negative? What is `h_min`?
- The front position/speed vs `2 c0 t`; the shape/width.
- The probes: the dry front arrival at `20/50/100/150 m`; `h` after the
  front (should be ~0, but may be negative/oscillatory).
- The run log: `FOAM FATAL ERROR`, `h` negative, `NaN`, or a crash.
- Compare with a positivity-preserving HLLC solver; this case is its
  first validation.

## Caveats

- `shallowWaterFoam` is PISO/PIMPLE, not HLLC; the `h` equation is
  elliptic (pressure-like), and positivity of `h` is not guaranteed.
- The dry bed is intentionally the stress case. A clean run is not
  expected; the failure mode is the result.
- The exact Ritter solution is for an infinite reservoir; the left
  reservoir is 160 m, the rarefaction head reaches `-78 m` at `t = 25 s`,
  so it does not leave the domain.
- `blockMesh` warns about `convertToMeters` (deprecated in 2406); the mesh
  is correct.
