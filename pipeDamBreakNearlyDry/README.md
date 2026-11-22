# pipeDamBreakNearlyDry — h_R = 1 mm pipe dam break

Nearly dry version: `h_R = 0.001 m` (1 mm) on the right, `h_L = 1 m` on
the left. The exact `h_R = 0` case makes `shallowWaterFoam` die at
initialisation (`U = hU/h` with `h = 0`); a 1 mm film avoids the division
by zero while keeping the front essentially dry. This is the case to look
for **negative `h`**, oscillations and front-speed errors.

## Geometry and mesh

- `x in [-160, 160] m`, one cell in `y`, one cell in `z` (0.1 m).
- `blockMesh`: `1600 x 1 x 1`; `dx = 0.2 m`; 1,600 cells.
- Side walls `slip`; ends `zeroGradient`; `frontAndBack` = `empty`.

## Initial condition

- `h_L = 1.0 m` for `x < 0`; `h_R = 0.001 m` for `x > 0`.
- `hU = 0`; `h0 = 0`; `g = (0 0 -9.81)`.
- `setFields`: default `h = 0.001`; left box `h = 1`.

## Reference

The exact solution is the wet-bed Stoker solution with `h_R = 0.001 m`
(it approaches Ritter as `h_R -> 0`). Values: `h_2 = 0.066830 m`,
`u_2 = 4.644801 m/s`, `c_s = 4.715359 m/s`. The front is at
`x = 77.72 m` at `t = 16.5 s`, and would reach `x = 117.9 m` at
`t = 25 s`. The 1 mm film makes the front a very strong bore
(`h_2/h_R ~ 67`), not a dry front.

## Observed

At `t = 16.5 s` the field is still clean: `h_min = 0.001`, 0 negative
cells, `h2 = 0.06681` (0.03 % vs Stoker), front `77.72 m`,
speed `4.710 m/s` (0.11 % low), front width `1.176 m = 5.88 cells`.

At `t = 16.725 s` the run dies with `FOAM_SIGFPE` in `Foam::sqrt`:
`h` became negative in the height equation. `Co_max` grew
`0.999 -> 3.24 -> 11.04` in the last three steps. The native PISO/height
solver has no positivity limiter; this is the failure the custom
HLLC/positivity solver must fix.


## Run

```sh
blockMesh
setFields
shallowWaterFoam
```

## What to look for

- `h_min`: does it go negative?
- `h(x,t)` at the front: oscillations, undershoot, front width.
- The front speed vs the Stoker value.
- The probes at `x = 20/50/100/150 m`.
- The log: `FOAM FATAL ERROR`, `NaN`, or a crash.

## Caveats

- `shallowWaterFoam` is PISO/PIMPLE, not positivity-preserving; the 1 mm
  film is a numerical device, not a physical dry bed.
- The exact solution assumes an infinite reservoir; the left reservoir is
  160 m and the rarefaction head is at `-78 m` at `t = 25 s`, so it is
  not drained.
- `blockMesh` warns about `convertToMeters`; the mesh is correct.
