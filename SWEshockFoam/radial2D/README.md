# radial2D — 2D radial dam break with SWEshockFoam

The 2D Cartesian counterpart of `../radial1D`, on the same 400 x 400
mesh and initial condition as `../../radialDamBreak`:

- domain `[-10,10] x [-10,10] m`, one cell in `z` (2D);
- `h = 1 m` inside `r < 1 m`, `h = 0.1 m` outside;
- `hU = 0`; flat bed;
- closed tank: `h` zero-gradient, `hU` slip on the four outer patches;
- `t = 2 s`, `cfl 0.1`, HLLC, MUSCL/Barth-Jespersen.

The mesh and initial condition are set up by

```sh
./Allrun        # blockMesh + setFields + SWEshockFoam
```

The run takes roughly 10-15 min and writes `h`/`hU` every `0.1 s`.
`../../SWEshockFoam/analyze_radial2d.py` produces the radial average,
the symmetry/mass checks and `radial_1d_vs_2d.png`.

At `t = 2 s`:

| Quantity | Value |
|---|---|
| `h_center` | 0.014110 |
| `h_max` | 0.233014 at `r = 5.575` |
| `r_front` | 5.625 |
| `h_front` | 0.160661 |
| `u_max` | 1.30635 |
| mass error | `1.2e-10` relative |
| symmetry error | `0` (x mirror, y mirror, transpose) |
