# Plan: 1D pipe -> custom solver -> 1D radial -> 2D radial

This is the working plan. Everything else in the project is either a
supporting case, a reference, or out of scope.

## Goal

A verified **2D radial dam break** computed with **our own SWE solver**,
built in four steps. The 2D radial step is deliberately last.

## Step 1 — native 1D OpenFOAM SWE: baseline and failure modes

**What:** run the native `shallowWaterFoam` on the 1D (one cell in `y`)
pipe dam break and document where and how it breaks.

**Cases:**

- `pipeDamBreak` — 2D channel, 10 cells in `y`, wet 10:1; the "standard"
  baseline (effectively 1D, uniform in `y`).
- `pipeDamBreakLong` — 1 cell in `y`, 200 m long, wet 10:1; the long 1D
  pipe; shows the finite-reservoir effect.
- `pipeDamBreakDry` — 1 cell in `y`, `h_R = 0`; the dry stress case.
- `pipeDamBreakNearlyDry` — 1 cell in `y`, `h_R = 1 mm`; the near-dry
  positivity stress case.

**Reference:** Ritter (dry), Stoker (wet); the front speed, `h_2`, `u_2`,
the shock shape.

**Observed:**

- wet 10:1 — Stoker matches well: `c_s = 3.094` vs `3.105 m/s`,
  `h_2 = 0.39625` vs `0.39617 m`, front width ~2.5 cells;
- `h_R = 0` — `FOAM_SIGFPE` at `U = hU/h` (`0/0`) at initialisation; the
  solver cannot even start;
- `h_R = 1 mm` — runs to `t ≈ 16.5 s`, then `FOAM_SIGFPE` in `sqrt`
  because `h` goes negative; `Co_max` grows `0.999 -> 3.24 -> 11.04`;
- the wet long pipe — the plateau decays because the 20 m reservoir is
  drained; a finite-reservoir effect, not a Stoker failure.

**Deliverables:** the failure modes, the baseline numbers, the logs, the
analysis. This is the motivation for Steps 2–3.

## Step 2 — our 1D Cartesian solver

**What:** implement a 1D Cartesian SWE solver that fixes the Step-1
failures.

**Numerics:** HLL (1D) / HLLC (2D) flux, MUSCL reconstruction with a
Barth-Jespersen limiter, positivity floor and dry tolerance, dry-aware
first-order fallback, CFL based on `|u| + c`.

**Validation:** the same 1D pipe cases — wet, dry, nearly dry — against
Ritter/Stoker; front speed, shock shape, `h_min`, positivity, grid
convergence, limiter/CFL sensitivity; a direct comparison with
`shallowWaterFoam`.

**Deliverable:** a working 1D Cartesian solver and a quantitative
improvement table over Step 1.

**Status (done):** implemented in
`solvers/SWEshockFoam/SWEshockFoam.C` as a face loop that also
runs on 2D meshes.  The wet case matches Stoker to `0.001 %` in `h_2`;
the nearly-dry case to `0.09 %`; the dry case stays positive and
conserves mass, although the HLL front is smeared by ~6 %.  See
`SWEshockFoam/README.md` and `SWEshockFoam/all_cases.png`.

## Step 3 — our 1D axisymmetric ("imitation radial") solver

**What:** add the radial geometry to the Step-2 solver and compute the
1D axisymmetric radial dam break.

**Equations (r-weighted conservative form):**

```
d_t (r h)     + d_r (r h u)                 = 0
d_t (r h u)   + d_r (r h u^2 + (g/2) r h^2) = (g/2) h^2
```

with `u(0) = 0` at the symmetry axis and a transmissive outer boundary.
The same HLL/MUSCL/positivity machinery as Step 2.

**Validation:** `h_L = 1`, `h_R = 0.1`, `R0 = 1`, `g = 9.81`; the 1D
radial reference; the Clawpack 2:1 wet benchmark after scaling; the
`h(r,t)`, `u(r,t)`, the front, the shock, the rarefaction.

**Prototype:** `radialDamBreak/radial1d_hllc.py` is the first draft. It
shows the equation and the expected wave structure, but it is
oscillatory near the centre/front; it is superseded by
`SWEshockFoam/radial1D` (see the status below) and is kept only as a
comparison. Do not treat the draft as the final reference.

**Deliverable:** the 1D radial solution and the input for Step 4.

**Status (done):** the `radial` switch in `SWEshockFoam` implements the
`r`-weighted fluxes and the `g/2 h^2` geometric source; `radial1D` runs
to `t = 2 s` with `h_center = 0.013742`, `h_max = 0.238120`,
`r_front = 5.6275`, `u_max = 1.389`.  Mass is conserved and the
lake-at-rest balance is exact.  A small sawtooth remains in `u_r` on the
outer plateau; the Python `radial1d_hllc.py` is much more oscillatory and
is kept only as a prototype.  See `SWEshockFoam/README.md`.

## Step 4 — our 2D radial solver (first validation done)

**What:** the 2D Cartesian SWE solver on the radial initial condition.

**Validation:** axial symmetry; the Step-3 1D radial solution; the
Clawpack 2:1 wet benchmark; the `h`/`u` fields; the front; the shock.

**Status (done, first validation):** `SWEshockFoam/radial2D` uses the
same 400 x 400 mesh and initial condition as `radialDamBreak` with the
2D HLLC/MUSCL solver.  At `t = 2 s` the field is exactly axisymmetric
(`max|Δh| = 0` for both reflections and the transpose), mass is
conserved to `1.2e-10` relative, and the radial profile matches
`radial1D` (`r_front = 5.625` vs `5.6275`; `h_center = 0.01411` vs
`0.01374`).  The 2D mesh is 10x coarser in the radial direction, so the
remaining differences are resolution, not symmetry.  See
`SWEshockFoam/radial_1d_vs_2d.png`.

## Current status

| Step | Status |
|---|---|
| 1. native 1D OpenFOAM | done: wet baseline, dry init FPE, nearly-dry `sqrt` FPE, long-pipe finite-reservoir effect |
| 2. our 1D Cartesian | done: `SWEshockFoam` HLL/MUSCL, wet `h_2` +0.001 %, nearly dry `h_2` −0.09 %, dry positive |
| 3. our 1D axisymmetric | done: `radial1D` r-weighted source, mass conserved, small `u_r` sawtooth |
| 4. our 2D radial | done, first validation: `radial2D` exact symmetry, mass `1.2e-10`, matches `radial1D` |

## Out of scope / side tracks

The following are **not** part of this plan and should not drive the work:

- the nonhydrostatic `Model__` / OpenFOAM Boussinesq, KdV,
  Gardner–Ostrovsky, Sturm–Liouville (old project options);
- the long-time dry self-similar gravity current (a separate asymptotic
  limit; only a short note is kept in the radial analytics);
- finite-reservoir tuning of the long pipe (a Step-1 caveat, not a goal);
- tuning `shallowWaterFoam` itself (it is the baseline, not the target);
- comparing the Clawpack benchmark to our current draft solver (the
  benchmark is a reference in its own right; see the radial analytics).
