# Radial dam break: what is used for comparison, and what analytics exists

This note answers two questions:

1. Which initial states (depths, wet/dry) are used in the 2D radial
   simulation we planned to compare against?
2. How hard is it to write the radial "formulas", and what exactly can be
   written in closed form?

Short answer: the **wet, finite-radius** radial dam break has **no simple
Ritter-like closed form**. What *can* be written is (a) the 1D radial
characteristic form and the method of characteristics, (b) the
self-similar **dry** axisymmetric gravity current, and (c) the exact
Riemann problem at the shock. The wet finite-radius case is a 1D radial
**numerical reference**, not a formula. The Clawpack benchmark is a
**numerical** benchmark, not an analytic solution.

---

## 1. The comparison cases

### 1.1 Clawpack `radial_dam_break` benchmark

The Clawpack example we referenced uses:

| Parameter | Value |
|---|---|
| Domain | `[-2.5, 2.5] x [-2.5, 2.5]` |
| Grid | `150 x 150` |
| Inner depth | `h_in = 2.0` |
| Outer depth | `h_out = 1.0` |
| Dam radius | `R0 = 0.5` |
| Gravity | `g = 1.0` (dimensionless) |
| Final time | `tfinal = 2.5` |
| Riemann solver | Roe or HLLE |
| Limiter | MC (TVD) |
| BCs | top/right: wall; bottom/left: outflow |

So the Clawpack case is **wet**, with a **2:1** depth ratio (`2.0` inside,
`1.0` outside). It is **not** a dry-bed case. The shock is relatively
weak; there is no `h -> 0`; no positivity issue. The Clawpack example is a
**numerical benchmark** (a reference solution computed by a well-tested
code), not an analytic formula.

### 1.2 Our `radialDamBreak` case

| Parameter | Value |
|---|---|
| Domain | `[-10, 10] x [-10, 10] m` |
| Grid | `400 x 400 x 1` |
| Inner depth | `h_L = 1.0 m` |
| Outer depth | `h_R = 0.1 m` |
| Dam radius | `R0 = 1.0 m` |
| Gravity | `g = 9.81 m/s^2` |
| Final time | `t = 2 s` |
| BCs | closed tank (`hU` slip) |

Our case is also **wet**, but with a **10:1** depth ratio. It is a stronger
bore than the Clawpack case. The 1D radial reference must be computed for
**these** numbers, not for the Clawpack numbers.

### 1.3 Dry radial case (stress test)

`h_R = 0` is a separate stress case. It has no ambient water; the front is
a **dry front**; the 1D radial problem approaches the **self-similar
axisymmetric gravity current** at long times (Section 4). `shallowWaterFoam`
cannot even initialise `h = 0` (`U = hU/h` gives `0/0`), so the dry radial
case is a target for the custom positivity-preserving solver, not for the
native PISO run.

---

## 2. 1D radial shallow-water equations

### 2.1 Equations

Axisymmetric shallow water, no rotation, flat bed:

```
h_t + (1/r) (r h u)_r = 0
(h u)_t + (1/r) (r h u^2)_r + g h h_r = 0
```

Non-conservative form:

```
h_t + u h_r + h u_r + h u / r = 0
u_t + u u_r + g h_r = 0
```

With `c = sqrt(g h)`:

```
c_t + u c_r + (c/2) u_r + (c u) / (2 r) = 0
u_t + u u_r + 2 c c_r = 0
```

### 2.2 Characteristics and Riemann invariants

Define

```
J_+ = u + 2c,     J_- = u - 2c.
```

Along the characteristics

```
C_+:  dr/dt = u + c,     dJ_+/dt = - c u / r
C_-:  dr/dt = u - c,     dJ_-/dt = + c u / r
```

In the Cartesian limit `1/r = 0`, `J_±` are constant and the Ritter fan
follows. In the radial case the **geometric source** `∓ c u / r` changes
the invariants along the characteristics: the characteristics are not
straight, and there is no simple `r/t` similarity for the finite-radius
problem.

### 2.3 Why there is no Ritter-like closed form

The Cartesian Ritter solution works because:

- there is no length scale except the initial discontinuity;
- the solution depends on `x/t`;
- `J_±` are constant;
- the rarefaction is a straight fan;
- the shock speed is constant.

The radial case has:

- a length scale `R0`;
- the `1/r` source, so `J_±` vary;
- curved characteristics;
- a shock whose strength and speed **decay** with `r`;
- the ambient depth `h_R`, which adds a second state and a bore.

There is no closed-form `h(r,t)`, `u(r,t)` for the wet finite-radius case.
The practical reference is a **1D radial numerical solution** (method of
characteristics or a high-resolution finite-volume solver).

---

## 3. 1D radial method of characteristics (the practical reference)

### 3.1 Characteristic ODEs

Given `h(r,t)`, `u(r,t)` on a previous time level, integrate:

```
dr/dt = u ± c
dJ_±/dt = ∓ c u / r
```

with `J_+ = u + 2c`, `J_- = u - 2c`. At each new point, recover

```
u = (J_+ + J_-)/2
c = (J_+ - J_-)/4
h = c^2/g.
```

The `r = 0` axis is a symmetry axis: `u(0,t) = 0`, `h_r(0,t) = 0`. The
`C_-` characteristic is degenerate there; use the symmetry condition.

### 3.2 Shock conditions

For a shock at `r_s(t)` with speed `s = dr_s/dt`, the Rankine-Hugoniot
conditions in the radial frame are

```
s (h_- - h_+) = h_- u_- - h_+ u_+
s (h_- u_- - h_+ u_+) = h_- u_-^2 + (g/2) h_-^2
                         - [h_+ u_+^2 + (g/2) h_+^2]
```

where `-` is behind and `+` ahead. These are the same algebraic relations
as in the Cartesian case; the radial geometry enters through the
characteristics and through the varying `h_+` ahead of the shock. The
shock speed is **not constant**: as `r` grows, the shock weakens.

### 3.3 Rarefaction

In a rarefaction, one family of characteristics is a fan. In the radial
case the invariants are not constant along the fan (because of
`∓ c u / r`), so the fan is not a straight-line solution. It can be
integrated as a set of `C_+` or `C_-` curves with the source term.

### 3.4 Algorithm and boundary conditions

A practical 1D radial reference:

1. Initialise `h(r,0)`, `u(r,0)` (the dam break step).
2. Compute `c`, `J_±`.
3. Advance the `C_±` characteristics with a small `dt`; integrate the
   source `∓ c u / r`.
4. Detect and fit the shock (RH); apply the entropy condition.
5. Apply `u(0)=0`, `h_r(0)=0`; use a transmissive or fixed reservoir
   boundary at `r = R_max`.
6. Output `h(r,t)`, `u(r,t)` on a regular `r` grid for comparison with
   the 2D OpenFOAM fields.

This is the "not difficult" part: it is a 1D code, not a new closure. It is
also the reference that should be used for our `h_L = 1`, `h_R = 0.1`,
`R0 = 1`, `g = 9.81` case.

---

## 4. Dry limit (out of scope for the main plan)

For `h_R = 0` the long-time axisymmetric gravity current is self-similar: `R ~ t^{1/2}`, `h ~ t^{-1}`, `u ~ t^{-1/2}`; the profiles solve a second-order ODE system. This is a separate asymptotic limit, not needed for the wet radial case and not part of the four-step plan.

## 5. The wet radial case: what can and cannot be written

### 5.1 No simple closed form

For `h_R > 0` the flow is a **rarefaction + bore**. The shock is a bore
into the ambient; its strength decays as `r` grows. The initial radius
`R0` is a length scale, and the `1/r` source makes the Riemann invariants
vary. There is no closed-form `h(r,t)`, `u(r,t)` analogous to Ritter or
Stoker.

### 5.2 What the 1D radial reference gives

The 1D radial method of characteristics / finite-volume solution gives the
reference profiles `h(r,t)`, `u(r,t)`, the shock position `r_s(t)`, the
shock speed `dr_s/dt`, the depth jump `h_2(r_s) - h_R`, and the
rarefaction. This is the practical "analytic-like" reference for the 2D
simulation.

### 5.3 What the Clawpack benchmark gives

The Clawpack `radial_dam_break` is a **numerical** benchmark for the wet
2:1 case (`h_in=2`, `h_out=1`, `R0=0.5`, `g=1`, `tfinal=2.5`). It is not an
analytic solution; it is a well-tested reference. It is useful as a
cross-check of the 2D solver and the 1D radial reference after scaling the
parameters.

---

## 6. Where this fits the project plan

This document supports Step 3 (1D axisymmetric) and Step 4 (2D radial) of
[PLAN.md](../PLAN.md):

1. Step 1 — native 1D OpenFOAM `shallowWaterFoam` on the pipe: baseline
   and failure modes (dry `U=hU/h` FPE; nearly dry negative `h`/`sqrt`
   FPE).
2. Step 2 — our 1D Cartesian HLL/HLLC/MUSCL/positivity solver
   (`SWEshockFoam`).
3. Step 3 — our 1D axisymmetric solver (`r`-weighted SWE, geometric
   source); this document's equations and the 1D radial reference
   (`SWEshockFoam/radial1D`).
4. Step 4 — our 2D radial solver (first validation:
   `SWEshockFoam/radial2D`).

---

## 7. Summary

- **Clawpack benchmark:** wet, `h_in=2`, `h_out=1`, `R0=0.5`, `g=1`,
  `150 x 150`, `tfinal=2.5`; a **2:1 wet** numerical benchmark, not a dry
  case and not an analytic formula.
- **Our `radialDamBreak`:** wet, `h_L=1`, `h_R=0.1`, `R0=1`, `g=9.81`,
  `400 x 400`, `t=2`; a **10:1 wet** case; needs its own 1D radial
  reference.
- **Radial formulas:** the 1D radial **characteristic form** can be
  written; the **wet finite-radius** radial dam break has **no simple
  closed form**. The practical reference is a 1D radial numerical
  solution. The dry self-similar limit is a separate asymptotic and is
  out of scope for the main plan.
- **The "not difficult" part** is the 1D radial reference (a 1D code),
  not a Ritter-like closed formula for the wet radial case.
