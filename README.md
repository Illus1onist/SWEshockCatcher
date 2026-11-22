# SWEshockCatcher

A compact, explicit finite-volume solver for the two-dimensional
shallow-water equations, implemented as an OpenFOAM application, together
with its 1D/2D verification cases and the progress report.

**Author:** Grigory Guliy  
**Course:** Modeling with Open Source Software, Skoltech, Fall 2026  
**Repository:** <https://github.com/Illus1onist/SWEshockCatcher>  
**Progress report:** [report/report.pdf](report/report.pdf)

## What is inside

The repository contains only the SWEshockCatcher project: the solver,
its cases, the analysis scripts, the report, and the reference/baseline
cases used in the project.  Generated OpenFOAM meshes, time directories,
logs, and build artifacts are intentionally not committed; they are
recreated by `blockMesh`, `setFields`, and the `Allrun` scripts.

The solver is a single C++ file
[`solvers/SWEshockFoam/SWEshockFoam.C`](solvers/SWEshockFoam/SWEshockFoam.C).
It reads and writes the usual OpenFOAM fields `h` and `hU`, so
`blockMesh`, `setFields`, and ParaView work without changes.

## Numerical method

- Cell-centred finite volumes on the OpenFOAM polyhedral mesh.
- HLL (1D) / HLLC (2D) Riemann fluxes; HLL fallback at wet/dry faces.
- MUSCL reconstruction in primitive variables with a Barth--Jespersen
  limiter; first-order fallback near drying fronts.
- Forward Euler with a CFL time step; default `CFL = 0.1`.
- Positivity floor, dry tolerance, and a well-balanced 1D radial
  (`r`-weighted) mode with the geometric source `g h^2/2`.
- Exact wall flux `(0, p_w n)` for solid walls and symmetry.

The full method, formulas, and validation are described in the
[progress report](report/report.pdf) and in
[`SWEshockFoam/README.md`](SWEshockFoam/README.md).

## Cases

| Case | Domain | Mesh | Initial state | End time |
|---|---|---|---|---|
| `SWEshockFoam/pipe1D` | `x in [-10,10]` | `2000 x 1 x 1` | wet 10:1 | 1 s |
| `SWEshockFoam/pipe1D_dry` | `x in [-10,10]` | `2000 x 1 x 1` | `h_R=0` | 1 s |
| `SWEshockFoam/pipe1D_nearly` | `x in [-160,160]` | `2000 x 1 x 1` | `h_R=1` mm | 16.5 s |
| `SWEshockFoam/radial1D` | `r in [0,10]` | `2000 x 1 x 1` | `h_L=1`, `h_R=0.1`, `R_0=1` | 2 s |
| `SWEshockFoam/radial2D` | `[-10,10]^2` | `400 x 400 x 1` | `h_L=1`, `h_R=0.1`, `R_0=1` | 2 s |
| `pipeDamBreak*` | native OpenFOAM baseline | see case READMEs | Step 1 | see READMEs |
| `radialDamBreak` | native OpenFOAM baseline | `400 x 400 x 1` | radial dam break | 2 s |

## Key results

- Wet dam break: `h_2 = 0.396180` vs Stoker `0.396175` (`+0.001%`),
  `u_2` within `+0.005%`, shock speed within `-0.07%`.
- Nearly dry (`h_R=1` mm): `h_2` within `-0.09%`, `u_2` within
  `+0.06%`, `c_s` within `-0.12%`.
- Dry bed: the depth remains positive (`h_min = 10^-10` m), mass is
  conserved, and the HLL front is smeared ahead by about `5.8%` in the
  `h > 10^-6` m metric.
- Grid convergence (`N = 128 ... 2048`): observed `L_1` order about `0.9`,
  consistent with a shock-dominated solution.
- 1D radial: `h_centre = 0.013742` m, `r_front = 5.6275` m.
- 2D radial: exactly axisymmetric (symmetry error 0), mass error
  `1.2e-10`, `r_front = 5.625` m vs `5.6275` m in 1D.

## Repository layout

```
solvers/SWEshockFoam/     C++ solver source and Make files
SWEshockFoam/             custom-solver cases, analysis, figures
SWEshockFoam/convergence/ 1D grid-convergence study
report/                   LaTeX progress report and figures
docs/                     radial dam-break analytics
customSolver/             Python 1D SWE prototype and validation
pipeDamBreak*/            Step-1 native OpenFOAM baseline cases
radialDamBreak/           native 2D baseline + Python radial prototype
PLAN.md                   four-step project plan and status
```

## Build and run

A working OpenFOAM 2406 environment is assumed.  Build the solver:

```sh
cd solvers/SWEshockFoam
wmake
```

Run a custom-solver case (this creates the mesh and the initial fields):

```sh
cd SWEshockFoam/pipe1D
./Allrun
```

The other cases are analogous:

```sh
cd SWEshockFoam/pipe1D_dry && ./Allrun
cd SWEshockFoam/pipe1D_nearly && ./Allrun
cd SWEshockFoam/radial1D && ./Allrun
cd SWEshockFoam/radial2D && ./Allrun      # about 10-15 min
```

The grid-convergence study is reproduced with:

```sh
cd SWEshockFoam/convergence
./run_convergence.sh
python3 analyze_convergence.py
```

## Analysis and report

- [`SWEshockFoam/analyze_all.py`](SWEshockFoam/analyze_all.py) вЂ” wet,
  dry, nearly-dry, and radial-1D metrics and `all_cases.png`.
- [`SWEshockFoam/analyze_radial2d.py`](SWEshockFoam/analyze_radial2d.py) вЂ”
  2D radial average, symmetry, mass, and the 1D-vs-2D comparison.
- [`SWEshockFoam/convergence/analyze_convergence.py`](SWEshockFoam/convergence/analyze_convergence.py)
  вЂ” grid-convergence data and plot.
- [`report/report.pdf`](report/report.pdf) вЂ” the progress report;
  [`report/report.tex`](report/report.tex) вЂ” its LaTeX source.

The report follows the structure of the course template:
introduction, mathematical formulation, dimensional analysis, numerical
methodology, computational setup, preliminary results, conclusions,
references, and appendices.  The report PDF is built with
`pdflatex report.tex` (twice, from the `report/` directory).

## Status and limitations

This is a preliminary version of the final project report.  The main
limitations are the small stable CFL number, the diffusive dry front, a
residual oscillation in the radial velocity on the outer plateau, and the
coarse 2D mesh.  The final report will add the Clawpack 2:1 radial
benchmark, a 2D mesh-refinement study, a second-order SSP time
integrator, a smoother well-balanced radial reconstruction, and an
extension to bed slopes, friction, and wetting/drying.

