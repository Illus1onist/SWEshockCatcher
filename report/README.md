# Progress report

`report.pdf` is the compiled progress report; `report.tex` is its LaTeX
source.  It follows the structure of the course template (provided
separately by the course): introduction, mathematical formulation,
dimensional analysis, numerical methodology, computational setup, results,
conclusions, references, and appendices.

## Build

The report includes a source listing from
`../solvers/SWEshockFoam/SWEshockFoam.C`, so compile it from this directory:

```sh
pdflatex report.tex
pdflatex report.tex
```

The second pass resolves the cross-references and the bibliography.
The `Student Name` placeholder is defined once near the top of `report.tex`:

```tex
\newcommand{\studentname}{Student Name} % replace with your name
```

Change it before submission.

## Contents

- `report.tex`, `report.pdf` — the report.
- `figures/setup.png` — the 1D pipe and 2D radial setups.
- `figures/all_cases.png` — wet/dry/nearly-dry/radial-1D verification.
- `figures/convergence.png` — 1D grid-convergence study.
- `figures/radial_evolution.png` — 1D radial solution at three times.
- `figures/radial_1d_vs_2d.png` — 1D axisymmetric vs 2D radial.
- `make_setup_figure.py`, `make_radial_evolution.py` — figure generators.

The numerical data behind the figures are in `../SWEshockFoam`
(`all_summary.json`, `radial2d_summary.json`,
`convergence/convergence.json`).
