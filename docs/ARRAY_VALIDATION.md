# Array stage verification — 2026-10-14

Locally checked with Julia 1.10.12, Gridap 0.20.10, GridapGmsh 0.7.4,
Gmsh.jl 0.3.1 on macOS / Apple Silicon. Two independent reviewers checked
the new stage against the previous working-tree snapshot; no commits/pushes
were part of verification.

## Independent checks

- `test/array_fem.jl`: 95 assertions passed.
- Previous `test/cylinder_fem.jl`: 42 assertions passed.
- Reviewer-only matrix checks: 26 assertions passed; four inclusions have
  separate DOF, no direct inter-body M/K coupling, bidirectional fluid coupling
  on each body, and no volume source applied inside solid bodies.
- Companion `meta-demo/test/array_ui.jl`: 297 assertions passed, including
  288 geometry specifications at control-range corners and actual solver/cache
  tests. These corners were specification-validated, not all accuracy-studied.
- Previous demo tests: 31 passed; detector-only library tests: 18 passed.
- Four targeted checks of post-review UI fixes passed independently.
- `git diff --check`: clean in both repositories.

Commands from sibling `meta-demo`:

```bash
julia --startup-file=no --project=demo ../meta-ndt/test/array_fem.jl
julia --startup-file=no --project=demo ../meta-ndt/test/cylinder_fem.jl
julia --startup-file=no --project=demo test/array_ui.jl
julia --startup-file=no --project=demo test/runtests.jl
```

## Numerical observations

Test array: four cylinders, radius 4 mm, centers at x/y = ±9 mm, 60 kHz,
source at (−45,0) mm, duration 120 μs, h = 2 mm, dt = 0.25 μs.

| Check | Observed result |
|---|---:|
| Near-fluid array versus empty water, relative trace L2 difference | 1.4524% |
| Refinement to h=1.5 mm, dt=0.125 μs, relative trace L2 difference | 4.6849% |
| Maximum displacement over all integration steps | 0.161997 nm |

The paired sweep additionally passed exact replay of a stored point's traces,
reference ratio = 1, invariance of ratios under doubling the source amplitude,
common mesh/time settings and full-step diagnostics despite only two frames.

## Review fixes and browser check

A test initially expected OCC CAD type `Line` after `.msh` reload. Reloaded
entities are discrete; the regression now checks actual outer-node coordinates
on the rectangle's four sides. The mesher identifies those sides geometrically,
so a translated cylinder near an edge is not accidentally tagged as an absorber.

The implementation reviewer found two minor issues: the sweep plot displayed
the requested, not effective, time step, and the beginner guide used the old
button label. Both were corrected and rechecked; no open review findings remain.

The actual Pluto browser session was exercised with six cylinders, then nine.
Unconfirmed control changes retained the old computed sample; confirmation
recomputed the geometry and enabled sweep. A five-point 40–80 kHz sweep
rendered actual results and per-point TOML download buttons. Exported TOML
contained all nine centers and the common actual h=1.3214286 mm; the sweep CSV
contained the correct sensor/window/settings and five rows; sensor CSV contained
601 time samples. Empty/rigid/elastic controls, the scattered-signal selector,
and the 45 μs frame were checked in the browser.

These checks establish the functioning MVP and regression behavior for the
stated presets. They do not establish engineering accuracy over the whole
parameter range, an exact nonreflecting boundary, a band gap, absorption or
metamaterial performance. See [ARRAY_MODEL.md](ARRAY_MODEL.md) for the metric
and limitations. One combined h/dt refinement is a smoke test, not a proven
convergence rate or absolute accuracy bound.
