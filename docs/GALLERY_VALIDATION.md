# Six-example / immersed STG verification — 2026-10-15

Local environment: Julia 1.10.12, Gridap 0.20.10, GridapGmsh 0.7.4,
Gmsh.jl 0.3.1, macOS / Apple Silicon. The stage was reviewed against a
working-tree snapshot taken before the gallery changes; older uncommitted MVP
work was preserved. No commit, push or remote CI run was performed.

## Independent review

The implementation reviewer checked repository standards, gallery controls,
four-case rendering, cache keys, schema compatibility and exports. The physics
reviewer checked the requested geometry and actual FEM domains/coupling.
No open blocking findings remain in those reviews.

| Suite | Passed assertions |
|---|---:|
| Library `test/gallery_fem.jl` | 229 |
| Previous cylinder FEM suite | 42 |
| Previous array FEM suite | 95 |
| Demo `test/gallery_ui.jl` | 20 |
| Updated array client suite | 307 |
| Previous demo suite | 31 |
| Detector-only library suite, run by the main agent | 18 |
| Total retained suites | 742 |

Additional reviewer-only checks, not included in that total:

- 3208 evaluations against the original corrected `depth_at_x` and
  `top_depth_at_x` functions, for N=1,2,3,6 including the end transitions.
- Six actual-mesh checks: 2582 nondegenerate STG solid triangles have area
  86.60000000000002 mm²; solid/water centroids lie in their respective domains;
  solid plus water areas reproduce the complete rectangle.
- Two interface-normal checks: closed-boundary normal integral is approximately
  zero and `∫ x·n ds = 2A`, with the normal outward from the solid.
- 24 rendered-description checks across all six notebook cards.

These confirm the meandering H−A geometry and equal-mass straight control,
not symmetric H−2A constrictions or a substituted striped medium. The six
FEM presets, including the three-body STG row, were actually solved.

Run from sibling `meta-demo`:

```bash
julia --startup-file=no --project=demo ../meta-ndt/test/gallery_fem.jl
julia --startup-file=no --project=demo ../meta-ndt/test/cylinder_fem.jl
julia --startup-file=no --project=demo ../meta-ndt/test/array_fem.jl
julia --startup-file=no --project=demo test/gallery_ui.jl
julia --startup-file=no --project=demo test/array_ui.jl
julia --startup-file=no --project=demo test/runtests.jl
```

Run the detector-only suite from `meta-ndt` with `--project=.`.

## Numerical observations, not performance claims

For the default immersed STG, 100 kHz, h=1 mm, dt=0.25 μs, duration 150 μs:

| Observation | Result |
|---|---:|
| Near-fluid limit versus empty water, relative trace L2 difference | 6.17135% |
| Joint refinement to h=0.75 mm / dt=0.125 μs, relative trace L2 difference | 1.13678% |
| STG versus equal-mass straight control, relative difference of all three traces | 33.87595% |
| Maximum STG displacement over all integration steps | 0.160175 nm |
| Maximum straight-control displacement | 0.063499 nm |

The paired 60/100/140 kHz sweep solves four cases per point. At sensor A the
STG ratios were 1.12204 / 0.99231 / 0.85889; the straight-control ratios were
1.01109 / 0.95297 / 0.82896. These are pressure L2 norms divided by the
corresponding empty-water norm, not transmitted powers or useful echo gains.
Stored point specifications replay their actual geometry and numerical preset.

The transparent limit still has a several-percent discretization discrepancy.
One combined h/dt refinement does not establish an absolute accuracy bound,
an asymptotic convergence rate or accuracy for every preset/frequency. Profile
polygon resolution and outer-boundary placement also need refinement studies
before quantitative scientific conclusions. No focusing, cloak, band gap,
absorption or negative effective parameters are established by these checks.

## Review fixes

The old array-client test checked the literal `data.centers` renderer text.
After shared outline rendering replaced circle-specific drawing, that test was
updated to inspect the actual embedded centers and both circular boundaries in
all three cases; the 307 assertions then passed independently.

The Julia Markdown string macro misinterpreted two property interpolations
in the new notebook's specimen description as inline mathematics. The card now
uses `Markdown.parse` on an ordinary interpolated string, with rendered-text
regressions. Download filenames are built outside the Markdown macro.

See [GALLERY_MODEL.md](GALLERY_MODEL.md) for the geometry and measurement
contract, and [ARRAY_MODEL.md](ARRAY_MODEL.md) for the local sweep observable.

## Live Pluto check

A fresh gallery session was exercised with the default one-body STG and the
three-body STG row. An unconfirmed selection kept the previously computed
specimen; confirmation updated the field, deformation inset, signals and
downloads. The fourth button switches to the separately solved equal-mass
straight geometry. The shared outline renderer showed the actual meandering
profiles in both water-domain and solid-motion views.

The default 60/100/140 kHz sweep displayed all four curves and a three-row
table. Its exported CSV had 12 columns, including `straight_ratio`, and the
same actual h=0.75510204 mm / dt=0.25 μs at each point. The 100 kHz point TOML
contained `[shape] kind="stg"`, N=2 and 80 samples, not the straight reference;
its frame count of two is the documented sweep storage optimization.

One older Pluto session stalled while displaying the new output; both server
and worker were idle in one-second diagnostic samples. A standalone full
notebook run succeeded, and warm HTML construction took 0.881 s for a 52.1 MB
payload. A fresh Pluto process displayed the same specimen, then completed
the row and sweep. The exact cause of that old-session stall was not isolated;
it is not evidence of a slow FEM solve or a universal UI timing guarantee.
