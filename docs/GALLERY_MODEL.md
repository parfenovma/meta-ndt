# Six laboratory examples, including the real STG brick

The gallery is a client of the same transient water/solid solver. No wave
equation is implemented in Pluto and no frozen echo multiplier is involved.
`lab_examples()` returns six named questions; `example_experiment(id)` returns
a reproducible specification with fixed geometry and numerical presets.

| Example | Physical question | Additional control |
|---|---|---|
| One cylinder | What changes when the object can move? | Empty water and rigid object |
| Nine-cylinder lattice | How do neighboring inclusions interact? | One-cylinder example |
| Missing central cell | How does a structural defect change the field? | Complete lattice, same positions |
| Straight brick | What is the mechanical/shape reference for STG? | Same length, area and material as STG |
| One sinusoidal STG | How does meandering change motion and scattering? | Straight brick of equal mass |
| Three STG bricks | How do immersed bricks interact? | Three straight bricks of equal mass and placement |

These questions are not labels of proven band gaps, negative parameters,
cloaking or focusing. The current materials are lossless; dissipation at the
outer boundary is not metamaterial absorption. A finite pulse and a pressure
minimum at one probe cannot establish those stronger properties.

## Geometry, not just an illustration

`CylinderExperiment` preserves its previous name for compatibility. Its
`shape=nothing` means circles with the existing `radius`. The two new immutable
shape descriptions are `RectangularBrick` and `SinusoidalBrick`. Locations stay
in `centers`. The same `specimen_outlines` generate the brick CAD polygons,
source/sensor distance checks, and client drawings. The FEM weak forms are
unchanged; the whole water domain and all solid inclusions are solved jointly.

The STG formula matches the project's corrected `SinusoidalProfile`, not its
legacy opposing-constriction branch. Default dimensions in metres are:
L=.017, H=.007, A=.0025, N=2, end margin m=.0004.

The bottom active span is `S=(L−2m)·2N/(2N+1)=.01296`, with half-period
shift `δ=S/(2N)=.00324`. Bottom indentation is
`db(x)=A(1−cos(2πN(x−m)/S))/2` within `[m,m+S]`, zero outside.
Top indentation starts at `m+δ`; over the common interval it is `dt=A−db`.
After the bottom active span ends, its first half-period is replayed to taper
the upper wall back to the straight right end. Thus

`ytop−ybottom = H−db−dt = H−A = .0045`

throughout the common span, with full-height end transitions. This is a
meandering brick, not symmetric narrowing with minimum `H−2A`.

The mesher uses a polygonal approximation with 80 bottom profile samples.
The count is divisible by 2N, so the half-period shift aligns the common
polyline knots. Nearly identical floating-point knot positions are merged
before CAD construction to avoid zero-length edges. Shape accuracy, as well
as the FE mesh, must be refined for quantitative studies.

STG area is `LH−AS=86.6 mm²`. The straight reference is L by 5.09411765 mm,
the same area and density, hence the same mass per axial length. It uses the
same source, probes, materials, time step and nominal mesh size; geometry
alone is varied. It is not the old rectangular waveguide experiment.

## Persistence and measurements

TOML schema 3 stores a typed `[shape]` description as well as coordinates.
Schemas 1/2 still load as circles. Cached results use the full serialized
specification, not the display name. Snapshot copying of the sensor and
placement lists is preserved.

`frequency_sweep(...; include_straight_reference=true)` adds a fourth actual
run at every frequency for an STG specimen. Ratios remain local pressure
L2 norms normalized to the paired empty-water trace, not power transmission.
The CSV gains a `straight_ratio` column; actual point specifications preserve
the common mesh/time settings. See [ARRAY_MODEL.md](ARRAY_MODEL.md).

```julia
using MetaNDT, Gridap, GridapGmsh, Gmsh
s = example_experiment(:stg_single; frequency=100e3)
r = simulate(s)
control = simulate(straight_reference(s))
scan = frequency_sweep(s, [60e3,100e3,140e3]; include_straight_reference=true)
save_experiment("stg.toml", r.experiment)
write("scan.csv", frequency_sweep_csv(scan))
```

This is immersed, linear 2D plane-strain elasticity coupled to water, not the
previous solid-only STG port/waveguide study. Its archived working-band or
amplitude-ratio results are not evidence for these immersed examples.
