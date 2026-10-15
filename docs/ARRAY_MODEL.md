# Finite arrays of elastic cylinders in water

`CylinderExperiment` now supports a list of `centers`, with one common radius
and material. Its default remains one cylinder at the origin, preserving the
single-cylinder benchmark. `cylinder_lattice(columns, rows, pitch)` produces
coordinates of a centered rectangular lattice; scripts may supply explicit
coordinates instead. Every entry is a separate, disjoint solid inclusion.

The water and **all** solid inclusions are assembled into one coupled FEM
system. Multiple scattering and fluid-mediated interactions come from this
solution, not an addition of independent cylinder responses. Each inclusion
has the same inviscid normal-velocity / pressure-traction coupling and plane
strain model described in [CYLINDER_MODEL.md](CYLINDER_MODEL.md). The water
boundary is the outer rectangle only, even for translated inclusions near it.

The current geometry requires at least two nominal mesh sizes of water between
neighboring cylinders, and clearance from the outer boundary. Sources/sensors
must lie in water, outside every inclusion. Geometry and sensor arrays are
copied when remaking or solving a specification. Current TOML schema 3 stores
centers and a typed shape; schema 2 loads as circular inclusions and schema 1
still loads as one cylinder at the origin. See [GALLERY_MODEL.md](GALLERY_MODEL.md)
for the brick shapes.

## Frequency sweep: a deliberately local observable

```julia
using MetaNDT, Gridap, GridapGmsh, Gmsh
spec = CylinderExperiment(radius=0.004,
    centers=cylinder_lattice(2, 2, 0.018))
scan = frequency_sweep(spec, [40e3, 60e3, 80e3];
    sensor=1, time_window=(0.0, spec.duration))
write("frequency_sweep.csv", frequency_sweep_csv(scan))
save_experiment("point.toml", scan.comparisons[2][3].experiment)
```

At every carrier frequency, empty water, the rigid array and the elastic array
are solved with the same source and probes. The returned metric is

`R = sqrt(∫window p_case(t)^2 dt / ∫window p_empty(t)^2 dt)`.

This is a **sensor pressure L2-norm ratio**, not transmitted acoustic energy,
an S-parameter, absorption or a band-gap measurement. It can exceed one and
depends on sensor position, source, measurement window and finite-domain
reflections. A two-cycle pulse has finite bandwidth: the x-axis is its carrier
frequency, not a monochromatic steady-state response.

The integral treats pressure as piecewise linear between integration steps,
including interpolation at window endpoints. A window with less than 1e-6 of
the empty-water full-trace squared norm is rejected instead of returning a
numerically meaningless ratio. A zero source is also rejected for this metric.

One common mesh/time preset is used for the whole sweep: the mesh is refined
to at least 14 nominal water elements per shortest carrier wavelength, and the
time step is shortened if necessary for 25 steps per highest carrier period.
These are presets, not accuracy guarantees; the pulse contains higher frequencies.
Every actual specification is stored in `scan.comparisons[j][case].experiment`.
Only two field frames are retained in sweep runs to limit memory; sensor traces,
energy histories and displacement maxima still use every integration step.

The CSV includes frequency, sensor, time window, reference squared norm,
three ratios, actual mesh/time settings and the maximum energy residual.
Saving an actual point specification allows its traces to be replayed through
`simulate` or `simulate_comparison`. Pluto contains visualization and controls,
not the metric or solver. Array geometry alone does not establish metamaterial
behavior; successful quantitative claims need refinement and boundary studies.
