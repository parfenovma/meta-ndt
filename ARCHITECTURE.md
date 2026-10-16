# Architecture: a scriptable virtual test stand

MetaNDT exposes signal detection, a transient acoustic–elastic experiment and
a connected-solid lens experiment through small interfaces. Clients share the same
interfaces; Pluto owns controls and visualization.

```text
ExperimentSpec / FEM export / hardware acquisition
              │
              ▼
      numerical MetaNDT core
  experiment → FEM → fields + traces + diagnostics
  arrays → detector → statistics
              │
              ▼
   Pluto demo / scripts / future dashboards
```

## Layer 1 — the library

The detector module has no dependency on Pluto, plotting packages, Gmsh or a
particular metamaterial geometry. Its detector functions accept ordinary Julia
arrays. Real and complex traces are supported by the matched-filter projection,
provided that a scientifically appropriate whitening/noise scale is supplied.

The detector deliberately accepts arrays instead of an experiment object.
Units, coordinate axes and acquisition provenance belong to explicit adapters;
the detector must never guess them.

## FEM module

`src/cylinder_experiment.jl` owns `CylinderExperiment`, validation, plain-array
`CylinderResult`, versioned TOML, CSV export and the sequential sweep interface.
`ext/MetaNDTGridapExt.jl` implements `simulate` when Gridap, GridapGmsh and Gmsh
are loaded. Thus users of detection do not load a mesher or FEM implementation.

The experiment specifies SI geometry/materials, a finite source, pressure
sensors, mesh/time parameters and the external boundary. The three initial
cases share one interface: empty water, rigid sample, elastic sample. A sample
is a list of disjoint equal shapes with one material: circles, rectangular
bricks or corrected sinusoidal STG bricks. The single-cylinder benchmark
remains the default. `sample_geometry.jl` owns immutable brick parameters,
outlines and geometric checks shared by meshing and drawing; `lab_examples.jl`
owns six fixed presets with physical questions and matched controls. A frequency sweep
uses the same interface and reports paired sensor pressure norms with explicit
time-window and numerical settings.
The coupled implementation assembles the weak forms with Gridap, including
both interface directions, then advances a monolithic system with Newmark.
The returned diagnostics include work/energy balance and interface skew symmetry.

The mesher writes water/solid/outer/interface physical groups. User-provided
geometry, mixed inclusion materials, alternative wave models, GPU implementations
and automated design search remain future work. They must preserve explicit
units, reproducible settings and numerical validity checks.

See [CYLINDER_MODEL.md](docs/CYLINDER_MODEL.md) for equations and limitations.
See [GALLERY_MODEL.md](docs/GALLERY_MODEL.md) for STG geometry, equal-mass
controls and the separation from the prior solid-only waveguide campaign.

`src/solid_lens.jl` owns the separate `SolidLensExperiment`, harmonic/pulse
plain-array results, reproducible geometry and six fixed presets.
`ext/solid_lens_fem.jl`, inside the same optional Gridap extension, owns genuine
solid connectivity, natural traction-free void boundaries, phase calibration
and the full-wall solve. `meta-demo/src/SolidLensUI.jl` owns bounded caches and
presentation only. Changing the UI does not change the weak form.
See [SOLID_LENS_MODEL.md](docs/SOLID_LENS_MODEL.md). A cell-design approximation
is never returned as if it were the calculated field of a complete wall.

`src/coated_cylinder.jl` owns a separate concentric tube experiment with three
cases and frequency sweeps. It uses mixed P2/P1 radial FEM for both elastic
layers and exact angular acoustic modes for water/air, with radiation into
infinite water. This intentionally small circular backend does not load the
full Gridap/Gmsh extension; it is not an arbitrary geometry substitute for it.
`meta-demo/src/CoatedCylinderUI.jl` owns bounded caches and presentation only.
See [COATED_CYLINDER_MODEL.md](docs/COATED_CYLINDER_MODEL.md).

`src/coated_pulse.jl` adds a direct Newmark concentric water/steel/coating/air
transient. `radial_elasticity.jl` supplies the shared P2/P1 solid forms for both
regimes; scalar radial acoustic FEM is advanced in time, not reconstructed from
harmonic fields. The retained every-step near-field pressure modes support
arbitrary water microphones independently of movie-frame sampling. The
`CoatedPulseExperiment` wrapper saves time/mesh settings and a fixed causal
Kelvin–Voigt loss reference alongside the underlying geometry/materials.
See [COATED_PULSE_MODEL.md](docs/COATED_PULSE_MODEL.md). The existing constant-loss
harmonic result and the finite-window pulse observable must not be conflated.

`src/point_observations.jl` adds backend-free postprocessing: P1
`solid_displacement(result, point)` and `harmonic_trace(A, frequency, times)`.
Together with the existing exact `coated_pressure` evaluator they support point
sensors over retained solutions, without resimulation or raster-pixel sampling.
The trace uses the solver's `exp(-iωt)` convention and is explicitly stationary,
not a synthesized transient arrival. Points in voids are rejected, not zeroed.

For transient STG, `simulate_pulse(spec;retain_history=true)` retains optional
Float32 nodal histories at every integration step (the solver stays Float64).
`solid_trace(result, point)` samples these histories by P1 interpolation,
independently of the sparse display frames. It rejects absent histories or
points in voids rather than fabricating a harmonic trace or rerunning FEM.
This costs 8 × nodes × time samples bytes; callers choose whether to retain it.

## Prepared canonical experiments

`canonical_experiments.jl` owns the tensioned axisymmetric membrane/port
model and a passive Helmholtz transmission-line network. `canonical_spatial.jl`
owns positive-tensor acoustic FEM for the ground cloak. These are explicitly identified model
reductions, not geometric replications of all cited metamaterials.

The client prepares immutable result bundles offline. Opening a lesson,
switching a physical control or sampling the stored measurement lattice does
not call a solver or mesher. The existing connected-solid STG and coated-tube
transients supply the other prepared examples. Rendering and bundle transport
remain outside MetaNDT. See `meta-demo/docs/CANONICAL_LAB.md` for the model
identities, limits, provenance and reproducible preparation command.

## Layer 2 — acquisition adapters

The STG what-if client supplies explicitly synthetic arrays. Future adapters can
provide the same array-level contract from:

- harmonic or transient FEM exports;
- an oscilloscope or DAQ;
- a multichannel transmitter/receiver matrix;
- stored calibration, noise and structural-clutter ensembles.

Adapters own units, channel identity, calibration and provenance. The detector
must not silently guess these.

## Layer 3 — clients

The Pluto notebooks in `meta-demo` are clients. They own the widgets,
explanations and SVG graphics, but none of the detection mathematics. A future
laboratory UI or batch-processing pipeline can call the same library without
Pluto.

## What is intentionally not implemented yet

The MVP adds a specific verified FEM benchmark to the equal-noise kernel.
Before scientific use on hardware, add and validate:

1. typed acquisition metadata and immutable provenance;
2. calibration and unit conversion at the adapter boundary;
3. measured noise/clutter ensembles and whitening;
4. complex/time-domain templates and timing synchronization;
5. repeated defect scenes for POD with confidence intervals;
6. localization and resolution tests, including spatial/temporal sampling;
7. archived raw data and end-to-end reproducibility checks.

This roadmap keeps the educational notebook honest while leaving a clean path
to the research bench.
