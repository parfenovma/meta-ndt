# MetaNDT.jl

Small Julia building blocks for metamaterial-assisted ultrasonic NDT research.

MetaNDT is the numerical layer between a signal source — FEM export or hardware
acquisition — and a client such as a Pluto notebook, batch pipeline or future
test-stand UI.

The package currently provides:

- a coupled air/steel/silicone/water concentric cylinder with mixed radial FEM,
  exact harmonic radiation, three controls and full-scattering frequency sweeps;
- a direct time-domain version of that tube: finite pressure pulse, coupled
  elastic motion, retained water-pressure histories and movable microphones;
- a connected solid STG wall with traction-free voids, actual FEM cell-phase
  calibration, harmonic fields and transient pulses;
- a real transient water–elastic-cylinder FEM benchmark through an optional
  Gridap/Gmsh extension, with bidirectional coupling and energy diagnostics;
- versioned experiment descriptions, pressure-trace export and parameter sweeps;
- deterministic Gaussian receiver-noise models for controlled A/B studies;
- paired measurements with one shared noise realization;
- known-template matched filtering for real and complex traces;
- fixed-threshold decisions with analytical `Pd`/`Pfa` under the stated white
  Gaussian model;
- array-level APIs with no dependency on Pluto, plotting or a geometry solver.

## Install from GitHub

### Run the companion Pluto application

The interactive workbench needs **both repositories cloned together into
sibling directories**. Run these commands from one common parent directory:

```sh
git clone https://github.com/parfenovma/meta-ndt.git
git clone https://github.com/parfenovma/meta-demo.git
cd meta-demo
julia +1.10 --startup-file=no --project=. start.jl
```

Do not clone one repository inside the other. Downloading two ZIP archives
is also supported. The [macOS / Windows installation guide](https://github.com/parfenovma/meta-demo#readme)
explains Julia installation, folder names and updates. Update both repositories
together. The [physical practicum and PDF](https://github.com/parfenovma/meta-demo/tree/main/docs/practicum)
live in `meta-demo`.

### Use only the numerical library

For scripts that do not use the application, install this package directly:

```julia
using Pkg
Pkg.add(url="https://github.com/parfenovma/meta-ndt")
```

## Hollow tube and scattering cancellation

```julia
using MetaNDT
spec = CoatedCylinderExperiment(frequency=10_000)
empty, bare, coated = coated_comparison(spec)
scan = coated_frequency_sweep(spec, 5_000:5_000:60_000; resolve_cavity=true)
save_experiment("coated.toml", spec)
write("coated_frequency.csv", coated_sweep_csv(scan))
```

This fixed circular backend needs neither Gridap nor Gmsh. Both solid layers
are elastic, the cavity contains acoustic air, and the exterior is infinite
water. The metric is genuine 2D scattering width, not a single receiver gain.
See [the physical model](docs/COATED_CYLINDER_MODEL.md) and
[verification](docs/COATED_CYLINDER_VALIDATION.md). It is a source-inspired
scattering-cancellation model, not a general-purpose transformation cloak.

For a finite pulse rather than a harmonic response:

```julia
pulse = CoatedPulseExperiment(model=spec)
water, bare, coated = coated_pulse_comparison(pulse)
microphone = coated_trace(coated, (0.015, 0.006)) # metres, every time step
save_experiment("coated_pulse.toml", pulse)
```

The transient solver advances the coupled water, solids and cavity in time;
it does not animate a phasor. Coating loss is causal Kelvin–Voigt shear damping,
matched to the harmonic model at a fixed 10 kHz reference. Its default finite
observation window does not include all late ringing. See
[COATED_PULSE_MODEL.md](docs/COATED_PULSE_MODEL.md) for equations and validation.

## Cylinder in water

Install `Gridap`, `GridapGmsh` and `Gmsh` in the calling environment to enable
the optional FEM extension, then:

```julia
using MetaNDT, Gridap, GridapGmsh, Gmsh
spec = CylinderExperiment()
water, rigid, elastic = simulate_comparison(spec)
save_experiment("experiment.toml", spec)
write("elastic_sensors.csv", sensor_csv(elastic))
elastic.diagnostics
```

Finite arrays use the same coupled solver, not a sum of single-cylinder traces:

```julia
spec = CylinderExperiment(radius=0.004, centers=cylinder_lattice(2, 3, 0.018))
results = simulate_comparison(spec)
scan = frequency_sweep(spec, [40e3, 60e3, 80e3]; sensor=1)
write("sweep.csv", frequency_sweep_csv(scan))
```

The sweep metric is a paired **sensor pressure norm ratio**, not energy
transmission. See [ARRAY_MODEL.md](docs/ARRAY_MODEL.md) for its definition,
reproduction, geometry restrictions and numerical limitations.

The six fixed laboratory examples include the corrected sinusoidal STG brick:

```julia
spec = example_experiment(:stg_single; frequency=100e3)
stg = simulate(spec)
straight = simulate(straight_reference(spec))  # same length and mass
scan = frequency_sweep(spec, [60e3,100e3,140e3]; include_straight_reference=true)
```

See [GALLERY_MODEL.md](docs/GALLERY_MODEL.md). Brick outlines drive both actual
Gmsh geometry and the client drawing; the display is not a substitute geometry.

This solves water velocity potential and plane-strain solid displacement on
a conforming P1 mesh, coupled through normal velocity and pressure traction.
`p = -ρ ∂tφ`. Newmark average acceleration advances the actual transient fields.
The finite outer domain has a first-order absorbing boundary.
See [the physical model and verification](docs/CYLINDER_MODEL.md).

## STG lens in solid aluminium

This is a different physical experiment: the STG ligaments are joined to solid
continua at their ends, not immersed in water. No legacy amplitude gain is used.

```julia
using MetaNDT, Gridap, GridapGmsh, Gmsh
design = design_stg_lens(SolidLensExperiment())
spec = solid_preset(design, :focused)
field = simulate(spec)
pulse = simulate_pulse(spec; retain_history=true)
sensor = solid_trace(pulse, (0.040, 0.001)) # seconds/metres, every time step
save_experiment("lens.toml", spec)
```

Six fixed examples compare bulk, equal-mass straight ligaments, uniform STG,
the calibrated wall, permuted cells and a missing ligament. See
[SOLID_LENS_MODEL.md](docs/SOLID_LENS_MODEL.md) for the equations, accepted-power
normalization and limits of the focusing claim. FWHM is not NDT resolution.

For local development, place `meta-ndt` and `meta-demo` in sibling directories;
the latter's `demo/` environment loads the local package. From `meta-demo`:

```bash
julia --startup-file=no --project=demo ../meta-ndt/test/cylinder_fem.jl
julia --startup-file=no --project=demo ../meta-ndt/examples/cylinder_comparison.jl
```

## Minimal example

```julia
using MetaNDT

primary = [1.0, 2.0, 1.0]
control = [0.7, 1.2, 0.7]
noise = GaussianReceiverNoise(0.8)

measurement = simulate_paired_measurement(
    primary,
    control;
    noise_model=noise,
    signal_present=true,
    measurement_index=3,
)

decisions = detect_pair(primary, control, measurement)
```

The caller owns units, calibration, channel identity and whitening. The library
does not infer any of them.

## Research path

The public API is intentionally geometry-neutral. Planned test-stand layers
include acquisition provenance, calibration adapters, measured clutter
ensembles, whitening, POD confidence intervals, localization and resolution
metrics. See [ARCHITECTURE.md](ARCHITECTURE.md).

The companion interactive cylinder benchmark and STG what-if demo live in
[`parfenovma/meta-demo`](https://github.com/parfenovma/meta-demo).

## Test

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. test/runtests.jl
```

## License

Apache-2.0. See [LICENSE](LICENSE).
