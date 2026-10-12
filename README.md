# MetaNDT.jl

Small Julia building blocks for metamaterial-assisted ultrasonic NDT research.

MetaNDT is the numerical layer between a signal source — FEM export or hardware
acquisition — and a client such as a Pluto notebook, batch pipeline or future
test-stand UI.

The package currently provides:

- deterministic Gaussian receiver-noise models for controlled A/B studies;
- paired measurements with one shared noise realization;
- known-template matched filtering for real and complex traces;
- fixed-threshold decisions with analytical `Pd`/`Pfa` under the stated white
  Gaussian model;
- array-level APIs with no dependency on Pluto, plotting or a geometry solver.

## Install from GitHub

```julia
using Pkg
Pkg.add(url="https://github.com/parfenovma/meta-ndt")
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

The companion teaching demo and frozen FEM dataset live in
[`parfenovma/meta-demo`](https://github.com/parfenovma/meta-demo).

## Test

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. test/runtests.jl
```

## License

Apache-2.0. See [LICENSE](LICENSE).
