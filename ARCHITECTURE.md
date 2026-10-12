# Architecture: from a classroom demo to a research test stand

MetaNDT is deliberately split into three layers.

```text
FEM export or hardware acquisition
              │
              ▼
      numerical MetaNDT core
  arrays → conditioning → detector → statistics
              │
              ▼
   Pluto demo / scripts / future dashboards
```

## Layer 1 — the library

`src/MetaNDT.jl` has no dependency on Pluto, plotting packages, Gmsh or a
particular metamaterial geometry. Its detector functions accept ordinary Julia
arrays. Real and complex traces are supported by the matched-filter projection,
provided that a scientifically appropriate whitening/noise scale is supplied.

The core deliberately accepts arrays instead of a universal experiment object.
Units, coordinate axes and acquisition provenance belong to explicit adapters;
the detector must never guess them.

## Layer 2 — adapters

The first adapter lives in the separate `meta-demo` repository and reads two
tiny CSV files produced by paired full-wave FEM runs. Future adapters can
provide the same array-level contract from:

- harmonic or transient FEM exports;
- an oscilloscope or DAQ;
- a multichannel transmitter/receiver matrix;
- stored calibration, noise and structural-clutter ensembles.

Adapters own units, channel identity, calibration and provenance. The detector
must not silently guess these.

## Layer 3 — clients

The Pluto notebook in `meta-demo` is the first client. It owns the widgets,
explanations and SVG graphics, but none of the detection mathematics. A future
laboratory UI or batch-processing pipeline can call the same library without
Pluto.

## What is intentionally not implemented yet

The current release contains a controlled equal-noise kernel, not a complete
test-stand stack. Before scientific use on hardware, add and validate:

1. typed acquisition metadata and immutable provenance;
2. calibration and unit conversion at the adapter boundary;
3. measured noise/clutter ensembles and whitening;
4. complex/time-domain templates and timing synchronization;
5. repeated defect scenes for POD with confidence intervals;
6. localization and resolution tests, including spatial/temporal sampling;
7. archived raw data and end-to-end reproducibility checks.

This roadmap keeps the educational notebook honest while leaving a clean path
to the research bench.
