# Coated cylinder: validation of this implementation

Local validation, 2026-10-15, Julia 1.10.12, Apple ARM64, one BLAS thread.
This document reports **our model**, not a measured physical cloak and not
a reproduction of the gain reported by Rohde et al.

## Reproduce

From the sibling `meta-demo` checkout:

```bash
julia --startup-file=no --project=demo ../meta-ndt/test/coated_fem.jl
julia --startup-file=no --project=demo test/coated_ui.jl
```

116 numerical/API assertions and 33 client assertions passed. Basic detector
and what-if suites also passed. The notebook successfully included as a script
and ran in actual Pluto, with all three fields and both graphs displayed.
The previous cylinder FEM suite (42 assertions) and solid-interface review
regressions (6 assertions) also passed. This is not a rerun of every older
Gridap array/gallery/solid numerical suite.

## Default coating, independently calculated results

Same 4 mm steel outer radius, air-filled wall, 1 mm coating, plane wave:

| Frequency, kHz | σcoated / σbare | 10 log10 ratio, dB |
|---|---:|---:|
| 5 | 0.002640 | −25.78 |
| 10 | 0.003485 | −24.58 |
| 20 | 0.007300 | −21.37 |
| 40 | 0.04069 | −13.90 |
| 60 | 0.1391 | −8.57 |

At 10 kHz bare scattering width is 2.4197172e−5 m; coated width is
8.4321048e−8 m. Empty water is exactly zero scattering by construction.
It is a visual reference, not an independent numerical test.

## Independent strong-form reference and convergence

`test/coated_fem.jl` contains a separate elastic potential/Bessel matching
solution for both annuli. It does **not** call the radial FEM assembler or
its constitutive helper. It matches acoustic air, elastic steel/coating and
outgoing water through displacement and traction conditions, using complex
material parameters. FEM coefficients n=0…5 are compared on five frequencies
for both bare and coated geometries.

For the default coated case, maximum relative modal-vector discrepancy in
this check was about 2.0e−7. Doubling radial elements from 12/48 to 24/96
changed the scattering width by at most 6.4e−8 relative at the five tested
frequencies. Both regression tolerances are 1e−4 to allow cross-platform
roundoff; they do not assert these observed precisions for arbitrary inputs.
Increasing modes 0…12 to 0…16 leaves width within 1e−8 relative.
This is the quantitative locking/convergence check at the actual near-
incompressible silicone parameters, not merely an algebraic residual check.

An independent review found narrow acoustic-air features which the original
5-kHz-spaced curve missed. The final Pluto scan contains 435 points with
0.5-Hz spacing around the three rigid-cavity sampling seeds (26.897242,
44.618332, 55.976119 kHz). The complex air amplitude is retained in both FEM
and Bessel reference, avoiding division by Jn′ at its zeros. Additional
modal/power checks at the first seed and offsets +1.5/+5.5 Hz passed.

At 26.898742 kHz the default coated/bare ratio is ≈1.02694: slight worsening,
not cancellation. The enriched scan's sampled maximum is ≈1.31751 near
44.624832 kHz. A default coated peak around 26.9027 kHz is only about 2.5 Hz
wide (independent review with 0.1 Hz sampling). These are model air-cavity
resonances, not a blanket broadband performance guarantee. The UI marks
computed points, offers local zoom and explicitly warns about finite sampling.
Independent review also checked the second seed's +6.5-Hz feature and the
third seed: Bessel discrepancy ≤3.9e−6, radial-doubling width change ≤4.9e−6,
and relative power-balance error <6e−11 in those additional checks.

## Independent review

Separate Standards and Spec reviewers examined only this stage against a
saved pre-cloak working-tree snapshot; neither changed files. Three P2
findings were addressed: a supported νsteel=0 limit breaking the mixed λ
form (replaced with mixed bulk K), numerical suite missing from standalone
library CI (included in its test entry), and a coarse frequency curve hiding
air-cavity resonances (targeted enrichment, seed zoom and explicit limitations).
Both reviewers independently reran the final numerical/client suites and
reported no remaining actionable implementation findings.

The final client uses a native submit form to read visible settings into one
submitted snapshot, rather than asynchronously combining/confirming child
widgets. A separate real headless Pluto session verified cold initialization,
bond registration and 10→20 kHz propagation to the calculated experiment and
rendered HTML with no errored cells. Numeric option restoration also avoids
the Julia `1.0` / JavaScript `'1'` string-format mismatch. Headless reactivity
is separate from browser event delivery and is not a substitute for visual QA.
The reproducible headless test is `meta-demo/test/coated_pluto.jl` (9 checks),
also included in demo CI.

No finite outer boundary exists, so no PML distance/absorber refinement is
needed for this circular backend. A separate radial **flux** integral on
10, 30 and 60 mm circles matches the width with 1e−10 tolerance. Simply
changing the result's observation radius leaves width unchanged but changes
the sampled pressure. The angular far-field integral independently matches
the modal width by Parseval.

## Passivity and controls

Extinction (optical theorem), scattered power and independent elastic
absorption balance to relative tolerance 1e−5 for bare and coated cases at
all five frequencies. The observed default coated discrepancy was <2e−8.
Material absorption is nonnegative to numerical tolerance; a lossless
coating has zero absorption to 1e−12 m width tolerance.

At 10 kHz, setting silicone shear loss to zero still gives coated/bare
scattering <2% (actual ratio ≈0.00353). Thus the effect is not only absorption.
Changing thickness to 0.5 or 1.5 mm spoils the default cancellation by >5×;
the output is not a stored or multiplied improvement factor.

Other checks: amplitude scales complex fields linearly while normalized width
is unchanged; rotated incidence rotates displacement/pressure; cavity and
solid points are not accepted as water samples; invalid layers/sampling/
wavelength settings fail early; TOML round-trips; reversed sweeps fail;
cache keys preserve geometry/material provenance and have bounded size;
mismatched cases/sweeps cannot be displayed as a paired comparison.
The νsteel=0 physical limit is supported through mixed **bulk K** with
3D deviatoric trace subtraction, and checked against Bessel for both cases.
The numerical suite is included in the standalone library's test entry point.

The three UI thicknesses (0.8/1.0/1.2 mm) build at 60 kHz, and their sweeps
remain independent of the selected-frequency slider. Maps have one linear
amplitude scale; angular curves have one linear intensity scale. Phase controls
animate an already calculated harmonic solution, not an arriving pulse.

## Limits of the evidence

Tests cover the prescribed circular geometry/materials and stated frequencies.
They do not validate arbitrary CAD, all user-selectable material extremes,
3D finite-length tubes, fabrication, laboratory sound speed, adhesives or
experimental detection performance. The specialized radial FEM is not a
full-domain arbitrary-mesh Gridap backend. Published gains are not imported.
