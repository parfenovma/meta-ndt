# Coated hollow cylinder: the fixed circular test stand

## Contract

`CoatedCylinderExperiment → simulate → CoatedCylinderResult` is independent
of Pluto. `coated_comparison` returns **empty water, bare tube, coated tube**
in that order. `coated_frequency_sweep` changes only frequency and returns
both absolute scattering widths and their ratio. All inputs use SI units.
The experiment is saved as TOML (`physics = "coated_cylinder"`, schema 1).

The default is a sealed, air-filled steel tube in water, outer steel radius
4 mm, inner radius 3.736832623 mm, perfectly bonded silicone thickness 1 mm.
The inner radius matches the static area-averaged steel/air density to water.
The coating makes the composite density slightly different from water again.
Water: 1000 kg/m³, 1480 m/s; air: 1.2 kg/m³, 343 m/s. Steel: 7850 kg/m³,
210 GPa, ν = 0.3. Silicone: 1020 kg/m³, cP = 980 m/s, cS = 30 m/s,
shear loss tangent 0.15. These are explicitly selected model parameters,
not a temperature-calibrated laboratory specification.

The geometry and steel/silicone properties are inspired by
[Rohde et al. (2015)](https://www.nature.com/articles/srep13175).
The authors' silicone shear speed and loss are estimates. This is **not**
an exact reproduction of that paper, nor its numerical gains. Our water
speed is 1480, not the research note's tentative 1482 m/s. The preset thickness
is 1 mm, not a frequency-optimized thickness.

## Equations and coupling

Linear isotropic **plane strain**, infinitely long cylindrical object,
complex harmonic convention `exp(−iωt)`. Fluids obey Helmholtz acoustics.
Both steel and silicone obey `div σ + ρ ω² u = 0` with the 3D isotropic
constitutive law and εzz = 0. The silicone is not a scalar fluid.

For each angular mode n, displacement is
`ur = U(r) cos(n(θ−α))`, `uθ = V(r) sin(n(θ−α))`, where α is incidence angle.
Radial strains are:

```text
εrr = U′,   εθθ = (U+nV)/r
2εrθ = V′ − (V+nU)/r
div u = εrr+εθθ
```

The mixed radial weak form, after cancelling the common angular integral, is:

```text
∫ r [2μ (ε(u):ε(v) − div(u)div(v)/3) − ρ ω² u·v + q div(v)] dr
    = b pair v_r(b) − R pwater v_r(R)
∫ r [div(u) τ − q τ/K] dr = 0
```

P2 displacement is continuous through the bonded steel/silicone interface;
P1 volumetric stress q = K div(u) is independent on each material side.
The deviatoric strain uses the **3D** trace subtraction with εzz=0, not 2D
trace subtraction. Positive bulk modulus K keeps the formulation defined
at λ = 0 / ν = 0. The strain
inner product includes `2 εrθ(u) εrθ(v)`. Five-point Gaussian quadrature.
This addresses nearly incompressible silicone (ν ≈ 0.49953); the independent
elastic Bessel reference and radial refinements check it quantitatively.

Inner air pressure is regular `A Jn(kair r)`. Water pressure is the known
plane wave plus `cn Hn^(1)(kw r) cos(n(θ−α))`. Both fluid interfaces enforce
`∂r p = ρ ω² ur`, `σrr = −p`, `σrθ = 0`. Steel/coating displacement and traction
are continuous. Neither tube nor coating is clamped. Air is not replaced
with water or a traction-free void.

The outer water field is eliminated as an exact boundary impedance in each
modal FEM solve. The regular cavity pressure amplitude is retained explicitly:
there is no division by Jn′ at the cavity's Neumann poles. The outer Hankel
function radiates into **infinite water**:
there is no absorbing rectangle, PML or finite-domain reflection to tune.
There is no full 2D Gridap mesh in this backend. It is a specialized
Fourier/radial FEM backend for concentric objects, not an arbitrary-CAD solver.
Existing Gridap water arrays and connected-solid STG backends remain intact.

μsteel and λsteel follow E,ν. For silicone:
`μ0=ρ cS²`, `K=ρ(cP²−4cS²/3)`, `μ=μ0(1−iη)`, `λ=K−2μ/3`.
K is real, shear loss is passive for our convention. This frequency-domain
loss law is not a causal transient material model. A λ with positive imaginary
part alone does not indicate gain; physical bulk/shear losses determine passivity.

## Results and metric

`coefficients[n+1]` includes incident pressure and the plane-wave factor
`εn i^n`, where ε0=1 and εn=2 otherwise. Exterior fields are calculated from
these coefficients, not interpolated between separately normalized plots.
`coated_pressure(result, point; scattered=true)` returns **complex** scattered
water pressure. The total field adds the same exact incident wave. Points
inside solids/cavity are rejected. `coated_displacement` interpolates actual
P2 solid motion; it is not used as water pressure inside the drawing.

The full scattering width (metres) is outgoing scattered power per unit tube
length divided by incident intensity:

```text
σ = (2/kw) [2|c0|² + Σ(n≥1)|cn|²] / |p0|²
f(θ) = sqrt(2/(πkw)) Σ cn exp(−inπ/2) cos(n(θ−α))
σ = ∫ |f(θ)|² dθ / |p0|²
```

It is radius-independent and includes every retained angular mode. The UI
ratio is `10 log10(σcoated/σbare)`, not single-microphone amplitude and not
a finite near-field pressure integral. Empty water has zero scattering and
is never the denominator. The angular graph shows `|f|²/|p0|²` on one linear
scale; the maps also share a scale, without magnifying the coated case.

Diagnostics include scaled algebraic residual, last-two-mode contribution,
material absorption, extinction and their power balance. Extinction follows
the optical theorem; absorption is calculated independently from the imaginary
part of the mixed elastic form. A small residual is not a mesh-convergence proof.

## Deliberate limits

Fixed concentric layers only; no arbitrary CAD, eccentricity, finite tube ends,
adhesive layer, roughness, flow, bubbles, nonlinear deformation, thermal losses
or experimental fabrication tolerance. The frequency range shown in Pluto is
5–60 kHz, fixed incident plane wave. Phase animation is **not pulse travel**.
Default numerical settings: 12 radial steel elements, 48 coating elements,
angular modes 0–12. Refinement and alternative materials are script inputs,
not a collection of meshing widgets for beginners. An arbitrary high-frequency
or extreme-material experiment needs additional convergence checks.

For actual time propagation see [COATED_PULSE_MODEL.md](COATED_PULSE_MODEL.md).
It uses direct Newmark with causal Kelvin–Voigt coating losses, not this
frequency-independent loss tangent or the harmonic phase animation.

`coated_cavity_resonances(spec, lower, upper)` returns rigid-wall air-cavity
Neumann frequencies as **sampling seeds**, not exact resonances of the
coupled elastic structure. `coated_frequency_sweep(spec, frequencies)` uses
exactly the requested points. With `resolve_cavity=true`, it also samples
±30 Hz around each seed at 0.5 Hz spacing and ±1000 Hz at 100 Hz spacing.
Pluto uses this targeted enrichment and marks computed points; the default
5–60 kHz scan contains 435 points, not only twelve. The actual peaks are
shifted relative to the seeds. This finite scan cannot certify the absence
of other narrower resonances or resolve arbitrary changed materials.

The lossless acoustic-air model has extremely narrow resonances. Near
26.899 and 44.625 kHz, the default coating can **increase** scattering.
The UI supplies a resonant-frequency example and zooms around each cavity
seed. These limitations are part of the experiment, not removed by a
broadband-cloak claim.

This is scattering cancellation, not universal transformation cloaking.
The coating reduces, not erases, the object's scattering. Disabling shear
loss is a control: cancellation survives, so it is not merely absorption.

## Minimal script

```julia
using MetaNDT
s = CoatedCylinderExperiment(frequency=20_000)
cases = coated_comparison(s)
scan = coated_frequency_sweep(s, 5_000:5_000:60_000; resolve_cavity=true)
save_experiment("coated.toml", s)
restored = load_experiment("coated.toml")
write("sweep.csv", coated_sweep_csv(scan))
```

SpecialFunctions is the only new numerical dependency; this fixed circular
backend does not require Gridap or Gmsh. See COATED_CYLINDER_VALIDATION.md.
