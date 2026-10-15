# Direct transient: a coated air-filled tube in water

The common workbench now has a real time-domain tube experiment. This is
direct Newmark integration, not phase rotation, Fourier synthesis of harmonic
fields, or a prescribed echo multiplier. The fixed concentric geometry is the
same as [the harmonic experiment](COATED_CYLINDER_MODEL.md).

## Interface

```julia
using MetaNDT
s = CoatedPulseExperiment(model=CoatedCylinderExperiment(frequency=10_000))
empty, bare, coated = coated_pulse_comparison(s)
mic = coated_trace(coated, (0.015, 0.006)) # times [s], pressure/scattered [Pa]
save_experiment("pulse.toml", s)
```

The core uses plain arrays and SparseArrays, not Pluto or the optional
Gridap/Gmsh extension. `CoatedPulseExperiment` saves BOTH the physical model and
the time/mesh/material-loss settings; loading does not run a solver. Remaking
the duration requires explicitly extending the outer radius if necessary.

## Coupled equations

In the exterior water use scattered velocity potential φw, plus an exact
known incident plane pulse. In the cavity use total potential φa. Pressure is
`p = −ρ ∂tφ`. The incident pressure is

`pinc(x,t) = p0 h(t − (x·eα − source_plane)/cw)`.

`h(t) = sin(2π f t) sin²(πt/T)` for `0<t<T`, zero otherwise, `T=3/f`.
The source plane is before the entire retained region. The solid, air and
scattered field start at rest, and the shown incident field is also zero at t=0.
This is a scattered-field formulation; the analytic background exists outside
the retained region, not an all-space initial rest condition.

Expand in angular cos(n(θ−α)) modes; solid tangential motion uses sin(n(θ−α)).
Scalar acoustic potentials use P1 radial FEM. Both perfectly bonded elastic
layers use the SAME mixed P2 displacement/P1 volumetric-stress forms as the
harmonic backend. Eliminate the static stress using its positive bulk compliance
matrix: `Ksolid = Kdev + B Q⁻¹ Bᵀ`, not an incompressibility-locking P1 substitute.
The radial integration carries weight r; angular energy weights are 2π for n=0
and π otherwise. At the origin, n>0 potentials vanish; n=0 has natural symmetry.

The global system for each mode is

`M d̈ + (Cvisc + Couter + Cinterface) ḋ + K d = f_inc(t)`.

At the outer solid surface of radius R, water normal-velocity coupling is
`+ρw R q u̇r`, and solid traction coupling is `−ρw R vr φ̇w` on the left.
Incident terms on the right are `ρw R q ∂rφinc` and `−R vr pinc`, respectively.
At the inner radius b, air/solid signs reverse: `−ρa b q u̇r` and
`+ρa b vr φ̇a`. Thus the interface matrix is skew; both directions are included.
Internally φ and u are scaled by `F=cw R` and `L=R`, without changing the physics.

Newmark uses β=1/4, γ=1/2 with zero displacement, velocity and acceleration.
The energy/work check is for scattered-water + solid + air energy and the
generalized incident forcing. It is NOT total acoustic energy, which also has
incident/scattered interference terms.

## A causal coating, not constant harmonic loss

`μ0 = ρcoat cS²`, positive bulk K as before. Kelvin–Voigt shear stress is
`2μ0 dev ε + 2μ0 τ dev ε̇`, with `τ = ηref/(2π fref)` and fixed fref=10 kHz.
K stays elastic; the existing 3D deviatoric trace subtraction uses εzz=0.
This is passive and causal. Its loss tangent is `η(f)=2π f τ`, so changing the
pulse carrier does NOT silently retune the material. It matches the old
constant-loss harmonic model at 10 kHz only. This is an explicit model choice,
not a measured broadband silicone specification.

## Domain, histories, observation

Water is locally refined at the interface and extends far beyond the 50×50 mm
picture. The outer boundary has an approximate outgoing impedance. Validation
requires `2(outer_radius−retained_radius)/cw > duration`; an additional check
actually moves the boundary out and compares traces, since implicit FEM does
not have a strict discrete characteristic cone.

Only the near-field pressure histories out to retained_radius=40 mm are kept:
Float32 pressure coefficients at EVERY time step. Newmark arithmetic is Float64.
`coated_trace` P1-interpolates these coefficients in r and evaluates the angular
modes; it adds the exact incident pulse for total pressure. It never reads
pixels, interpolates sparse movie frames or reruns FEM. Points inside the tube,
air or coating are not exterior-water microphone positions; empty water is a
valid reference there. Points beyond the retained region are rejected.

The processed angular observable is `∫p_scattered²dt / ∫p_incident(origin)²dt`
at observation_radius=30 mm. The displayed coated/bare ratio integrates that
observable around the contour. It is finite-window microphone pressure, NOT
stationary far-field scattering width, radiated energy or universal cloaking.
Changing the time window can change it. In particular, the default 460 μs
window at 10 kHz does not contain the complete coating ringdown.

## Verification and limits

`test/coated_pulse.jl`: zero start, finite drive, all three cases, actual air
and solid motion, two-way skew coupling, energy accounting, point/step queries,
TOML round trip and invalid microphone locations.

`test/coated_pulse_validation.jl`: Fourier response versus the independently
Bessel-validated harmonic model at 10 kHz; this needs a SEPARATE 4 ms history,
not the short movie window, because the tube rings after the pulse. Also checks
dt/2, simultaneous radial/angular refinement, and moving the outer boundary.
These checks are targeted, not a proof for every parameter/material or a
validation of manufactured hardware. Narrow air resonances can ring far longer.

Fixed concentric geometry, 2D plane strain, linear isotropic materials, perfect
bonding, lossless acoustic fluids; no CAD, flow, thermal/viscous fluid losses,
finite tube ends or calibrated microphone aperture/noise.
