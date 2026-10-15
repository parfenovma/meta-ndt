# Transient water–elastic-cylinder benchmark

This is a real monolithic, conforming P1 finite-element calculation, assembled
with Gridap on a tagged Gmsh mesh. It is a linear 2D cross-section of an infinite
cylinder, with **plane strain**, not a finite 3D specimen. SI units are used.
Water defaults: density 1000 kg/m³, sound speed 1480 m/s. Aluminium defaults:
density 2700 kg/m³, Young modulus 70 GPa, Poisson ratio 0.33.

## Equations and signs

The water unknown is velocity potential φ, with fluid velocity ∇φ and pressure
`p = −ρf ∂tφ`. The solid unknown is displacement u. On the interface, ns points
out of the solid into water; nf = −ns. The inviscid interface conditions are:

```
∂nf φ = −ns ⋅ ∂t u             normal velocity continuity
σ(u) ns = −p ns                pressure traction; zero tangential traction
```

The weak forms, with scalar test q and solid vector test v, are:

```
∫water (ρf/cf²) q φ̈ + ρf ∇q⋅∇φ
    + ∫interface ρf q ns⋅u̇ + ∫outer (ρf/cf) q φ̇ = ∫water −ρf q Q(t) g(x)

∫solid ρs v⋅ü + ε(v):σ(u) − ∫interface ρf (v⋅ns) φ̇ = 0
```

Here σ = λ tr(ε) I + 2με, with the 3D Lamé constants and εzz = 0.
The two interface terms are transposes with opposite signs. Their instantaneous
power cancels; neither side receives a one-way imposed or synthetic response.
The assembled system is `M ẍ + (Couter + Cinterface) ẋ + K x = f(t)`, with
symmetric mass/stiffness matrices and a skew-symmetric interface matrix.
The implementation applies a congruent variable scaling `φ = cf L xφ`,
`u = L xu`, with L the domain width, to improve conditioning without changing
the physical equations.

The source is a normalized Gaussian of fixed physical width, driven by a
Hann-windowed sinusoidal volume-velocity pulse Q(t), in m²/s per axial metre.
It is a finite spatial source, not an unresolved mesh-node delta. Pressure
sensors interpolate the FE solution in water; they do not alter the solution.

## Time, energy and external boundary

Newmark average acceleration (`β = 1/4`, `γ = 1/2`) uses one sparse LU
factorization per experiment. All frames are sampled from the actual transient
solution. A common pressure scale must be used when comparing the three cases.
The peak-displacement diagnostic is accumulated at every integration step,
independently of how many frames are retained for visualization.

The outer rectangle uses a **first-order local absorbing boundary**. It is
exact only for normal plane-wave incidence; oblique waves can reflect. The
finite domain and this boundary are part of the experiment, not a perfect PML.
Use a larger domain or early-time signals to check sensitivity before making
scientific conclusions. The cylinder is traction-free except for fluid loading;
its rigid translations are not artificially clamped.

Each result carries kinetic+strain/acoustic energy, source work and accumulated
outgoing boundary dissipation (J per axial metre). The discrete work uses
midpoint velocity and averaged endpoint force, consistent with Newmark. Check:
`E(t) + Wout(t) − Win(t) ≈ 0`. This checks integration and passive exchange;
it does not by itself prove spatial accuracy.

## Verification

`test/cylinder_fem.jl` exercises the public solver: zero forcing, nonzero solid
motion and feedback on pressure, interface power cancellation, energy balance,
the near-fluid solid limit, and mesh/time refinement. The near-fluid case uses
ρs = ρf, ν → 0.5 and E = 3ρf cf²(1−2ν), giving negligible shear modulus and
approximately matching compressional speed. It should approach empty water.

P1 elements and the default mesh are a responsive MVP preset. A quantitative
study must refine mesh and time step, check the outer boundary and source
width, and account for the 2D assumption. The displayed refinement test reports
sensor-trace differences rather than inferring accuracy from energy alone.

Useful primary references:

- [COMSOL's acoustic–structure interface conditions](https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_multiphysics_couplings.14.03.html)
- [Coupled elasto-acoustic displacement/potential formulation](https://arxiv.org/abs/1907.05405)
- [Gridap transient FEM framework](https://gridap.github.io/Tutorials/stable/pages/t017_transient_linear/)

## Loading and scripting

FEM is an optional package extension. Load `MetaNDT`, `Gridap`, `GridapGmsh`
and `Gmsh`; the detector-only library remains lightweight. Use
`simulate(CylinderExperiment())` or `simulate_comparison(spec)`.
`save_experiment` / `load_experiment` round-trip a versioned TOML description.
`sensor_csv` exports time in seconds and pressure in Pa. Results contain plain
mesh/snapshot arrays, so plotting and Pluto are not dependencies of the solver.
The solver snapshots the sensor list; changing the original input afterwards
does not rewrite a calculated result or a sweep entry. Meshing refuses an
already active caller-owned Gmsh session rather than modifying or closing it.
Finish that session first, or run the calculation in a separate process.
