# Local MVP verification — 2026-10-14

Environment: Julia 1.10.12 (LTS), Gridap 0.20.10, GridapGmsh 0.7.4,
Gmsh.jl 0.3.1 on macOS / Apple Silicon.

Command, from sibling `meta-demo`:

```
julia --startup-file=no --project=demo ../meta-ndt/test/cylinder_fem.jl
```

42 assertions passed: experiment serialization/validation (7), bidirectional
coupling and energy checks (23), frame-independent metrics / sweep provenance
(4), preservation of a caller-owned Gmsh session (4), zero forcing /
transparent-fluid limit (3), mesh+time refinement (1).

Test preset: 10 mm cylinder, 60 kHz two-cycle Hann pulse, 160 × 100 mm domain,
2 mm mesh size, 0.25 μs time step, 120 μs duration, three sensors.

| Check | Observed result |
|---|---:|
| Relative interface skew residual | 0 |
| Maximum relative work/energy imbalance, elastic | 1.24 × 10⁻¹⁵ |
| Maximum displacement over all time steps, elastic | 0.10494 nm |
| Near-fluid cylinder vs empty water, relative trace L2 difference | 1.33% |
| Refinement to 1.5 mm and 0.125 μs, relative trace L2 difference | 6.30% |

The trace comparisons include all three pressure sensors and all recorded
times. The transparent-fluid case changes the solid material, not the solver:
ρs=ρf, ν=0.4999, E=3ρf cf²(1−2ν).

Independent code/physics review found and fixed three regressions: mutable
sensor lists could rewrite result/cache provenance; the displacement peak
depended on retained frames; meshing could close a caller-owned Gmsh session.
Targeted tests with the actual solver/cache passed 10/10 assertions in addition
to the suite above. The interface-normal orientation was also checked by its
outward flux integral. The existing demo's 31 tests and detector-only library's
18 tests passed separately.

Cold-run assembly includes Julia compilation: about 16–17 s for each newly
compiled scalar/coupled branch in this run. A warm rigid case took about 0.49 s
including meshing/assembly and time stepping; its time loop was 0.16 s.
These are local observations, not performance promises for other hardware.

This validates passive exchange, signs, integration and a basic refinement
behavior for the stated preset. It does not validate every slider setting,
an exact nonreflecting boundary, finite 3D cylinders or hardware experiments.
One combined mesh/time refinement is a stability smoke test, not an established
convergence rate or an accuracy bound. Sixteen endpoint combinations of the
GUI's radius/frequency/source controls passed specification validation; they
were not all solved or independently accuracy-validated.
