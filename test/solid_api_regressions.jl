using Test, LinearAlgebra, MetaNDT, Gridap, GridapGmsh, Gmsh
BLAS.set_num_threads(1)
@testset "solid review regressions" begin
    b=SinusoidalBrick(length=.020,height=.008,amplitude=.003,end_margin=.0006)
    s=SolidLensExperiment(bricks=fill(b,13))
    d=LensDesign(s,[b],[1+0im],ones(Int,13),zeros(13))
    uniform=solid_preset(d,:uniform)
    @test all(c->c.length==b.length && c.height==b.height && c.end_margin==b.end_margin,uniform.bricks)
    @test first(uniform.bricks).amplitude/first(uniform.bricks).height ≈ 5/14
    pressure_design=LensDesign(remake_experiment(s;polarization=:pressure),d.candidates,
        d.transmission,d.selected,d.phase_error)
    @test solid_preset(pressure_design,:focused).polarization == :pressure
    small=remake_experiment(SolidLensExperiment();width=.04,source_x=-.014,
        source_width=.0005,focus_distance=.004,mesh_size=.0007,near_mesh_size=.0005)
    r=simulate(small)
    @test all(x->first(small.bricks).length/2 < x < small.width/2,r.scan_x)
    @test all(isfinite,r.scan_peak)
    @test r.diagnostics["power_balance_relative"] < 1e-8
end
