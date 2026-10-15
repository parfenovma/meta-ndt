using Test, MetaNDT

@testset "P1 point displacement without loading a FEM backend" begin
    # An affine vector field on two disconnected triangles: interpolation must
    # be exact and the space between them must not become a zero-valued sensor.
    xy=[0.0 1 0 3 4 3; 0 0 1 0 0 1]
    triangles=[1 4;2 5;3 6]
    ux=ComplexF64[1+2x-3y+im*(2-x+y) for (x,y) in eachcol(xy)]
    uy=ComplexF64[-2+x+4y+im*(1+2x-y) for (x,y) in eachcol(xy)]
    r=SolidLensResult(SolidLensExperiment(),xy,triangles,ux,uy,Float64[],Float64[],Float64[],Float64[],Float64[],Dict{String,Float64}())
    for p in ((.2,.3),(0.,0.),(.5,.5),(3.2,.3))
        u=solid_displacement(r,p); x,y=p
        @test u[1]≈1+2x-3y+im*(2-x+y)
        @test u[2]≈-2+x+4y+im*(1+2x-y)
    end
    for p in ((2.,.2),(-.1,0.),(0.,1.1),(NaN,0.),(0.,Inf),(1.,))
        @test_throws ArgumentError solid_displacement(r,p)
    end
    @test Base.get_extension(MetaNDT,:MetaNDTGridapExt)===nothing
end

@testset "Harmonic trace convention, units and validation" begin
    f=2000.; t=collect((0,1/(4f),1/(2f),3/(4f),1/f))
    @test harmonic_trace(2+3im,f,t)≈[2,3,-2,-3,2]
    @test harmonic_trace(2e-9+3e-9im,f,t)≈[2,3,-2,-3,2]*1e-9
    @test harmonic_trace(0im,f,t)==zeros(5)
    @test harmonic_trace(1+2im,f,Float64[])==Float64[]
    for (a,freq,ts) in ((NaN+im,f,t),(1,-1,t),(1,0,t),(1,Inf,t),(1,f,[Inf]))
        @test_throws ArgumentError harmonic_trace(a,freq,ts)
    end
    water=simulate(CoatedCylinderExperiment(scenario=:empty))
    p=(.012,.004); A=coated_pressure(water,p); ω=2π*water.experiment.frequency
    times=collect(range(0,2/water.experiment.frequency;length=241))
    expected=water.experiment.incident_pressure*cos.(ω*p[1]/water.experiment.water_speed.-ω*times)
    @test harmonic_trace(A,water.experiment.frequency,times)≈expected atol=1e-13
end
