using Test, LinearAlgebra, TOML, MetaNDT, Gridap, GridapGmsh, Gmsh

polygon_area(p)=abs(sum(p[i][1]*p[mod1(i+1,length(p))][2]-
    p[mod1(i+1,length(p))][1]*p[i][2] for i in eachindex(p)))/2

@testset "correct shifted STG geometry" begin
    b=SinusoidalBrick()
    span=MetaNDT.stg_active_span(b); shift=span/(2b.periods)
    @test span ≈ .01296
    @test shift ≈ .00324
    @test b.height-b.amplitude ≈ .0045
    for x in range(b.end_margin+shift,b.end_margin+span;length=101)
        bottom=MetaNDT.stg_bottom_depth(b,x)
        top=b.height-MetaNDT.stg_top_depth(b,x)
        @test top-bottom ≈ .0045 atol=1e-14
    end
    p=MetaNDT.brick_polygon(b)
    @test polygon_area(p) ≈ 86.6e-6 rtol=1e-12
    @test all(i -> hypot((p[i].-p[mod1(i+1,length(p))])...) > 1e-8,eachindex(p))
    @test extrema(first.(p)) == (-.0085,.0085)
    @test extrema(last.(p)) == (-.0035,.0035)
    @test MetaNDT.polygon_signed_distance((0.0,0.0),p) < 0
    @test MetaNDT.polygon_signed_distance((.020,0.0),p) > 0
    @test_throws ArgumentError MetaNDT.brick_polygon(SinusoidalBrick(amplitude=.007))
    @test_throws ArgumentError MetaNDT.brick_polygon(SinusoidalBrick(samples=81))
    zero_margin=MetaNDT.brick_polygon(SinusoidalBrick(end_margin=0.0))
    @test all(i -> hypot((zero_margin[i].-zero_margin[mod1(i+1,length(zero_margin))])...) > 1e-8,eachindex(zero_margin))
end

@testset "frozen presets and schema" begin
    examples=lab_examples()
    @test length(examples) == 6
    @test length(unique(e.id for e in examples)) == 6
    for e in examples, f in (40e3,100e3,160e3)
        s=example_experiment(e.id;frequency=f)
        @test validate_experiment(s) === s
        mktempdir() do dir
            path=joinpath(dir,"example.toml")
            save_experiment(path,s)
            @test experiment_toml(load_experiment(path)) == experiment_toml(s)
        end
    end
    s=example_experiment(:stg_single)
    ref=straight_reference(s)
    @test specimen_area(ref) ≈ specimen_area(s) rtol=1e-12
    @test ref.shape.length == s.shape.length
    @test ref.centers == s.centers
    @test ref.centers !== s.centers
    @test_throws ArgumentError remake_experiment(s;sensors=[(0.0,0.0)])
    @test_throws ArgumentError remake_experiment(s;mesh_size=.002)
    @test_throws ArgumentError remake_experiment(s;centers=[(0.0,0.0),(.010,0.0)])
    @test_throws ArgumentError straight_reference(CylinderExperiment())
    @test_throws ArgumentError frequency_sweep(CylinderExperiment(),[60e3];include_straight_reference=true)
    @test length(example_experiment(:missing_cell).centers) == 8
    @test length(example_reference(:missing_cell,example_experiment(:missing_cell)).centers) == 9
    mktempdir() do dir
        path=joinpath(dir,"legacy.toml")
        old=TOML.parse(experiment_toml(CylinderExperiment()))
        old["schema_version"]=2; delete!(old,"shape")
        open(path,"w") do io; TOML.print(io,old); end
        @test isnothing(load_experiment(path).shape)
    end
end

@testset "actual FEM for six examples" begin
    for e in lab_examples()
        s=remake_experiment(example_experiment(e.id);frame_count=7)
        r=simulate(s)
        @test r.diagnostics["inclusion_count"] == length(s.centers)
        @test r.diagnostics["energy_balance_relative"] < 1e-7
        @test r.diagnostics["interface_skew_relative"] < 1e-12
        @test r.diagnostics["peak_displacement_m"] > 1e-12
        @test all(isfinite,r.sensor_pressure)
        @test !isempty(r.solid_triangles)
        println("gallery ",e.id," ",r.diagnostics)
    end
end

@testset "STG coupled limit, refinement and matched control sweep" begin
    s=remake_experiment(example_experiment(:stg_single);frame_count=7)
    water,rigid,elastic=simulate_comparison(s)
    @test norm(elastic.sensor_pressure-rigid.sensor_pressure)/norm(rigid.sensor_pressure) > .002
    zero=simulate(remake_experiment(s;source_strength=0.0,frame_count=2))
    @test iszero(norm(zero.energy))
    K=s.water_density*s.water_speed^2; ν=.4999
    fluidlike=simulate(remake_experiment(s;solid_density=s.water_density,
        young_modulus=3K*(1-2ν),poisson_ratio=ν))
    transparent_error=norm(fluidlike.sensor_pressure-water.sensor_pressure)/norm(water.sensor_pressure)
    println("STG transparent trace error = ",transparent_error)
    @test transparent_error < .12
    refined=simulate(remake_experiment(s;mesh_size=.00075,time_step=.125e-6))
    refinement_error=norm(elastic.sensor_pressure-refined.sensor_pressure[:,1:2:end])/norm(refined.sensor_pressure[:,1:2:end])
    println("STG refinement trace difference = ",refinement_error)
    @test refinement_error < .15
    scan=frequency_sweep(s,[60e3,100e3,140e3];include_straight_reference=true)
    @test size(scan.ratios) == (4,3)
    @test scan.ratios[1,:] == ones(3)
    @test all(isfinite,scan.ratios)
    for cases in scan.comparisons
        @test cases[3].experiment.shape isa SinusoidalBrick
        @test cases[4].experiment.shape isa RectangularBrick
        @test specimen_area(cases[3].experiment) ≈ specimen_area(cases[4].experiment) rtol=1e-12
        @test cases[3].times == cases[4].times
        @test cases[3].experiment.source == cases[4].experiment.source
        @test cases[3].experiment.mesh_size == cases[4].experiment.mesh_size
        @test cases[4].diagnostics["energy_balance_relative"] < 1e-7
    end
    @test length(split(split(chomp(frequency_sweep_csv(scan)),'\n')[2],',')) == 12
    point=scan.comparisons[2][3]
    replay=simulate(point.experiment)
    @test replay.sensor_pressure == point.sensor_pressure
    println("STG sweep rows empty/rigid/STG/straight = ",scan.ratios)
end
