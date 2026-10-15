using Test, LinearAlgebra, TOML, MetaNDT, Gridap, GridapGmsh, Gmsh

@testset "array geometry and schema" begin
    centers=cylinder_lattice(2,3,0.018)
    @test length(centers) == 6
    @test sum(first,centers) ≈ 0 atol=1e-16
    @test sum(last,centers) ≈ 0 atol=1e-16
    @test cylinder_lattice(1,1,0.02;center=(0.01,-0.01)) == [(0.01,-0.01)]
    @test_throws ArgumentError cylinder_lattice(0,2,0.02)
    @test_throws ArgumentError cylinder_lattice(1,2,NaN)
    s=CylinderExperiment(radius=.004,centers=centers)
    @test validate_experiment(s) === s
    @test_throws ArgumentError remake_experiment(s;centers=NTuple{2,Float64}[])
    @test_throws ArgumentError remake_experiment(s;centers=[(0.0,0.0),(0.005,0.0)])
    @test_throws ArgumentError remake_experiment(s;centers=[(0.079,0.0)])
    @test_throws ArgumentError remake_experiment(s;centers=[(NaN,0.0)])
    @test_throws ArgumentError remake_experiment(s;centers=[s.source])
    @test_throws ArgumentError remake_experiment(s;centers=[s.sensors[1]])
    snapshot=remake_experiment(s)
    @test snapshot.centers !== s.centers
    mktempdir() do dir
        path=joinpath(dir,"array.toml")
        save_experiment(path,s)
        @test experiment_toml(load_experiment(path)) == experiment_toml(s)
        legacy=TOML.parse(experiment_toml(CylinderExperiment()))
        legacy["schema_version"]=1; delete!(legacy,"centers")
        open(path,"w") do io; TOML.print(io,legacy); end
        @test load_experiment(path).centers == [(0.0,0.0)]
    end
    @test MetaNDT.pressure_integral([0.,1.,2.],[1.,2.,3.],(.25,1.75)) ≈ (2.75^3-1.25^3)/3
    @test MetaNDT.pressure_integral([0.,1.],[2.,2.],(.25,.75)) ≈ 2.0
    @test_throws ArgumentError frequency_sweep(s,[60e3,40e3])
    @test_throws ArgumentError frequency_sweep(s,[Inf])
    @test_throws ArgumentError frequency_sweep(s,[60e3];sensor=0)
    @test_throws ArgumentError frequency_sweep(s,[60e3];time_window=(0.0,1.0))
    @test_throws ArgumentError frequency_sweep(remake_experiment(s;source_strength=0.0),[60e3])
end

@testset "translated cylinder is not an outer absorber" begin
    mktempdir() do dir
        path=joinpath(dir,"translated.msh")
        moved=CylinderExperiment(scenario=:rigid,radius=.004,centers=[(.068,0.0)],
            duration=80e-6,frame_count=2)
        simulate(moved;mesh_path=path)
        gmsh.initialize()
        try
            gmsh.option.setNumber("General.Terminal",0)
            gmsh.open(path)
            group=only([tag for (_,tag) in gmsh.model.getPhysicalGroups(1)
                if gmsh.model.getPhysicalName(1,tag) == "outer"])
            outer=gmsh.model.getEntitiesForPhysicalGroup(1,group)
            @test length(outer) == 4
            # Reloaded .msh curves are discrete entities, not OCC CAD lines.
            @test all(outer) do tag
                _,coordinates,_=gmsh.model.mesh.getNodes(1,tag,true,false)
                xyz=reshape(coordinates,3,:)
                all(axes(xyz,2)) do i
                    isapprox(abs(xyz[1,i]),moved.width/2;atol=1e-10) ||
                        isapprox(abs(xyz[2,i]),moved.height/2;atol=1e-10)
                end
            end
        finally
            gmsh.finalize()
        end
    end
end

s=CylinderExperiment(radius=.004,centers=cylinder_lattice(2,2,.018),duration=120e-6,frame_count=7)
empty_water,rigid,elastic=simulate_comparison(s)
@testset "all inclusions are coupled" begin
    @test elastic.diagnostics["cylinder_count"] == 4
    @test rigid.diagnostics["cylinder_count"] == 4
    @test empty_water.diagnostics["cylinder_count"] == 0
    @test elastic.diagnostics["interface_skew_relative"] < 1e-12
    @test elastic.diagnostics["energy_balance_relative"] < 1e-7
    @test norm(elastic.sensor_pressure-rigid.sensor_pressure)/norm(rigid.sensor_pressure) > .002
    for center in s.centers
        ids=findall(i -> hypot(elastic.solid_nodes[1,i]-center[1],elastic.solid_nodes[2,i]-center[2]) <= 1.001s.radius,
            axes(elastic.solid_nodes,2))
        @test !isempty(ids)
        @test maximum(hypot.(elastic.displacement_x[ids,:],elastic.displacement_y[ids,:])) > 1e-12
    end
    original_centers=copy(s.centers)
    s.centers[1]=(-0.025,-0.025)
    @test elastic.experiment.centers == original_centers
    s.centers[1]=original_centers[1]
end

@testset "array near-fluid limit and refinement" begin
    K=s.water_density*s.water_speed^2; ν=.4999
    fluidlike=simulate(remake_experiment(s;solid_density=s.water_density,
        young_modulus=3K*(1-2ν),poisson_ratio=ν))
    transparent_error=norm(fluidlike.sensor_pressure-empty_water.sensor_pressure)/norm(empty_water.sensor_pressure)
    println("array transparent-cylinder relative trace error = ",transparent_error)
    @test transparent_error < .12
    refined=simulate(remake_experiment(s;mesh_size=.0015,time_step=.125e-6))
    refinement_error=norm(elastic.sensor_pressure-refined.sensor_pressure[:,1:2:end])/norm(refined.sensor_pressure[:,1:2:end])
    println("array refinement relative trace difference = ",refinement_error)
    @test refinement_error < .15
end

@testset "paired frequency sweep" begin
    scan=frequency_sweep(s,[40e3,60e3,80e3];sensor=1,time_window=(.25e-6,119.75e-6))
    @test scan.frequencies == [40e3,60e3,80e3]
    @test scan.ratios[1,:] == ones(3)
    @test all(isfinite,scan.ratios)
    @test all(>(0),scan.reference_integrals)
    for (j,cases) in enumerate(scan.comparisons), r in cases
        @test r.experiment.frequency == scan.frequencies[j]
        @test r.experiment.centers == s.centers
        @test r.experiment.mesh_size == scan.comparisons[1][1].experiment.mesh_size
        @test length(r.frame_times) == 2
        @test r.diagnostics["energy_balance_relative"] < 1e-7
    end
    point=scan.comparisons[2][3]
    replay=simulate(remake_experiment(point.experiment;frame_count=7))
    @test replay.sensor_pressure == point.sensor_pressure
    double_source=frequency_sweep(remake_experiment(point.experiment;source_strength=2s.source_strength),[60e3];
        sensor=scan.sensor,time_window=scan.time_window)
    @test double_source.ratios[:,1] ≈ scan.ratios[:,2] rtol=1e-12
    csv=split(chomp(frequency_sweep_csv(scan)),'\n')
    @test length(csv) == 4
    @test length(split(csv[2],',')) == 11
    println("array frequency sweep rigid/elastic ratios = ",scan.ratios[2:3,:])
    println("array peak displacement m = ",elastic.diagnostics["peak_displacement_m"])
end
