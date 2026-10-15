using Test, LinearAlgebra, MetaNDT, Gridap, GridapGmsh, Gmsh

@testset "experiment contract and serialization" begin
    s=CylinderExperiment()
    @test validate_experiment(s) === s
    @test_throws ArgumentError validate_experiment(remake_experiment(s; radius=0.04))
    @test_throws ArgumentError remake_experiment(s; sensors=[(0.0,0.0)])
    @test_throws ArgumentError remake_experiment(s; frequency=Inf)
    copy_of_s=remake_experiment(s)
    @test copy_of_s.sensors !== s.sensors
    copy_of_s.sensors[1]=(0.060,0.0)
    @test s.sensors[1] == (0.040,0.0)
    mktempdir() do dir
        path=joinpath(dir,"experiment.toml")
        save_experiment(path,s)
        @test experiment_toml(load_experiment(path)) == experiment_toml(s)
    end
end

s=CylinderExperiment(duration=120e-6,frame_count=31)
results=simulate_comparison(s)
empty_water,rigid,elastic=results
@testset "real bidirectional coupling" begin
    original_sensors=copy(s.sensors)
    s.sensors[1]=(0.060,0.0)
    @test all(r.experiment.sensors == original_sensors for r in results)
    s.sensors[1]=original_sensors[1]
    @test elastic.diagnostics["interface_norm"] > 0
    @test elastic.diagnostics["interface_skew_relative"] < 1e-12
    @test elastic.diagnostics["peak_displacement_m"] > 1e-12
    @test maximum(elastic.solid_energy) > 0
    @test norm(elastic.sensor_pressure-rigid.sensor_pressure)/norm(rigid.sensor_pressure) > 0.005
    for r in results
        @test all(isfinite,r.pressure)
        @test r.diagnostics["energy_balance_relative"] < 1e-7
        @test minimum(r.energy) > -1e-12
        @test minimum(diff(r.absorbed_energy)) > -1e-12
        @test size(r.pressure,2) == length(r.frame_times)
    end
    @test isempty(rigid.displacement_x)
    @test isempty(empty_water.solid_nodes)
end

@testset "frame-independent displacement metric" begin
    sweep_spec=remake_experiment(s;frame_count=2)
    entry=only(sweep_experiments([sweep_spec]))
    sparse_frames=entry.result
    @test sparse_frames.sensor_pressure == elastic.sensor_pressure
    @test sparse_frames.diagnostics["peak_displacement_m"] == elastic.diagnostics["peak_displacement_m"]
    @test entry.experiment === sparse_frames.experiment
    sweep_spec.sensors[1]=(0.060,0.0)
    @test entry.experiment.sensors == s.sensors
end

@testset "preserve caller-owned Gmsh session" begin
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal",0)
        gmsh.model.add("caller_owned_geometry")
        disk=gmsh.model.occ.addDisk(0,0,0,0.01,0.01)
        gmsh.model.occ.synchronize()
        @test_throws ArgumentError simulate(s)
        @test gmsh.isInitialized() == 1
        @test gmsh.model.getCurrent() == "caller_owned_geometry"
        @test (2,disk) in gmsh.model.getEntities(2)
    finally
        gmsh.finalize()
    end
end

@testset "zero forcing and transparent-fluid limit" begin
    zero=simulate(remake_experiment(s; source_strength=0.0,frame_count=3))
    @test iszero(norm(zero.sensor_pressure))
    @test iszero(norm(zero.energy))
    K=s.water_density*s.water_speed^2
    ν=0.4999
    fluidlike=simulate(remake_experiment(s; solid_density=s.water_density,
        young_modulus=3K*(1-2ν),poisson_ratio=ν))
    error=norm(fluidlike.sensor_pressure-empty_water.sensor_pressure)/norm(empty_water.sensor_pressure)
    println("transparent-cylinder relative trace error = ",error)
    @test error < 0.12
end

@testset "mesh and time refinement" begin
    refined=simulate(remake_experiment(s;mesh_size=0.0015,time_step=0.125e-6))
    coarse=elastic.sensor_pressure
    fine=refined.sensor_pressure[:,1:2:end]
    error=norm(coarse-fine)/norm(fine)
    println("mesh+time refinement relative trace difference = ",error)
    @test error < 0.12
end

for r in results
    println(r.experiment.scenario," ",r.diagnostics)
end
